"""li_zhang_2009.py: Calculate the wind-rotation based monsoon onset and withdrawal indices of Li and Zhang 2009.
https://doi.org/10.1007/s00382-008-0465-8

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr

def rename_coords(var_in, pdim='plev', punits='Pa', latdim='lat', londim='lon'):
    
    var_in = var_in.rename({latdim: 'lat', londim: 'lon'})

    if pdim in var_in.dims:
        if punits == 'Pa':
            lev=85000.
        elif punits == 'hPa':
            lev=850.
        else:
            print('Error, punits not recognised')
            return
        var_in = var_in.rename({pdim:'plev'})
        try:
            var_in = var_in.sel(plev=lev)
        except:
            print('Warning, 850hPa level not found, looking for nearest level')
            var_in = var_in.sel(plev=lev, method='nearest')
            print('Nearest level found: ' + str(int(var_in.plev.values)))
    
    return var_in


def abs_angle(u, v, u_rel, v_rel):
    u_mag = np.sqrt(u**2. + v**2.)
    u_rel_mag = np.sqrt(u_rel**2. + v_rel**2.)
    return np.arccos( (u_rel*u + v_rel*v) / (u_mag * u_rel_mag) )


def smooth_data(var_in, timedim='time'):
    # Smooth data
    var_mean = var_in.mean(timedim) # save mean
    var_in = var_in - var_mean # subtract
    axisno = var_in.get_axis_num(timedim) # Find number of pentad dimension
    var_fft = np.fft.fft(var_in.values, axis=axisno) # Take fourer transform along pentad dimension - use fft not rfft to conserve axis length
    # Discard all but 1st 12 harmonics
    if axisno==0:
        var_fft[12:-12,:,:] = 0
    elif axisno==1:
        var_fft[:,12:-12,:] = 0
    elif axisno==2:
        var_fft[:,:,12:-12] = 0
    var_smooth = np.fft.ifft(var_fft, axis=axisno) # Take inverse fourier transform to recover smoothed timeseries
    var_smooth = xr.DataArray(var_smooth.real, coords=var_in.coords, dims=var_in.dims) # Make dataarray
    var_smooth = var_smooth + var_mean # Add mean back on
    return var_smooth
        

def li_zhang(u, v, p=15, areacell=None, latdim='lat', londim='lon', pdim='plev', punits='Pa', smooth=False):
    # Inputs:
    # u, v: lat-lon-datetime (+pressure) DataArrays of zonal and meridional wind speed. Input either only 850-hPa level, or specify pressure dimension to search for this level over
    # areacell: Grid of cell areas for spatial averaging
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    # pdim: Name of pressure dimension, default plev
    # punits: Pressure units, default Pa
    # smooth: If true, smooth data using mean + 12 harmonics
    
    # Returns:
        
    u = rename_coords(u, latdim=latdim, londim=londim, pdim=pdim, punits=punits)
    v = rename_coords(v, latdim=latdim, londim=londim, pdim=pdim, punits=punits)
    
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 0. and u.lat[i] <= 40.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= 40. and u.lon[i] <= 180.]
    u = u.sel(lon=lons, lat=lats);     v = v.sel(lon=lons, lat=lats)
    
    if smooth:
        u = smooth_data(u); v = smooth_data(v)
    
    try:
        u_jan = u.sel(time=u.time.dt.month.isin([1])).mean('time') 
        v_jan = v.sel(time=u.time.dt.month.isin([1])).mean('time') 
        u_julaug = u.sel(time=u.time.dt.month.isin([7,8])).mean('time')   
        v_julaug = v.sel(time=u.time.dt.month.isin([7,8])).mean('time')
    except:
        print('Error: Time axis is not datetime indexable')
        return
    
    beta = abs_angle(u, v, u_jan, v_jan)
    beta_julaug = abs_angle(u_julaug, v_julaug, u_jan, v_jan)
    beta_onset = beta.sel(time=u.time.dt.month.isin(range(7))) 
    beta_withdrawal = beta.sel(time=u.time.dt.month.isin(range(9,13)))
    
    # Condition a, find first time where beta remains above half the july/aug climatological value for p days or more
    cond_a_onset = beta_onset > beta_julaug/2.
    
    a = np.zeros(cond_a_onset.shape)
    for i in range(len(cond_a_onset.time)-p):
        a[:,:,i] = xr.where(cond_a_onset.isel(time=range(i,i+p)).sum('time') == p, i, np.nan).values
    a = xr.DataArray(a, dims=cond_a_onset.dims, coords=cond_a_onset.coords) 
    a = a.isel(time=range(1,len(a.time)-p))
    a = a.min('time')
    
    # Condition b, find the time where the slope of beta changes fastest
    b = np.zeros(beta_onset.shape)
    for i in range(10,len(beta_onset.time)-10):
        b1 = beta_onset.isel(time=range(i-10,i)).differentiate('time').mean('time')
        b2 = beta_onset.isel(time=range(i+1,i+10)).differentiate('time').mean('time')
        a_mask = i >= a.values
        b[:,:,i] = (np.arctan(b1)-np.arctan(b2)) * a_mask
    
    b = xr.DataArray(b, dims=beta_onset.dims, coords=beta_onset.coords) 
    onset = np.round((b.argmax('time')+1) /5.)
    
    
    # Condition a, find first time where beta remains above half the july/aug climatological value for p days or more
    cond_a_withdrawal = beta_withdrawal > beta_julaug/2.
    #print(cond_a_withdrawal[0,0,:])
    a = np.zeros(cond_a_withdrawal.shape)
    for i in range(len(cond_a_withdrawal.time)-1, p-1, -1):
        a[:,:,i] = xr.where(cond_a_withdrawal.isel(time=range(i,i-p,-1)).sum('time') == p, i, np.nan).values
    a = xr.DataArray(a, dims=cond_a_withdrawal.dims, coords=cond_a_withdrawal.coords) 
    a = a.isel(time=range(p,len(a.time)))
    a = a.max('time')
    
    # Condition b, find the time where the slope of beta changes fastest
    b = np.zeros(beta_withdrawal.shape)
    for i in range(10,len(beta_withdrawal.time)-10):
        b1 = beta_withdrawal.isel(time=range(i-10,i)).differentiate('time').mean('time')
        b2 = beta_withdrawal.isel(time=range(i+1,i+10)).differentiate('time').mean('time')
        a_mask = i <= a.values
        b[:,:,i] = (np.arctan(b1)-np.arctan(b2)) * a_mask
    
    b = xr.DataArray(b, dims=beta_withdrawal.dims, coords=beta_withdrawal.coords) 
    offset = np.where(beta.time.values==b.time.min().values)[0][0]
    withdrawal = np.round((b.argmax('time')+1 + offset) /5.)
    
    
    return onset, withdrawal
   
