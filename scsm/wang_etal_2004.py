"""wang_etal_2004.py: 
Calculate the South China Sea monsoon onset index of Wang et al. 2004:
https://doi.org/10.1175/2932.1

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr

def scsm_onset(u, areacell=None, pentaddim='pentad', pdim='plev', punits='Pa', latdim='lat', londim='lon'):
    # Inputs:
    # u: 1 year long lat-lon-pentad DataArray of zonal wind speed. Input either only 850-hPa level, or specify pressure dimension
    # areacell: Grid of cell areas for spatial averaging
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    # pdim: Name of pressure dimension, default plev
    # punits: Pressure units, default Pa
    
    # Returns:
    # Area average of u over relevant box
    # Onset pentad of SCSM for that year as per Wang et al. 2004 index, or None if none found
    
    # Rename coords for ease if needed
    u = u.rename({pentaddim: 'pentad', latdim: 'lat', londim: 'lon'})
    
    if pdim in u.dims:
        if punits == 'Pa':
            lev=85000.
        elif punits == 'hPa':
            lev=850.
        else:
            print('Error, punits not recognised')
            return
        u = u.rename({pdim:'plev'})
        try:
            u = u.sel(plev=lev)
        except:
            print('Warning, 850hPa level not found, looking for nearest level')
            u = u.sel(plev=lev, method='nearest')
            print('Nearest level found: ' + str(int(u.plev.values)))
    
    
    # Get specified lats and lons, select pentad range to look at
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 5. and u.lat[i] <= 15.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= 110. and u.lon[i] <= 120.]
    pentads = [u.pentad[i].values for i in range(len(u.pentad)) if u.pentad[i] >=  24]
    
    # Take weighted means of u over above box
    if areacell is not None:
        areacell = areacell.rename({latdim:'lat', londim:'lon'})
        if (len(areacell.lat) == len(u.lat)) and (len(areacell.lon) == len(u.lon)):
            latsame = (np.round(areacell.lat.values - u.lat.values,2) ==0.).all()
            lonsame = (np.round(areacell.lon.values - u.lon.values,2) ==0.).all()
        else:
            latsame = False;  lonsame = False
        if latsame and lonsame:
            u_wt = u * areacell # If cell area grid provided use this to area weight values before averaging
            u_mean = u_wt.sel(lon=lons, lat=lats, method='nearest').sum(('lat','lon')) / areacell.sel(lat=lats, lon=lons, method='nearest').sum(('lat','lon')) 
        else:
            print('Warning, cell area dimensions do not match those of u, defaulting to using cosine weighted averaging')
            areacell = None
            
    if areacell is None:
        coslat = np.cos(u.lat * np.pi/180.)
        u_wt = u * coslat # In absence of cell area grid weight average by latitude before averaging
        u_mean = u_wt.sel(lon=lons, lat=lats, method='nearest').sum('lat').mean('lon') / coslat.sel(lat=lats, method='nearest').sum('lat')
    
    
    for i in range(len(pentads)):  # For each pentad after pentad 24
        u_pos_mean = (u_mean.sel(pentad=pentads[i+1:i+5]).where(u_mean.sel(pentad=pentads[i+1:i+5]).values > 0.)).mean('pentad')
        if ((u_mean.sel(pentad=pentads[i]).values) > 0.  # Check if u_mean is greater than zero
             and (u_mean.sel(pentad=pentads[i:i+4]).sum('pentad') > 1.) # and if the u_mean over that pentad and next 3 is greater than 1
             and (np.sum(u_mean.sel(pentad=pentads[i:i+4]).values > 0.) >= 2.)): # and if u_mean is greater than zero in at least 3 out of 4 pentads
            onset_pentad = pentads[i]  # If all that is true, that's your onset pentad
            return u_mean, onset_pentad # Return it, and u
    
    return u_mean, None # Return u and an empty value if no onset found




    