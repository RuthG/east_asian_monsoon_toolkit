"""yim_etal_2014_vorticity.py: Calculate the vorticity indices of Yim at al. 2014
https://doi.org/10.1007/s00382-013-1956-9

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr


regions = {
 'ISM': [[[5,15],[40,80]], [[25,35],[70,90]]],
 'WNPSM': [[[5,15],[100,130]], [[20,35],[110,140]]],
 'NASM': [[[5,15],[230,260]], [[20,30],[250,280]]],
 'NAFSM': [[[0,15],[300,350]], [[None,None],[None,None]]],
 'SASM':[[[-20,-5],[290,320]], [[-35,-20],[290,320]]],
 'SAFSM':[[[-15,-5],[20,50]], [[-30,-20],[30,55]]],
 'AUSSM':[[[-15,0],[90,130]], [[-30,-20],[100,140]]]
  }


def yim_vort(u, areacell=None, latdim='lat', londim='lon', pdim='plev', punits='Pa', region='ISM'):
    # Inputs:
    # u: lat-lon(+time/pressure) DataArray of zonal wind speed. Input either only 850-hPa level, or specify pressure dimension to search for this level over
    # areacell: Grid of cell areas for spatial averaging
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    # pdim: Name of pressure dimension, default plev
    # punits: Pressure units, default Pa
    # region: Key for 'regions' dictionary
    
    # Returns:
    # Time series of Yim et al. 2014 index, where time axis matches that of input
    
    # Rename coords for ease if needed
    u = u.rename({latdim:'lat', londim:'lon'})
    
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
            
    [[[lat1a,lat1b],[lon1a,lon1b]], [[lat2a,lat2b],[lon2a,lon2b]]] = regions[region]
    
    if u.lon.min() < -10.:
        print('check')
        lon1a = lon1a - 360.; lon1b = lon1b - 360.; lon2a = lon2a - 360.; lon2b = lon2b - 360.
    
    # Identify grid points in region
    lons1 = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= lon1a  and u.lon[i] <= lon1b]
    lats1 = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= lat1a   and u.lat[i] <= lat1b]
    
    if region is not 'NAFSM':
        lons2 = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= lon2a  and u.lon[i] <= lon2b]
        lats2 = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= lat2a  and u.lat[i] <= lat2b]
    
    # Take weighted means of u over above boxes
    if areacell is not None:
        areacell = areacell.rename({latdim:'lat', londim:'lon'})
        if (len(areacell.lat) == len(u.lat)) and (len(areacell.lon) == len(u.lon)):
            latsame = (np.round(areacell.lat.values - u.lat.values,2) ==0.).all()
            lonsame = (np.round(areacell.lon.values - u.lon.values,2) ==0.).all()
        else:
            latsame = False;  lonsame = False
        if latsame and lonsame:
            u_wt = u * areacell # If cell area grid provided use this to area weight values before averaging
            u1 = u_wt.sel(lon=lons1, lat=lats1, method='nearest').sum(('lat','lon')) / areacell.sel(lat=lats1, lon=lons1, method='nearest').sum(('lat','lon')) 
            
            if region is not 'NAFSM':
                u2 = u_wt.sel(lon=lons2, lat=lats2, method='nearest').sum(('lat','lon')) / areacell.sel(lat=lats2, lon=lons2, method='nearest').sum(('lat','lon')) 
            else:
                u2=0.
        else:
            print('Warning, cell area dimensions do not match those of u, defaulting to using cosine weighted averaging')
            areacell = None
            
    if areacell is None:
        coslat = np.cos(u.lat * np.pi/180.)
        u_wt = u * coslat # In absence of cell area grid weight average by latitude before averaging
        u1 = u_wt.sel(lon=lons1, lat=lats1, method='nearest').sum('lat').mean('lon') / coslat.sel(lat=lats1, method='nearest').sum('lat')
        
        if region is not 'NAFSM':
            u2 = u_wt.sel(lon=lons2, lat=lats2, method='nearest').sum('lat').mean('lon') / coslat.sel(lat=lats2, method='nearest').sum('lat')
        else:
            u2=0.
            
    return u1-u2 #Return shear
   
   
