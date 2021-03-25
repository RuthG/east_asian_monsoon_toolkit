"""wang_fan_1999.py: Calculate the vorticity index of Wang and Fang 1999: 
https://doi.org/10.1175/1520-0477(1999)080%3C0629:COSASM%3E2.0.CO;2

U850 in (5°–15°N, 90°–130°E) minus U850 in (22.5°–32.5°N, 110°–140°E)

This was demonstrated by Wang et al. 2008 (https://doi.org/10.1175/2008JCLI2183.1) 
to well-approximate PC1 of the East Asian Summer Monsoon and be a useful indicator 
of monsoon strength.

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr

def wang_fan(u, areacell=None, latdim='lat', londim='lon', pdim='plev', punits='Pa'):
    # Inputs:
    # u: lat-lon DataArray of zonal wind speed. Input either only 850-hPa level, or specify pressure dimension to search for this level over
    # areacell: Grid of cell areas for spatial averaging
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    # pdim: Name of pressure dimension, default plev
    # punits: Pressure units, default Pa
    
    # Returns:
    # Time series of Wang and Fan 1999 index, where time axis matches that of input
    
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
            
    
    # Identify grid points in (5°–15°N, 90°–130°E)
    lons1 = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= 90.  and u.lon[i] <= 130.]
    lats1 = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 5.   and u.lat[i] <= 15.]
    
    # Identify grid points in (22.5°–32.5°N, 110°–140°E)
    lons2 = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= 110.  and u.lon[i] <= 140.]
    lats2 = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 22.5  and u.lat[i] <= 32.5]
    
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
            u2 = u_wt.sel(lon=lons2, lat=lats2, method='nearest').sum(('lat','lon')) / areacell.sel(lat=lats2, lon=lons2, method='nearest').sum(('lat','lon')) 
        else:
            print('Warning, cell area dimensions do not match those of u, defaulting to using cosine weighted averaging')
            areacell = None
            
    if areacell is None:
        coslat = np.cos(u.lat * np.pi/180.)
        u_wt = u * coslat # In absence of cell area grid weight average by latitude before averaging
        u1 = u_wt.sel(lon=lons1, lat=lats1, method='nearest').sum('lat').mean('lon') / coslat.sel(lat=lats1, method='nearest').sum('lat')
        u2 = u_wt.sel(lon=lons2, lat=lats2, method='nearest').sum('lat').mean('lon') / coslat.sel(lat=lats2, method='nearest').sum('lat')
    
    return u1-u2 #Return shear
   
   
