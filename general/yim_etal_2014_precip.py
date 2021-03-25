"""yim_etal_2014_precip.py: Calculate the precip region means of Yim at al. 2014
https://doi.org/10.1007/s00382-013-1956-9

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr


regions = {
 'IN': [[10,30],[70,105]],
 'WNP': [[12.5,22.5],[110,150]],
 'EA':  [[22.5,45],[110,135]],
 'NAM': [[7.5,22.5],[250,280]],
 'NAF': [[5,15],[330,30]],
 'SAM': [[-25,-5],[290,320]],
 'SAF': [[-25,-7.5],[25,70]],
 'AUS': [[-20,-5],[110,150]],
  }


def yim_precip(p, areacell=None, latdim='lat', londim='lon', pdim='plev', punits='Pa', region='IN'):
    # Inputs:
    # p: lat-lon(+time/pressure) DataArray of precipitation.
    # areacell: Grid of cell areas for spatial averaging
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    # region: Key for 'regions' dictionary
    
    # Returns:
    # Time series of Yim et al. 2014 precip mean, where time axis matches that of input
    
    # Rename coords for ease if needed
    p = p.rename({latdim:'lat', londim:'lon'})
    
      
    [[lata,latb],[lona,lonb]] = regions[region]
    
    if p.lon.min() < -10.:
        if region is not 'NAF':
            lona = lona - 360.; lonb = lonb - 360.
        else:
            lona = lona - 360.;
    
    # Identify grid points in region
    if lonb > lona:
        lons = [p.lon[i].values for i in range(len(p.lon)) if p.lon[i] >= lona  and p.lon[i] <= lonb]
    else:
        lons = [p.lon[i].values for i in range(len(p.lon)) if p.lon[i] >= lona or p.lon[i] <= lonb]
    lats = [p.lat[i].values for i in range(len(p.lat)) if p.lat[i] >= lata   and p.lat[i] <= latb]
    
    # Take weighted means of u over above boxes
    if areacell is not None:
        areacell = areacell.rename({latdim:'lat', londim:'lon'})
        if (len(areacell.lat) == len(p.lat)) and (len(areacell.lon) == len(p.lon)):
            latsame = (np.round(areacell.lat.values - p.lat.values,2) ==0.).all()
            lonsame = (np.round(areacell.lon.values - p.lon.values,2) ==0.).all()
        else:
            latsame = False;  lonsame = False
        if latsame and lonsame:
            p_wt = p * areacell # If cell area grid provided use this to area weight values before averaging
            p_mean = p_wt.sel(lon=lons, lat=lats, method='nearest').sum(('lat','lon')) / areacell.sel(lat=lats, lon=lons, method='nearest').sum(('lat','lon')) 
            
        else:
            print('Warning, cell area dimensions do not match those of p, defaulting to using cosine weighted averaging')
            areacell = None
            
    if areacell is None:
        coslat = np.cos(p.lat * np.pi/180.)
        p_wt = p * coslat # In absence of cell area grid weight average by latitude before averaging
        p_mean = p_wt.sel(lon=lons, lat=lats, method='nearest').sum('lat').mean('lon') / coslat.sel(lat=lats, method='nearest').sum('lat')

            
    return p_mean #Return area mean precip
   
