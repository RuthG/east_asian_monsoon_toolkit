"""li_etal_2018.py
Identify Meiyu front based on equivalent potential temperature criteria presented in Li et al. 2018:
https://doi.org/10.1007/s00382-017-3975-4

* NB Method has been coded up and initial tests done but should be checked on a hi-res reanalysis dataset. 
Resolution dependence of thresholds is likely to need adjusting.

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr

cp_air = 1004.6
L = 2.507e6
rdgas = 287.04
rvgas = 461.50 
tfreeze = 273.16
a = 6371.e3

def sat_vap_pres(t):
    # Calculate saturation vapor pressure 
    es = 610.78 * np.exp(-1.* L/rvgas * (1/t - 1/tfreeze))
    # e = rp/0.622 -> r = 0.622e/p
    rs = 0.622 * es / (t.plev)
    return es, rs
    
    
def equiv_pot_t_ams(t, q):
    # Calculate: theta_equiv = T(po/pd)^[Rd/(cpd + rtc)] * H^[-rvRv/(cpd + rtc)] * e^[Lvrv / ((cpd + rtc)T)]
        
    # Assume water only present in air in vapour form, so rt = rv = r:
    r = q/(1-q)
    c = 4217. # heat capacity of liquid water at 0 degrees
    denom = cp_air + r*c

    # H = q/qs = e/es  relative humidity
    es, rs = sat_vap_pres(t)
    rh = r/rs
    denom_s = cp_air + rs*c
    
    theta_equiv = t * (100000./t.plev)**(rdgas/denom) * rh**(-r * rvgas/denom) * np.exp(L * r/denom/t)
    
    theta_equiv_s = t * (100000./t.plev)**(rdgas/denom_s) * np.exp(L * rs/denom_s/t)
     
    return theta_equiv, theta_equiv_s
    


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
    

def li_meiyu(t, q, areacell=None, pdim='plev', punits='Pa', latdim='lat', londim='lon'):
    # Inputs:
    # t, q: lat-lon (+pressure, time) DataArrays of temperature and specific humidity
    # areacell: Grid of cell areas for spatial averaging
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    # pdim: Name of pressure dimension, default plev
    # punits: Pressure units, default Pa
    
    # Returns: 
    # mbf_lats: Meiyu-Baiu front latitudes for each input time
    # dtheta/dy: Meridional gradient of equivalent potential temperature for checking
    
    # Rename coords for ease if needed
    q = rename_coords(q, pdim=pdim, punits=punits, latdim=latdim, londim=londim)
    t = rename_coords(t, pdim=pdim, punits=punits, latdim=latdim, londim=londim)
    
    
    lonres = (t.lon[1]-t.lon[0]).values
    latres = (t.lat[1]-t.lat[0]).values
        
    # If resolution is lower than original study, relax the threshold on dthetady as band may be blurred
    if latres > 0.5:
        dthdy_threshold = 0.04/latres
        continuity_threshold = latres
    else:
        dthdy_threshold = 0.04
        continuity_threshold = 1.
        
    lons = [t.lon[i].values for i in range(len(t.lon)) if t.lon[i] >= 105. and t.lon[i] <= 145.] # 80 * 36 = 2880
    lats = [t.lat[i].values for i in range(len(t.lat)) if t.lat[i] >= 22.  and t.lat[i] <= 40.]
    
    # 0.5 degree grid in original study has 2880 cells in study area
    # sets threshold of 200 cells must meet threshold for MBF day 
    # 200/2880 simplifies to 5./72.
    nlons = len(lons)
    nlats = len(lats)
    cellno_threshold = 5./72. * nlons * nlats
        
    theta_equiv, theta_equiv_s = equiv_pot_t_ams(t, q)
    
    dthetady = theta_equiv.differentiate('lat') / a * 180./np.pi * 1000.
    
    lons = [dthetady.lon[i].values for i in range(len(dthetady.lon)) if dthetady.lon[i] >= 105. and dthetady.lon[i] <= 145.]
    lats = [dthetady.lat[i].values for i in range(len(dthetady.lat)) if dthetady.lat[i] >= 22.  and dthetady.lat[i] <= 40.]
    
    dthetady = np.abs(dthetady.sel(lon=lons, lat=lats))
    
    mbf_lats = dthetady.lat.where(dthetady > dthdy_threshold).mean('lat')
    
    mbfno = (dthetady>dthdy_threshold).sum(('lat','lon'))
    mbf_lats = mbf_lats.where(mbfno > cellno_threshold)
    
    continuous = np.abs(mbf_lats.diff('lon')).sum('lon')/(nlats-1) < continuity_threshold
    mbf_lats = mbf_lats.where(continuous==True)
    
    return mbf_lats, dthetady

 