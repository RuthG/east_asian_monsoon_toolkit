"""gao_etal_2001.py: 
Calculate the South China Sea monsoon onset index of Gao et al. 2001 (https://doi.org/10.3969/j.issn.1674-7097.2001.03.012)
NB. original paper is in Chinese, method is used as described in Martin et al. 2019 (https://doi.org/10.1007/s00376-018-8100-z), 
but returning to 335K (rather than 340K) threshold (based on tests with JRA-55 and CMIP6 climatologies)


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
    


def rename_coords(var_in, pentaddim='pentad', pdim='plev', punits='Pa', latdim='lat', londim='lon'):
    var_in = var_in.rename({pentaddim: 'pentad', latdim: 'lat', londim: 'lon'})

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


def area_mean(var_in, areacell=None):
    
    # Get specified lats and lons, select pentad range to look at
    lats = [u.lat[i].values for i in range(len(u.lat)) if u.lat[i] >= 10. and u.lat[i] <= 20.]
    lons = [u.lon[i].values for i in range(len(u.lon)) if u.lon[i] >= 110. and u.lon[i] <= 120.]
    
    # Take weighted means over above box
    if areacell is not None:
        areacell = areacell.rename({latdim:'lat', londim:'lon'})
        if (len(areacell.lat) == len(var_in.lat)) and (len(areacell.lon) == len(var_in.lon)):
            latsame = (np.round(areacell.lat.values - var_in.lat.values,2) ==0.).all()
            lonsame = (np.round(areacell.lon.values - var_in.lon.values,2) ==0.).all()
        else:
            latsame = False;  lonsame = False
        if latsame and lonsame:
            var_in_wt = var_in * areacell # If cell area grid provided use this to area weight values before averaging
            var_in_mean = var_in_wt.sel(lon=lons, lat=lats, method='nearest').sum(('lat','lon')) / areacell.sel(lat=lats, lon=lons, method='nearest').sum(('lat','lon')) 
        else:
            print('Warning, cell area dimensions do not match those of u, defaulting to using cosine weighted averaging')
            areacell = None
            
    if areacell is None:
        coslat = np.cos(var_in.lat * np.pi/180.)
        var_in_wt = var_in * coslat # In absence of cell area grid weight average by latitude before averaging
        var_in_mean = var_in_wt.sel(lon=lons, lat=lats, method='nearest').sum('lat').mean('lon') / coslat.sel(lat=lats, method='nearest').sum('lat')
    
    return var_in_mean
    



def gao_onset(t, q, u, areacell=None, pdim='plev', punits='Pa', pentaddim='pentad', latdim='lat', londim='lon'):
    # Inputs:
    # t, q, u: Single years of pentad mean temperature, specific humidity and zonal wind speed
    # areacell: Grid of cell areas for spatial averaging
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    # pdim: Name of pressure dimension, default plev
    # punits: Pressure units, default Pa
    # pentaddim: Name of pentad dimension, default 'pentad'
    
    # Returns:
    # SCSM onset pentad
    
    # Rename coords for ease if needed
    q = rename_coords(q, pdim=pdim, punits=punits, pentaddim=pentaddim, latdim=latdim, londim=londim)
    t = rename_coords(t, pdim=pdim, punits=punits, pentaddim=pentaddim, latdim=latdim, londim=londim)
    u = rename_coords(u, pdim=pdim, punits=punits, pentaddim=pentaddim, latdim=latdim, londim=londim)
    
    # Get equivalent potential temperature
    th_e, th_es = equiv_pot_t_ams(t, q)
    
    # Get area means from 10-20N, 110-120E
    th_e_mean = area_mean(th_e, areacell=areacell)
    u_mean = area_mean(u, areacell=areacell)
    u_pos = (u_mean > 0.) * 1.
    
    for i in range(len(u_mean.pentad)):  # For each pentad:
        if (th_e_mean.isel(pentad=i) > 335.) and (u_pos.isel(pentad=i) > 0.) and (u_pos.isel(pentad=i+1) > 0.): 
            # Check first if u and theta_eq exceed their thresholds, and u remains positive for at least 2 pentads
            # Now check steadiness. Is u either steady for 3 pentads with a break of 2 or fewer?
            if (u_pos.isel(pentad=i+2) > 0.) and (u_pos.isel(pentad=range(i+3,i+6)).sum('pentad') > 0.):
                return int(u_pos.pentad[i].values)
            # Or steady for 2 pentads with a break of only 1?
            elif (u_pos.isel(pentad=i+2) < 0.) and (u_pos.isel(pentad=[i+3]) > 0.):
                return int(u_pos.pentad[i].values)
            # If above are met, return pentad, otherwise continue to next pentad.
                
    return  None
