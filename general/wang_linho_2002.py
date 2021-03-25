"""wang_linho_2002.py:  https://doi.org/10.1175/1520-0442(2002)015%3C0386:RSOTAP%3E2.0.CO;2
Calculate the precipitation index of Wang and LinHo 2002:
Climatological monsoon season defined as pentad at which precipitation is 5mm/day > January (NH) or July (SH) mean
Onset is pentad when this first occurs, withdrawal when it ends, 
duration difference between these, peak pentad of max difference.

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr

def wang_linho(p_pentad, p_month, latdim='lat', londim='lon', pentaddim='pentad'):
    # Inputs:
    # p_pentad: lat-lon-time DataArray of precipitation, units mm/day, time units pentad
    # p_month: lat-lon-time DataArray of precipitation, units mm/day, average over Jan or July
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    
    # Returns:
    # Onset, peak and withdrawal pentads, and season duration in pentads
    
    # Smooth data
    p_mean = p_pentad.mean(pentaddim) # save mean
    p_pentad = p_pentad - p_mean # subtract
    axisno = p_pentad.get_axis_num(pentaddim) # Find number of pentad dimension
    p_fft = np.fft.fft(p_pentad.values, axis=axisno) # Take fourer transform along pentad dimension - use fft not rfft to conserve axis length
    # Discard all but 1st 12 harmonics
    if axisno==0:
        p_fft[12:-12,:,:] = 0
    elif axisno==1:
        p_fft[:,12:-12,:] = 0
    elif axisno==2:
        p_fft[:,:,12:-12] = 0
    p_smooth = np.fft.ifft(p_fft, axis=axisno) # Take inverse fourier transform to recover smoothed timeseries
    p_smooth = xr.DataArray(p_smooth.real, coords=p_pentad.coords, dims=p_pentad.dims) # Make dataarray
    p_pentad = p_smooth + p_mean # Add mean back on
    
    rain_rel = p_pentad - p_month
    
    rain_rel_masked = np.ma.masked_less(rain_rel.values, 5)  # Mask relative rainfall where it is less than 5mm/day

    onset_index = np.ma.notmasked_edges(rain_rel_masked, axis=0) # Look along pentad axis to find edges of mask
    onset = np.zeros((len(rain_rel[latdim]),len(rain_rel[londim]))) # Create array of nans to load onsets into
    onset[:] = np.nan
    withdrawal = np.zeros((len(rain_rel[latdim]),len(rain_rel[londim]))) # Create array of nans to load onsets into
    withdrawal[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])): # Extract onsets from mask edges
        onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1
        withdrawal[ onset_index[1][1][i], onset_index[1][2][i] ] = onset_index[1][0][i]+1
        
    onset_pentad = xr.DataArray(onset, [(latdim, rain_rel[latdim]), (londim, rain_rel[londim])]) # Make dataarray
    withdrawal_pentad = xr.DataArray(withdrawal, [(latdim, rain_rel[latdim]), (londim, rain_rel[londim])]) # Make dataarray
    peak_pentad = xr.DataArray(rain_rel_masked.argmax(axis=0) + 1., [(latdim, rain_rel[latdim]), (londim, rain_rel[londim])]) # Make dataarray
    peak_pentad = peak_pentad.where(peak_pentad!=1.)
    
    duration = withdrawal_pentad - onset_pentad
    
    return onset_pentad, withdrawal_pentad, peak_pentad, duration
   
   
