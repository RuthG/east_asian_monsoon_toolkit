"""kitoh_uchiyama_2006.py: 
Calculate the precipitation index of Kitoh and Uchiyama 2006: https://doi.org/10.2151/jmsj.84.247
Climatological monsoon season defined as pentad at which Normalised Pentad Precipitation Index exceeds 0.618

Onset is pentad when this first occurs, withdrawal when it ends, duration difference between these, peak pentad of max difference.

Requires xarray, numpy"""
__author__      = "Ruth Geen"


import numpy as np
import xarray as xr

def kitoh_uchiyama(p_pentad, latdim='lat', londim='lon', pentaddim='pentad'):
    # Inputs:
    # p_pentad: lat-lon-time DataArray of precipitation, units mm/day, time units pentad
    # latdim: Name of latitude dimension, default lat
    # londim: Name of longitude dimension, default lon
    
    # Returns:
    # Onset, peak and withdrawal pentads, and season duration in pentads
        
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
    p_smooth = p_smooth + p_mean # Add mean back on
        
    npi = (p_smooth - p_smooth.min(pentaddim)) / (p_smooth.max(pentaddim) - p_smooth.min(pentaddim)) # Evaluate Normalised Pentad Precipitation Index
        
    npi_masked = np.ma.masked_less(npi.values, 0.618)  # Mask npi where it is less than 0.618

    onset_index = np.ma.notmasked_edges(npi_masked, axis=0) # Look along pentad axis to find edges of mask
    onset = np.zeros((len(npi[latdim]),len(npi[londim]))) # Create array of nans to load onsets into
    onset[:] = np.nan
    withdrawal = np.zeros((len(npi[latdim]),len(npi[londim]))) # Create array of nans to load onsets into
    withdrawal[:] = np.nan
    
    for i in range(0,len(onset_index[0][1])): # Extract onsets from mask edges
        onset[ onset_index[0][1][i], onset_index[0][2][i] ] = onset_index[0][0][i]+1
        withdrawal[ onset_index[1][1][i], onset_index[1][2][i] ] = onset_index[1][0][i]+1
        
    onset_pentad = xr.DataArray(onset, [(latdim, npi[latdim]), (londim, npi[londim])]) # Make dataarray
    withdrawal_pentad = xr.DataArray(withdrawal, [(latdim, npi[latdim]), (londim, npi[londim])]) # Make dataarray
    peak_pentad = xr.DataArray(npi_masked.argmax(axis=0) + 1., [(latdim, npi[latdim]), (londim, npi[londim])]) # Make dataarray
    peak_pentad = peak_pentad.where(peak_pentad!=1.)
    
    duration = withdrawal_pentad - onset_pentad
    
    return onset_pentad, withdrawal_pentad, peak_pentad, duration
   
