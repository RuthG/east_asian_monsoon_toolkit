# Readme for East Asian Monsoon Toolkit

Project contains code for a range of indices/metrics for monsoon intensity, onset, withdrawal and duration.
Code requires numpy and xarray

Code is organised by region/monsoon feature:

## easm: East Asian Summer Monsoon
- *wang_fan_1999.py*: Intensity index from https://doi.org/10.1175/1520-0477(1999)080%3C0629:COSASM%3E2.0.CO;2

## scsm: South China Sea Monsoon
- *wang_etal_2004.py*: Onset index described in https://doi.org/10.1175/2932.1
- *gao_etal_2001.py*: Onset index described in https://doi.org/10.1007/s00376-018-8100-z

## meiyu: Meiyu-Baiu rainband
- *li_etal_2018.py*: Front identification as described in https://doi.org/10.1007/s00382-017-3975-4

## general: Indices broadly applicable across Asia/globe
- *wang_linho_2002.py*: Identify pentads of monsoon onset, withdrawal and peak based on value relative to Jan climatological mean https://doi.org/10.1175/1520-0442(2002)015%3C0386:RSOTAP%3E2.0.CO;2
- *kitoh_uchiyama_2006.py*: Identify pentads of monsoon onset, withdrawal and peak based on Normalised Pentad Precipitation Index https://doi.org/10.2151/jmsj.84.247
- *li_zhang_2009.py*: Identify pentads of monsoon onset and withdrawal based on wind rotation criteria https://doi.org/10.1007/s00382-008-0465-8
- *yim_etal_2014_vorticity.py*: Evaluate vorticity/wind metrics as a proxy for intensity of various regional monsoons https://doi.org/10.1007/s00382-013-1956-9
- *yim_etal_2014_precip.py*: Calculate regional monsoon precipitation intensities as per https://doi.org/10.1007/s00382-013-1956-9


## Work in progress:
- Clean up code further, aim to make input data format needed consistent
- Add further tools as useful, e.g. list in https://doi.org/10.1175/2008JCLI2183.1 
- Expand to South Asian and other regional monsoons