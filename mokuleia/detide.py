''' This code is used to detide pressure/surface elevation/ current data using the UTide package '''

from utide import solve, reconstruct
import numpy as np
import pandas as pd



def detide(
    p: np.ndarray,
    t: np.ndarray,
    LAT, 
    avg = '15min'
    ):

    ''' Removes the tidal signal given a pressure/surface elevation/current array. Averages are computed over specifief time windows. This
    is necessary because the computation can be computationally expensive for longer time series. 

    Parameters
    -----------------
    p: ndarray
        pressure/surface elevation/ current array
    t: ndarray
        time array
    LAT: float
        latitude of sensor
    avg: string
        time average needed. Becker (2014) computes 15 min averages

    Returns
    -------------
    tide_full: ndarray
        tidal signal
    p_prime: ndarray
        pressure from waves. That is, the full array with tide subtracted
    p_full: ndarray
        the original array with any NaNs interpolated

    '''


    LAT = LAT

    idx_1hz = pd.to_datetime(t)
    p_interp= (
        pd.Series(p, index=idx_1hz)
        .interpolate(method='linear', limit_direction='both')
        .astype('float32') # make 32 bit to stop pooter from crashing
    )   

    #  15min means to help with computation of tides
    pmean = p_interp.resample(avg).mean().dropna()
    idx_mean = pmean.index

    # tideal computation on 30s grid
    coef = solve(
        idx_mean.to_numpy(),
        pmean.to_numpy(),
        lat=LAT,
        method="ols",
        conf_int="linear", # 'none' skips uncertainty calc, 'linear' uses quick lin. estimate, 'MC' run monte-carlo (mad spendy)
        nodal=True,
        trend=True,
        verbose=False)

    # # Reconstruct on the 1 Hz grid (chunking to avoid crashes)
    chunks = np.array_split(idx_1hz, 8)  # adjust slice count for memory
    tide_parts = []
    for chunk in chunks:
        # chunk is a DatetimeIndex
        recon_chunk = reconstruct(chunk.to_numpy(), coef)
        tide_parts.append(pd.Series(recon_chunk.h.astype("float32"), index=chunk))
    tide_full = pd.concat(tide_parts).sort_index()
    p_wave = (p_interp - tide_full) # water level chane due to wave effects only
    p_detide = p_wave.resample('1s').mean().dropna()

    return tide_full, p_detide, p_interp
