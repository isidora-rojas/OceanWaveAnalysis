''' This code is used to detide pressure/surface elevation/ current data using the UTide package '''

from utide import solve, reconstruct
import numpy as np
import pandas as pd



def detide(
    df: pd.DataFrame, 
    col: str, 
    LAT, 
    avg = '15min', 
    inplace: bool = False
    ):

    ''' Removes the tidal signal given a pressure/surface elevation/current dataframe. Averages are computed over specifief time windows. This
    is necessary because the computation can be computationally expensive for longer time series. 

    Adds three columns:
        <prefix>_tide     : tidal signal reconstructed on 1 Hz (or native) grid
        <prefix>_interp   : original series linearly interpolated (float32)
        <prefix>_detided  : <prefix>_interp - <prefix>_tide (wave-only signal)

    Parameters
    -----------------
    df: pd.DataFrame
        Dataframe with datetime index and pressure/surface elevation/current column to detide 
    col: str
        column name of pressure/surface elevation/current to detide 
    LAT: float
        latitude of sensor
    avg: string
        time average needed. Becker (2014) computes 15 min averages
    inplace: bool
        write into df (True) or return a copy (False)

    Requirements
    ------------
    - df.index must be a DatetimeIndex at (about) 1 Hz (or at least regular).
    - df[col] is the pressure / surface elevation / current to detide.

    Returns
    -------------
    tide_full: ndarray
        tidal signal
    p_prime: ndarray
        pressure from waves. That is, the full array with tide subtracted
    p_full: ndarray
        the original array with any NaNs interpolated

    '''
    # --- check inputs ---
    if not isinstance(df.index, pd.DatetimeIndex):
        raise TypeError("detide: df.index must be a DatetimeIndex.")
    if col not in df.columns:
        raise KeyError(f"detide: column '{col}' not found in DataFrame.")
    if not df.index.is_monotonic_increasing:
        # keep it simple: just sort
        df = df.sort_index()
    out = df if inplace else df.copy()


    LAT = LAT

    # fill gaps where needed. Assumes regular 1hz data but there should not be any nans unless tide is really low   
    s = out[col].astype("float32")
    s_interp = (
        s.to_frame(name=col)
         .assign(_tmp=s)
         .pop("_tmp")
         .interpolate(method="linear", limit_direction="both")
         .astype("float32")
    )
    out[f"{col}_interp"] = s_interp 

    #  15min means to help with computation of tides
    s_mean = s_interp.resample(avg).mean().dropna()
    t_mean = s_mean.index.to_numpy()

    # tideal computation on grid averaged to 'avg'
    coef = solve(
        t_mean,
        s_mean.to_numpy(),
        lat=LAT,
        method="ols",
        conf_int="linear",   # quick CI; set 'none' if you want even faster
        nodal=True,
        trend=True,
        verbose=False,
    ) 

    # # Reconstruct on the 1 Hz grid (chunking to avoid crashes)
    # --- Chunked UTide reconstruction directly on DataFrame index ---

    tide_parts = []
    for idx_chunk in np.array_split(out.index, 8):   # 8 chunks, adjust as needed
        rc = reconstruct(idx_chunk.to_numpy(), coef)  # UTide reconstruct call
        tide_parts.append(pd.Series(rc.h.astype("float32"), index=idx_chunk))

    tide_full = pd.concat(tide_parts).sort_index()
    out[f"{col}_tide"] = tide_full

    # detided signal (wave + setup + IG preserved)
    out[f"{col}_detided"] = (out[f"{col}_interp"] - out[f"{col}_tide"]).astype("float32")

    return out