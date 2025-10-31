from utide import solve, reconstruct
import numpy as np
import pandas as pd

def detide_df(
    df: pd.DataFrame,
    col: str,
    lat: float,
    *,
    avg: str = "15min",         # Becker(2014) style
    chunks: int = 8,            # chunked reconstruct to save RAM
    out_prefix: str = None,     # default = col name
    inplace: bool = False,      # write into df or return a copy
):
    """
    Detide a DataFrame column using UTide.
    
    Adds three columns:
      <prefix>_tide     : tidal signal reconstructed on 1 Hz (or native) grid
      <prefix>_interp   : original series linearly interpolated (float32)
      <prefix>_detided  : <prefix>_interp - <prefix>_tide (wave-only signal)

    Requirements
    ------------
    - df.index must be a DatetimeIndex at (about) 1 Hz (or at least regular).
    - df[col] is the pressure / surface elevation / current to detide.
    """
    if out_prefix is None:
        out_prefix = col

    target = df if inplace else df.copy()

    # --- 1) downsize to  32 & fill small gaps (linear)
    if not isinstance(target.index, pd.DatetimeIndex):
        raise TypeError("detide_df: DataFrame index must be a DatetimeIndex.")
    s = target[col].astype("float32")
    s_interp = (
        s.reindex(pd.DatetimeIndex(target.index))  # ensure proper index dtype
         .interpolate(method="time", limit_direction="both")
         .astype("float32")
    )

    # --- 2) Fit UTide on downsampled means (to keep it fast/stable)
    s_mean = s_interp.resample(avg).mean().dropna()
    t_mean = s_mean.index.to_numpy()

    coef = solve(
        t_mean,
        s_mean.to_numpy(),
        lat=lat,
        method="ols",
        conf_int="linear",
        nodal=True,
        trend=True,
        verbose=False,
    )

    # --- 3) Reconstruct on full grid in chunks (avoid memory spikes)
    t_full = target.index
    parts = []
    for chunk in np.array_split(t_full, chunks):
        if len(chunk) == 0:
            continue
        rc = reconstruct(chunk.to_numpy(), coef)
        parts.append(pd.Series(rc.h.astype("float32"), index=chunk))
    tide_full = pd.concat(parts).sort_index()

    # --- 4) Assemble outputs
    target[f"{out_prefix}_interp"]  = s_interp
    target[f"{out_prefix}_tide"]    = tide_full
    target[f"{out_prefix}_detided"] = (s_interp - tide_full).astype("float32")

    return target




from utide import solve, reconstruct
import numpy as np
import pandas as pd

def detide(
    df: pd.DataFrame,
    col: str,
    LAT: float,
    avg: str = "15min",
    inplace: bool = False,
) -> pd.DataFrame:
    """
    Remove tidal signal from a pressure/surface elevation/current time series.

    Adds:
        <col>_interp   : linearly interpolated series on native grid (float32)
        <col>_tide     : UTide reconstruction on native grid (float32)
        <col>_detided  : <col>_interp - <col>_tide (float32)

    Requirements:
      - df.index is a DatetimeIndex (regular-ish; UTide can handle datetimes)
      - df[col] is the series to detide
    """
    if not isinstance(df.index, pd.DatetimeIndex):
        raise TypeError("detide: df.index must be a DatetimeIndex.")
    if col not in df.columns:
        raise KeyError(f"detide: column '{col}' not found in DataFrame.")
    if not df.index.is_monotonic_increasing:
        # keep it simple: just sort
        df = df.sort_index()

    out = df if inplace else df.copy()

    # 1) Light interpolation to get a continuous series on the native grid
    s = out[col].astype("float32")
    s_interp = (
        s.to_frame(name=col)
         .assign(_tmp=s)
         .pop("_tmp")
         .interpolate(method="linear", limit_direction="both")
         .astype("float32")
    )
    out[f"{col}_interp"] = s_interp

    # 2) Coarsen to 'avg' for UTide coefficient solve (faster & stable)
    s_mean = s_interp.resample(avg).mean().dropna()
    t_mean = s_mean.index.to_numpy()

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

    # 3) Reconstruct tide on the native grid
    t_full = out.index.to_numpy()
    recon = reconstruct(t_full, coef)
    tide_full = pd.Series(recon.h.astype("float32"), index=out.index, name=f"{col}_tide")
    out[f"{col}_tide"] = tide_full

    # 4) Wave-only / detided signal on native grid
    out[f"{col}_detided"] = (out[f"{col}_interp"] - out[f"{col}_tide"]).astype("float32")

    return out
