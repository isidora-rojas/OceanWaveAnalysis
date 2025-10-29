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

    # --- 1) Coerce to float32 & fill small gaps (linear)
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
