import numpy as np
import pandas as pd
import xarray as xr
from scipy.signal import spectrogram
from BulkWaveStats import Spp_to_Seta

# Assumes you have:
# from BulkWaveStats import Spp_to_Seta
# from your_module import compute_eta_f   # (the Becker 2014 Eq.4 implementation)

def make_sensor_ds(
    df: pd.DataFrame,
    *,
    pressure_detided_col: str = "pressure_pa_detided",
    depth_col: str = "depth_m",
    fs: float = 1.0,
    window_len: str = "15min",         # 15-min products
    win: str = "hann",
    band: tuple = (0.05, 0.33),        # sea/swell by default; use (0.004, 0.04) for IG
):
    """
    Build a tidy xarray.Dataset of 15-min spectral products for ONE sensor:
      - Spp(freq,time), Seta(freq,time), h_mean(time)
      - eta_f(time), H_rms(time), k_f(time) from Becker (2014) Eq. (4)

    Inputs come directly from the detided DataFrame `df`:
      df[pressure_detided_col] : 1 Hz detided pressure
      df[depth_col]            : 1 Hz depth (m, positive down)
    """
    # --- 15-min mean depth, broadcast to 1 Hz (piecewise constant)
    h15 = df[depth_col].resample(window_len, label="left", closed="left").mean()
    h15_1Hz = h15.reindex(df.index, method="pad")

    # --- Spectrogram (non-overlapping 15-min windows)
    nperseg = int(pd.to_timedelta(window_len).total_seconds() * fs)
    noverlap = 0
    x = df[pressure_detided_col].astype(float).to_numpy()

    freqs, t_spec_sec, Spp = spectrogram(
        x,
        fs=fs,
        window=win,
        nperseg=nperseg,
        noverlap=noverlap,
        detrend="linear",
        scaling="density",
        mode="psd",
    )

    # Map spectrogram times (seconds from start) -> datetimes (bin centers)
    t0 = df.index[0]
    t_eta = pd.to_datetime(t0.value + (t_spec_sec * 1e9).astype("int64"))

    # Convert Spp -> Seta, forcing 15-min mean depth per window via depth_interp
    Seta, t_eta_chk, h_eta = Spp_to_Seta(
        Spp,
        freqs,
        t_eta,
        df.index,
        df[depth_col].to_numpy(),
        depth_interp=h15_1Hz.to_numpy(),
    )

    # Package spectra
    ds_spec = xr.Dataset(
        data_vars={
            "Spp":  (("freq", "time"), Spp),
            "Seta": (("freq", "time"), Seta),
            "h_mean": ("time", h_eta),
        },
        coords={
            "freq": ("freq", freqs, {"units": "Hz"}),
            "time": ("time", pd.to_datetime(t_eta), {"long_name": f"{window_len} bin center"}),
        },
        attrs={
            "description": f"{window_len} spectral products (detided); Spp→Seta via linear wave theory",
            "window": win,
            "nperseg": int(nperseg),
            "noverlap": int(noverlap),
            "sampling_hz": fs,
            "Seta_units": "m^2/Hz",
            "Spp_units": "Pa^2/Hz",
            "depth_units": "m (positive down)",
            "depth_definition": f"{window_len} mean depth broadcast to 1 Hz before interpolation",
        },
    )

    return ds_spec


import numpy as np
import pandas as pd
import xarray as xr

def combine_eta_i(
    ds_face: xr.Dataset,
    ds_inner: xr.Dataset,
    *,
    depth_face_var: str = "h_f",     # fall back to 'h_mean' if not present
    depth_inner_var: str = "h_mean",
    tolerance: str = "2min",
    round_to: str | None = None,     # e.g. "1min" if your centers are messy
    drop_unmatched: bool = True,
) -> xr.Dataset:
    """
    Compute inner-reef setup on the face time grid:
        eta_i = eta_f + (h_i - h_f)

    - Aligns h_i to face times with nearest+tolerance using REINDEX (no KeyError).
    - Forces depths positive.
    - Optionally rounds time coords before alignment.
    """

    # choose face depth var
    if depth_face_var in ds_face:
        h_f = ds_face[depth_face_var]
    else:
        h_f = ds_face["h_mean"]
    h_f = np.abs(h_f).rename("h_f")

    # inner depth
    h_i = np.abs(ds_inner[depth_inner_var]).rename("h_i")

    # (optional) round times to improve matches
    tf = pd.to_datetime(ds_face.time.values)
    ti = pd.to_datetime(ds_inner.time.values)
    if round_to is not None:
        tf = pd.to_datetime(pd.DatetimeIndex(tf).round(round_to))
        ti = pd.to_datetime(pd.DatetimeIndex(ti).round(round_to))

    # reindex inner depth onto face times by nearest with tolerance (no exception on misses)
    tol = pd.to_timedelta(tolerance) if tolerance is not None else None
    h_i_on_face = h_i.copy()
    h_i_on_face = h_i_on_face.assign_coords(time=ti).reindex(
        time=tf, method="nearest", tolerance=tol
    )

    # build output on face grid
    out = xr.Dataset(coords={"time": tf})
    out["eta_f"] = ds_face["eta_f"].assign_coords(time=tf)
    out["h_f"]   = h_f.assign_coords(time=tf)
    out["h_i"]   = h_i_on_face

    # optionally drop times where inner didn't match within tolerance
    if drop_unmatched:
        valid = out["h_i"].notnull()
        out = out.sel(time=out.time.where(valid, drop=True))

    # eta_i = eta_f + (h_i - h_f)
    out["eta_i"] = (out["eta_f"] + (out["h_i"] - out["h_f"])).assign_attrs(
        {"units": "m", "long_name": "inner-reef setup"}
    )

    # keep attrs tidy
    out["h_f"].attrs.update({"units": "m", "positive": "down"})
    out["h_i"].attrs.update({"units": "m", "positive": "down"})
    out["eta_f"].attrs.update({"units": "m", "long_name": "reef-face setdown"})
    return out

