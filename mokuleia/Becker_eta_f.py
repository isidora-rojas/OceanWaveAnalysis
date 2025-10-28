"""
Wave setdown at reef face based on equation (4) in Becker (2014)
----------------------------------------------------------------

eta_f = - (H_f^2 * k_f) / (8 * sinh(2 * k_f * h_f))

Inputs (per time slice):
- Seta: surface elevation spectrogram, shape (n_freqs, n_times), units m^2/Hz
- freqs: frequency axis [Hz], shape (n_freqs,)
- t_eta: time centers for Seta's columns, shape (n_times,)
- h_eta: water depth at the reef face [m, positive], shape (n_times,)
- band: (fmin, fmax) [Hz] for band-limited computations

Outputs:
- xr.Dataset (default) or pandas.DataFrame with variables:
  eta_f, H_rms, k_f, h_f, f_rep  (indexed by time)
"""

from __future__ import annotations
import math
import numpy as np
import pandas as pd
import xarray as xr

from BulkWaveStats import wavenumber  # wavenumber(omega_array, depth_scalar) -> k_array


# --------------------------- helpers ---------------------------

def _as_2d_array(Seta) -> np.ndarray:
    """Return Seta as a 2D ndarray with shape (n_freqs, n_times)."""
    if isinstance(Seta, xr.DataArray):
        arr = np.asarray(Seta.values, dtype=float)
    else:
        arr = np.asarray(Seta, dtype=float)
    if arr.ndim != 2:
        raise ValueError("Seta must be 2D with shape (n_freqs, n_times).")
    return arr


def _representative_frequency(
    Seta: np.ndarray,
    freqs: np.ndarray,
    band: tuple[float, float],
) -> np.ndarray:
    """Energy-weighted mean frequency inside `band` for each time slice. Returns (n_times,)."""
    freqs = np.asarray(freqs, dtype=float)
    Seta = np.asarray(Seta, dtype=float)
    if Seta.shape[0] != freqs.size:
        raise ValueError("len(freqs) must equal Seta.shape[0] (n_freqs).")

    fmin, fmax = band
    mask = (freqs >= fmin) & (freqs <= fmax)
    if not np.any(mask):
        raise ValueError(f"Band {band} does not overlap provided frequencies.")

    fb = freqs[mask]
    Sb = Seta[mask, :]
    Sb = np.where(np.isfinite(Sb), Sb, 0.0)

    m0 = np.trapz(Sb, fb, axis=0)               # (n_times,)
    m1 = np.trapz(Sb * fb[:, None], fb, axis=0) # (n_times,)

    out = np.full(m0.shape, np.nan, dtype=float)
    good = m0 > 0.0
    out[good] = m1[good] / m0[good]
    return out


def _hrms_from_seta(
    Seta: np.ndarray,
    freqs: np.ndarray,
    band: tuple[float, float],
) -> np.ndarray:
    """
    Convert band-limited spectral energy to Hrms for each time slice.
    Hs = 4*sqrt(m0_band); Hrms = Hs/sqrt(2) = sqrt(8)*sqrt(m0_band).
    Returns (n_times,).
    """
    freqs = np.asarray(freqs, dtype=float)
    Seta = np.asarray(Seta, dtype=float)

    fmin, fmax = band
    mask = (freqs >= fmin) & (freqs <= fmax)
    if not np.any(mask):
        raise ValueError(f"Band {band} does not overlap provided frequencies.")

    fb = freqs[mask]
    Sb = Seta[mask, :]
    Sb = np.where(np.isfinite(Sb), Sb, 0.0)

    m0 = np.trapz(Sb, fb, axis=0)                # (n_times,)
    hrms = math.sqrt(8.0) * np.sqrt(np.maximum(m0, 0.0))
    hrms[~np.isfinite(hrms)] = np.nan
    return hrms


def _compute_kf(
    freq_rep: np.ndarray,
    depth_series: np.ndarray,
) -> np.ndarray:
    """
    Solve the linear dispersion relation for each time slice.
    freq_rep [Hz], depth_series [m] (positive). Returns (n_times,).
    """
    f = np.asarray(freq_rep, dtype=float)
    h = np.asarray(depth_series, dtype=float)
    if f.shape != h.shape:
        raise ValueError("freq_rep and depth_series must have the same shape (n_times,).")

    k_vals = np.full_like(f, np.nan, dtype=float)
    omega = 2.0 * np.pi * f

    for i, (om, depth) in enumerate(zip(omega, h)):
        if om > 0.0 and depth > 0.0 and np.isfinite(om) and np.isfinite(depth):
            k_vals[i] = wavenumber(np.array([om], dtype=float), float(depth))[0]
    return k_vals


# ------------------------- public API --------------------------

def compute_eta_f(
    Seta,
    freqs: np.ndarray,
    t_eta: np.ndarray,
    h_eta: np.ndarray,
    *,
    band: tuple[float, float] = (0.05, 0.33),
    output: str = "xarray",
):
    """
    Evaluate Becker (2014) Eq. (4) using Seta, freqs, t_eta, and h_eta (+band).

    Parameters
    ----------
    Seta : (n_freqs, n_times) ndarray or xarray.DataArray
    freqs : (n_freqs,) array-like [Hz]
    t_eta : (n_times,) array-like (datetime-like or numeric)
    h_eta : (n_times,) array-like [m, positive]
    band : (float, float), optional
    output : {"xarray","dataframe"}

    Returns
    -------
    xr.Dataset or pandas.DataFrame with variables:
        eta_f, H_rms, k_f, h_f, f_rep (indexed by time)
    """
    arr = _as_2d_array(Seta)             # (n_freqs, n_times)
    f = np.asarray(freqs, dtype=float)   # (n_freqs,)
    t = np.asarray(t_eta)
    h = np.asarray(h_eta, dtype=float)

    if arr.shape[0] != f.size:
        raise ValueError("len(freqs) must match Seta.shape[0] (n_freqs).")
    if arr.shape[1] != t.size or t.size != h.size:
        raise ValueError("Seta.shape[1], len(t_eta), and len(h_eta) must all match (n_times).")

    # core calculations
    Hrms = _hrms_from_seta(arr, f, band)                # (n_times,)
    f_rep = _representative_frequency(arr, f, band)     # (n_times,)
    k_vals = _compute_kf(f_rep, h)                      # (n_times,)

    # eta_f = -(Hrms^2 * k) / (8 * sinh(2 k h))
    denom = 8.0 * np.sinh(2.0 * k_vals * h)
    eta_vals = np.full_like(Hrms, np.nan, dtype=float)
    good = (denom != 0.0) & np.isfinite(k_vals) & np.isfinite(Hrms) & np.isfinite(h)
    eta_vals[good] = -((Hrms[good] ** 2) * k_vals[good]) / denom[good]

    # package results
    time_index = pd.to_datetime(t)

    attrs = {
        "equation": "eta_f = -(H_f^2 k_f) / (8 sinh(2 k_f h_f))",
        "band": tuple(float(x) for x in band),
        "freq_units": "Hz",
        "depth_units": "m (positive downward)",
        "Seta_units": "m^2/Hz",
        "notes": "Hrms from band-limited m0; f_rep is energy-weighted mean frequency in band.",
    }

    if output.lower() == "dataframe" or xr is None:
        df = pd.DataFrame(
            {"eta_f": eta_vals, "H_rms": Hrms, "k_f": k_vals, "h_f": h, "f_rep": f_rep},
            index=time_index,
        )
        df.index.name = "time"
        df.attrs = attrs
        return df

    ds = xr.Dataset(
        data_vars={
            "eta_f": ("time", eta_vals),
            "H_rms": ("time", Hrms),
            "k_f": ("time", k_vals),
            "h_f": ("time", h),
            "f_rep": ("time", f_rep),
        },
        coords={"time": time_index},
        attrs=attrs,
    )
    return ds


__all__ = ["compute_eta_f"]