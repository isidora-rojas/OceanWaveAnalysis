"""
Wave setdown at reef face based on equation (4) in Becker (2014)
-----------------------------------------------------------

eta_f = - (H_f^2 * k_f) / (8 * sinh(2 * k_f * h_f))

This module provides helpers to derive each term from a surface-elevation
spectrogram (`Seta`) together with the companion frequency and time axes.
"""

from __future__ import annotations
import math
from typing import Iterable, Tuple
import numpy as np
import pandas as pd
try:  # pragma: no cover - xarray is optional
    import xarray as xr
except ImportError:  # pragma: no cover - fall back to pandas
    xr = None
from BulkWaveStats import sig_wave_height, wavenumber


def _depth_at_centers(
    t_spec: np.ndarray,
    t1: np.ndarray,
    h1: np.ndarray,
    *,
    depth_interp: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate the depth record to the spectrogram time centres."""
    time_index = pd.to_datetime(t1)
    if depth_interp is None:
        depth_interp = (
            pd.Series(h1, index=time_index)
            .interpolate(method="time", limit_direction="both")
            .to_numpy()
        )
    else:
        depth_interp = np.asarray(depth_interp, dtype=float)
        if depth_interp.shape != np.asarray(h1).shape:
            raise ValueError("depth_interp must match the shape of h1")

    t0 = np.array(time_index[0], dtype="datetime64[ns]")
    seconds_full = (time_index.to_numpy() - t0) / np.timedelta64(1, "s")
    depth_at_centers = np.interp(t_spec, seconds_full, depth_interp)
    time_centers = t0 + (t_spec * 1e9).astype("timedelta64[ns]")
    return depth_at_centers, time_centers


def _representative_frequency(
    Seta: np.ndarray,
    freqs: np.ndarray,
    band: tuple[float, float],
) -> np.ndarray:
    """Energy-weighted mean frequency inside `band` for each time slice."""
    freqs = np.asarray(freqs, dtype=float)
    Seta = np.asarray(Seta, dtype=float)
    if Seta.shape[0] != freqs.size:
        raise ValueError("Seta must have shape (n_freqs, n_windows)")

    mask = (freqs >= band[0]) & (freqs <= band[1])
    if not np.any(mask):
        raise ValueError(f"Band {band} does not overlap provided frequencies")

    freqs_band = freqs[mask]
    spectra_band = Seta[mask, :]

    is_finite = np.isfinite(spectra_band)
    valid_cols = np.any(is_finite, axis=0)

    spectra_band = np.where(is_finite, spectra_band, 0.0)

    # Zeroth and first spectral moments in the selected band.
    m0 = np.trapz(spectra_band, freqs_band, axis=0)
    m1 = np.trapz(spectra_band * freqs_band[:, None], freqs_band, axis=0)

    freq_rep = np.full(m0.shape, np.nan, dtype=float)
    positive = (m0 > 0.0) & valid_cols
    freq_rep[positive] = m1[positive] / m0[positive]
    return freq_rep


def compute_hrms(
    Seta: np.ndarray,
    freqs: np.ndarray,
    t_spec: np.ndarray,
    t1: np.ndarray,
    h1: np.ndarray,
    band: tuple[float, float],
    *,
    depth_interp: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert the band-limited significant wave height to Hrms.
    """
    _, _, _, hs_band, _, time_centers = sig_wave_height(
        Seta,
        freqs,
        t_spec,
        t1,
        h1,
        band,
        depth_interp=depth_interp,
        dataviz=False,
    )
    hrms = 0.5 * math.sqrt(2.0) * hs_band
    return hrms, time_centers


def compute_kf(
    freq_rep: np.ndarray,
    depth_series: np.ndarray,
) -> np.ndarray:
    """
    Solve the dispersion relation for each time slice.

    Parameters
    ---------------
    freq_rep: np.ndarray
        median frequency array for each time slice
    depth_series: np.ndarray
        !!MUST BE POSITIVE VALUES!!
    """
    freq_rep = np.asarray(freq_rep, dtype=float)
    depth_series = np.asarray(depth_series, dtype=float)

    k_vals = np.zeros_like(freq_rep, dtype=float)
    omega_vals = 2.0 * np.pi * freq_rep

    

    for idx, (omega, depth) in enumerate(zip(omega_vals, depth_series, strict=True)):
        if omega <= 0.0 or depth <= 0.0 or not np.isfinite(omega) or not np.isfinite(depth):
            k_vals[idx] = np.nan
            print('There are negative values in the input arrays.')
        else:
            k_vals[idx] = wavenumber(np.array([omega]), float(depth))[0]
    return k_vals



def compute_eta_f(
    Seta: np.ndarray,
    freqs: np.ndarray,
    t_spec: np.ndarray,
    t1: np.ndarray,
    h1: np.ndarray,
    *,
    band: Tuple[float, float] = (0.05, 0.33),
    depth_interp: np.ndarray | None = None,
    output: str = "xarray",
):
    """
    Evaluate equation (4) from Becker (2014) for a given spectrogram.

    Parameters
    ----------
    Seta, freqs, t_spec, t1, h1
        Surface elevation spectra and metadata produced by `Spp_to_Seta`.
    band
        Frequency interval [Hz] used to compute Hrms and representative
        frequency. Defaults to the sea/swell band.
    depth_interp
        Optionally supply a pre-interpolated depth series to avoid repeating
        the interpolation performed in other routines.
    output
        ``"xarray"`` (default) returns an ``xr.Dataset`` when xarray is
        available. Set to ``"dataframe"`` to always receive a
        ``pandas.DataFrame``.
    """
    hrms, time_centers = compute_hrms(
        Seta,
        freqs,
        t_spec,
        t1,
        h1,
        band,
        depth_interp=depth_interp,
    )
    depth_at_centers, time_centers_check = _depth_at_centers(
        t_spec,
        t1,
        h1,
        depth_interp=depth_interp,
    )
    if not np.all(time_centers == time_centers_check):
        # Align arrays if sig_wave_height and our interpolation applied
        # slightly different rounding.
        time_centers = time_centers_check

    freq_rep = _representative_frequency(Seta, freqs, band)
    k_vals = compute_kf(freq_rep, depth_at_centers)

    denominator = 8.0 * np.sinh(2.0 * k_vals * depth_at_centers)
    eta_vals = np.full_like(hrms, np.nan, dtype=float)
    mask = (denominator != 0.0) & np.isfinite(k_vals) & np.isfinite(hrms)
    eta_vals[mask] = -((hrms[mask] ** 2) * k_vals[mask]) / denominator[mask]

    if output.lower() == "dataframe" or xr is None:
        df = pd.DataFrame(
            {
                "eta_f": eta_vals,
                "H_rms": hrms,
                "k_f": k_vals,
                "h_f": depth_at_centers,
                "f_rep": freq_rep,
            },
            index=pd.to_datetime(time_centers),
        )
        df.index.name = "time"
        df.attrs = {
            "equation": "eta_f = -(H_f^2 k_f) / (8 sinh(2 k_f h_f))",
            "band": band,
        }
        return df

    dataset = xr.Dataset(
        data_vars={
            "eta_f": ("time", eta_vals),
            "H_rms": ("time", hrms),
            "k_f": ("time", k_vals),
            "h_f": ("time", depth_at_centers),
            "f_rep": ("time", freq_rep),
        },
        coords={"time": pd.to_datetime(time_centers)},
        attrs={
            "equation": "eta_f = -(H_f^2 k_f) / (8 sinh(2 k_f h_f))",
            "band": band,
        },
    )
    return dataset


__all__ = [
    "compute_eta_f",
    "compute_hrms",
    "compute_kf",
]
