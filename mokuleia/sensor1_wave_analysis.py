#!/usr/bin/env python3
"""
End-to-end spectral processing for Mokuleia deployment A, sensor 1.

The script reproduces the wave bulk parameters shown in PILOTmokuleia.pdf
using the pressure record stored in MoA18411.mat.  It

1. loads the raw 12 h bursts and inserts the 20 s download gaps,
2. converts pressure from psi to gauge Pascals and then to water depth,
3. removes the low-frequency tide with a moving mean,
4. interpolates across the short gaps, and
5. computes sliding spectra of surface elevation along with bulk metrics.

Figures are written to mokuleia/figures/.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.signal import spectrogram

# physical constants
PSI_TO_PA = 6894.757293168  # [Pa/psi]
RHO_SEAWATER = 1025.0  # [kg/m^3]
G = 9.80665  # [m/s^2]


@dataclass
class SensorSeries:
    """Container for the cleaned 1 Hz time series."""

    time: np.ndarray  # datetime64[ns]
    seconds: np.ndarray  # seconds since start, float64
    pressure_pa: np.ndarray  # gauge pressure, Pa
    depth_m: np.ndarray  # instantaneous water depth, m
    pressure_hp: np.ndarray  # high-passed pressure, Pa
    depth_lp: np.ndarray  # low-passed depth (tide), m


@dataclass
class SpectralSummary:
    """Sliding spectral diagnostics for quick plotting."""

    freqs: np.ndarray
    surface_spectra: np.ndarray  # shape (n_windows, n_freqs)
    time_centers: np.ndarray  # datetime64[ns]
    hs_total: np.ndarray
    hs_sea_swell: np.ndarray
    hs_infragravity: np.ndarray
    tp_sea_swell: np.ndarray


def matlab_datenum_to_datetime64(tnums: np.ndarray) -> np.ndarray:
    """Convert MATLAB datenum to numpy datetime64[ns]."""
    tnums = np.asarray(tnums, dtype=np.float64)
    days_since_unix = tnums - 719529.0
    ns_since_unix = np.round(days_since_unix * 86400.0 * 1e9).astype("int64")
    return np.datetime64("1970-01-01T00:00:00") + ns_since_unix.astype("timedelta64[ns]")


def moving_average_nan(x: np.ndarray, window: int) -> np.ndarray:
    """Centered moving average that ignores NaNs."""
    if window < 1:
        raise ValueError("window must be >= 1")
    x = np.asarray(x, dtype=np.float64)
    kernel = np.ones(window, dtype=np.float64)
    valid = np.isfinite(x).astype(np.float64)
    filled = np.where(np.isfinite(x), x, 0.0)
    num = np.convolve(filled, kernel, mode="same")
    den = np.convolve(valid, kernel, mode="same")
    out = np.full_like(x, np.nan)
    ok = den > 0
    out[ok] = num[ok] / den[ok]
    return out


def nan_interp(x: np.ndarray) -> np.ndarray:
    """Interpolate across NaNs using first-order hold."""
    x = np.asarray(x, dtype=np.float64)
    out = x.copy()
    isnan = ~np.isfinite(out)
    if not np.any(isnan):
        return out
    idx = np.arange(out.size)
    valid = ~isnan
    if not np.any(valid):
        raise ValueError("cannot interpolate array of all NaNs")
    out[isnan] = np.interp(idx[isnan], idx[valid], out[valid])
    return out


def wavenumber(omega: np.ndarray, depth: float, tol: float = 1e-12, max_iter: int = 64) -> np.ndarray:
    """Solve the linear dispersion relation for k(ω)."""
    omega = np.asarray(omega, dtype=np.float64)
    depth = float(depth)
    k = np.zeros_like(omega)
    if depth <= 0.0:
        return k
    mask = omega > 0.0
    if not np.any(mask):
        return k
    k_mask = (omega[mask] ** 2) / G  # deep-water guess
    for _ in range(max_iter):
        kh = k_mask * depth
        tanh_kh = np.tanh(kh)
        cosh_kh = np.cosh(kh)
        sech_kh_sq = 1.0 / (cosh_kh ** 2)
        f = G * k_mask * tanh_kh - omega[mask] ** 2
        df = G * tanh_kh + G * depth * k_mask * sech_kh_sq
        step = np.divide(f, df, out=np.zeros_like(f), where=df != 0.0)
        k_next = k_mask - step
        if np.nanmax(np.abs(step)) < tol:
            break
        k_mask = np.where(np.isfinite(k_next), k_next, k_mask)
    k[mask] = k_mask
    return k


def load_sensor_series(mat_path: Path, patm_psi: float = 14.7, tide_window: int = 3600, wave_window: int = 600) -> SensorSeries:
    """
    Load the 12 h bursts, insert 20 s download gaps, and derive helper series.

    Parameters
    ----------
    mat_path : Path
        Location of MoA18411.mat.
    patm_psi : float
        Atmospheric pressure used to convert absolute to gauge pressure (psi).
    tide_window : int
        Samples for the low-pass (tide) moving mean. Defaults to 1 h.
    wave_window : int
        Samples for the high-pass (wave) moving mean. Defaults to 10 min.
    """
    data = loadmat(mat_path)
    if "pclip" not in data or "tclip" not in data:
        raise KeyError("MAT file must contain 'pclip' and 'tclip'")
    pclip = np.asarray(data["pclip"], dtype=np.float64)
    tclip = np.asarray(data["tclip"], dtype=np.float64)
    if pclip.shape != tclip.shape:
        raise ValueError("pclip and tclip must share the same shape")

    nsamp, nburst = pclip.shape
    gap = 20  # samples
    block_len = nsamp + gap

    # convert to gauge Pascals and add NaN-filled gap
    pgauge_psi = pclip - patm_psi
    pgauge_pa = pgauge_psi * PSI_TO_PA
    stacked = np.full((block_len, nburst), np.nan)
    stacked[:nsamp, :] = pgauge_pa
    pressure_pa = stacked.reshape(-1, order="F")[:-gap]

    # water depth inferred from hydrostatic balance
    depth_stack = np.full((block_len, nburst), np.nan)
    depth_stack[:nsamp, :] = pgauge_pa / (RHO_SEAWATER * G)
    depth_m = depth_stack.reshape(-1, order="F")[:-gap]

    # build continuous time vector (1 Hz sampling)
    dt_seconds = 1.0
    n_total = pressure_pa.size
    time_seconds = np.arange(n_total, dtype=np.float64) * dt_seconds
    time_days = tclip[0, 0] + time_seconds / 86400.0
    time = matlab_datenum_to_datetime64(time_days)

    # slow-varying tide (ignore NaNs) and high-pass pressure component
    depth_lp = moving_average_nan(depth_m, tide_window)
    pressure_trend = moving_average_nan(pressure_pa, wave_window)
    pressure_hp = pressure_pa - pressure_trend

    return SensorSeries(
        time=time,
        seconds=time_seconds,
        pressure_pa=pressure_pa,
        depth_m=depth_m,
        pressure_hp=pressure_hp,
        depth_lp=depth_lp,
    )


def compute_sliding_spectra(series: SensorSeries, win_len: int = 2048, step: int = 512) -> SpectralSummary:
    """Welch-like spectral analysis in sliding windows."""
    p_hp = nan_interp(series.pressure_hp)
    depth_lp = nan_interp(series.depth_lp)
    n = p_hp.size
    if n < win_len:
        raise ValueError("time series shorter than window length")
    starts = np.arange(0, n - win_len + 1, step, dtype=int)
    hann = np.hanning(win_len)
    hann *= math.sqrt(win_len / np.sum(hann ** 2))  # energy conserving

    freqs = np.fft.rfftfreq(win_len, d=1.0)
    omega = 2.0 * np.pi * freqs

    spectra = []
    hs_total = []
    hs_ss = []
    hs_ig = []
    tp_ss = []
    t_mid = []

    mask_total = freqs > 0.0
    mask_ss = (freqs >= 0.05) & (freqs <= 0.33)
    mask_ig = (freqs >= 0.004) & (freqs <= 0.04)

    for start in starts:
        stop = start + win_len
        seg = p_hp[start:stop]
        depth_seg = depth_lp[start:stop]
        if np.any(~np.isfinite(seg)):
            continue
        seg = seg - np.mean(seg)
        seg *= hann
        fft_vals = np.fft.rfft(seg)

        # one-sided PSD of pressure (Pa^2/Hz)
        scale = (1.0 / win_len) ** 2
        Spp = (2.0 * scale) * (np.abs(fft_vals) ** 2)
        Spp[0] = scale * (np.abs(fft_vals[0]) ** 2)
        if win_len % 2 == 0:
            Spp[-1] = scale * (np.abs(fft_vals[-1]) ** 2)
        Spp *= 1.0  # dt = 1 s already accounted for

        h_eff = float(np.nanmean(depth_seg))
        k = wavenumber(omega, h_eff)
        transfer = np.cosh(k * h_eff) / (RHO_SEAWATER * G)
        Seta = (transfer ** 2) * Spp
        spectra.append(Seta)

        def m0(mask: np.ndarray) -> float:
            if not np.any(mask):
                return 0.0
            return float(np.trapz(Seta[mask], freqs[mask]))

        m0_total = m0(mask_total)
        m0_ss = m0(mask_ss)
        m0_ig = m0(mask_ig)

        hs_total.append(4.0 * math.sqrt(max(m0_total, 0.0)))
        hs_ss.append(4.0 * math.sqrt(max(m0_ss, 0.0)))
        hs_ig.append(4.0 * math.sqrt(max(m0_ig, 0.0)))

        if np.any(mask_ss):
            ss_spec = Seta[mask_ss]
            if np.all(np.isfinite(ss_spec)) and np.nanmax(ss_spec) > 0.0:
                fp = freqs[mask_ss][np.nanargmax(ss_spec)]
                tp_ss.append(1.0 / fp if fp > 0.0 else np.nan)
            else:
                tp_ss.append(np.nan)
        else:
            tp_ss.append(np.nan)

        t_idx = start + win_len // 2
        t_mid.append(series.time[t_idx])

    if not spectra:
        raise RuntimeError("no valid segments found for spectral analysis")

    return SpectralSummary(
        freqs=freqs,
        surface_spectra=np.vstack(spectra),
        time_centers=np.array(t_mid),
        hs_total=np.array(hs_total),
        hs_sea_swell=np.array(hs_ss),
        hs_infragravity=np.array(hs_ig),
        tp_sea_swell=np.array(tp_ss),
    )


def compute_spectrogram_summary(
    series: SensorSeries,
    *,
    nperseg: int = 8192,
    noverlap: int | None = None,
    window: str = "hann",
    detrend: str | None = "constant",
) -> SpectralSummary:
    """
    Time-evolving autospectra using scipy.signal.spectrogram.

    Longer segment lengths (``nperseg``) push the frequency resolution to lower
    bands while retaining a rolling estimate of the sea/swell energy.
    """
    p_hp = nan_interp(series.pressure_hp)
    depth_lp = nan_interp(series.depth_lp)
    n_samples = p_hp.size
    if nperseg > n_samples:
        raise ValueError("nperseg exceeds available samples")
    if noverlap is None:
        noverlap = nperseg // 2
    if not 0 <= noverlap < nperseg:
        raise ValueError("noverlap must satisfy 0 <= noverlap < nperseg")

    fs = 1.0  # Hz
    freqs, times_sec, Spp = spectrogram(
        p_hp,
        fs=fs,
        window=window,
        nperseg=nperseg,
        noverlap=noverlap,
        detrend=detrend,
        scaling="density",
        mode="psd",
        boundary=None,
        padded=False,
    )
    if Spp.size == 0:
        raise RuntimeError("spectrogram returned empty PSD array")

    omega = 2.0 * np.pi * freqs
    seconds = np.asarray(series.seconds, dtype=np.float64)
    time_int = series.time.astype("datetime64[ns]").astype("int64")

    depth_interp = np.interp(times_sec, seconds, depth_lp)
    time_ns = np.interp(times_sec, seconds, time_int)
    time_centers = np.round(time_ns).astype("int64").astype("datetime64[ns]")

    Setas = np.empty_like(Spp)
    mask_total = freqs > 0.0
    mask_ss = (freqs >= 0.05) & (freqs <= 0.33)
    mask_ig = (freqs >= 0.004) & (freqs <= 0.04)

    hs_total = []
    hs_ss = []
    hs_ig = []
    tp_ss = []

    for j, h_eff in enumerate(depth_interp):
        k = wavenumber(omega, float(h_eff))
        transfer = np.cosh(k * h_eff) / (RHO_SEAWATER * G)
        Setas[:, j] = (transfer ** 2) * Spp[:, j]

        spec = Setas[:, j]

        def m0(mask: np.ndarray) -> float:
            if not np.any(mask):
                return 0.0
            return float(np.trapz(spec[mask], freqs[mask]))

        m0_total = m0(mask_total)
        m0_ss = m0(mask_ss)
        m0_ig = m0(mask_ig)

        hs_total.append(4.0 * math.sqrt(max(m0_total, 0.0)))
        hs_ss.append(4.0 * math.sqrt(max(m0_ss, 0.0)))
        hs_ig.append(4.0 * math.sqrt(max(m0_ig, 0.0)))

        if np.any(mask_ss):
            ss_spec = spec[mask_ss]
            if np.all(np.isfinite(ss_spec)) and np.nanmax(ss_spec) > 0.0:
                fp = freqs[mask_ss][np.nanargmax(ss_spec)]
                tp_ss.append(1.0 / fp if fp > 0.0 else np.nan)
            else:
                tp_ss.append(np.nan)
        else:
            tp_ss.append(np.nan)

    return SpectralSummary(
        freqs=freqs,
        surface_spectra=Setas.T,
        time_centers=time_centers,
        hs_total=np.array(hs_total),
        hs_sea_swell=np.array(hs_ss),
        hs_infragravity=np.array(hs_ig),
        tp_sea_swell=np.array(tp_ss),
    )


def to_datetime(objects: Iterable[np.datetime64]) -> list[datetime]:
    """Convert numpy datetime64 array to Python datetime list for plotting."""
    arr = np.asarray(list(objects), dtype="datetime64[ns]")
    return arr.astype("datetime64[ms]").astype(object).tolist()


def plot_results(series: SensorSeries, summary: SpectralSummary, out_dir: Path) -> None:
    """Create figures similar to PILOTmokuleia sensor summaries."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # Hs time series
    time_mid = to_datetime(summary.time_centers)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(time_mid, summary.hs_total, label="Hs (total)")
    ax.plot(time_mid, summary.hs_sea_swell, label="Hs (sea/swell)")
    ax.plot(time_mid, summary.hs_infragravity, label="Hs (IG)")
    ax.set_ylabel("Significant wave height [m]")
    ax.set_xlabel("Time")
    ax.set_title("Sensor 1 significant wave height")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    fig.tight_layout()
    fig.savefig(out_dir / "sensor1_hs_timeseries.png", dpi=200)
    plt.close(fig)

    # peak period
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(time_mid, summary.tp_sea_swell, color="#d95f02")
    ax.set_ylabel("Tp (sea/swell) [s]")
    ax.set_xlabel("Time")
    ax.set_title("Sea/swell peak period")
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    fig.tight_layout()
    fig.savefig(out_dir / "sensor1_tp_timeseries.png", dpi=200)
    plt.close(fig)

    # representative spectrum (median over time)
    S_eta_median = np.nanmedian(summary.surface_spectra, axis=0)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(summary.freqs[1:], S_eta_median[1:], lw=1.8)
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Surface elevation spectrum [m^2/Hz]")
    ax.set_title("Median surface elevation spectrum (sensor 1)")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "sensor1_median_spectrum.png", dpi=200)
    plt.close(fig)

    # raw gauge pressure vs time for context
    time_all = to_datetime(series.time)
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(time_all, series.pressure_pa / (RHO_SEAWATER * G), color="#1b9e77", lw=0.6)
    ax.set_ylabel("Water depth [m]")
    ax.set_xlabel("Time")
    ax.set_title("Sensor 1 hydrostatic water depth")
    ax.grid(True, alpha=0.2)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    fig.tight_layout()
    fig.savefig(out_dir / "sensor1_depth_series.png", dpi=200)
    plt.close(fig)


def main() -> None:
    base_dir = Path(__file__).resolve().parent
    mat_path = base_dir / "MoA18411.mat"
    series = load_sensor_series(mat_path)
    summary = compute_sliding_spectra(series)
    plot_results(series, summary, base_dir / "figures")


if __name__ == "__main__":
    main()
