# Calculating wave properties

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd


G = 9.81

def wavenumber(omega: np.ndarray, depth: float, tol: float = 1e-12, max_iter: int = 64) -> np.ndarray:
    """Solve the linear dispersion relation for k(ω) using Newton's method"""
    G = 9.81
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


def sig_wave_height(
    Spp: np.ndarray,
    freqs: np.ndarray,
    t_spec: np.ndarray,
    t1: np.ndarray,
    h1: np.ndarray,
    *,
    depth_interp: np.ndarray | None = None,
    dataviz: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute bulk wave statistics from a pressure spectrogram.

    Parameters
    ----------
    Spp : np.ndarray
        Pressure spectra (Pa^2/Hz) from ``scipy.signal.spectrogram`` with shape (n_freqs, n_windows).
    freqs : np.ndarray
        Frequency vector corresponding to rows of ``Spp`` [Hz].
    t_spec : np.ndarray
        Spectrogram time centers in seconds since the start of the record.
    t1 : np.ndarray
        Native datetime64 array for the 1 Hz series.
    h1 : np.ndarray
        Raw hydrostatic depth series (may contain NaNs).
    depth_interp : np.ndarray, optional
        Pre-interpolated depth series matching ``h1`` (1 Hz). If not supplied, it is built via linear interpolation.
    dataviz : bool, optional
        When True, generate a quick-look plot of the Hs time series.

    Returns
    -------
    hs_tot, hs_ig, hs_ss, tp_ss, time_centers : tuple of np.ndarray
        Bulk metrics for each spectrogram window and the corresponding datetime64 centers.
    """
    RHO = 1025.0

    if depth_interp is None:
        depth_interp = (
            pd.Series(h1, index=pd.to_datetime(t1))
            .interpolate(method="linear", limit=20, limit_direction="both")
            .to_numpy()
        )

    # spectrogram centers as datetime64
    t0 = np.array(t1[0], dtype="datetime64[ns]")
    time_offsets = (t_spec * 1e9).astype("timedelta64[ns]")
    time_centers = t0 + time_offsets

    seconds_full = (t1 - t1[0]) / np.timedelta64(1, "s")
    depth_at_centers = np.interp(t_spec, seconds_full, depth_interp)

    mask_total = freqs > 0.0
    mask_ss = (freqs >= 0.05) & (freqs <= 0.33)
    mask_ig = (freqs >= 0.004) & (freqs <= 0.04)

    omega = 2.0 * np.pi * freqs

    hs_tot: list[float] = []
    hs_ss: list[float] = []
    hs_ig: list[float] = []
    tp_ss: list[float] = []
    Seta = np.empty_like(Spp)

    for col, depth_val in enumerate(depth_at_centers):
        k = wavenumber(omega, float(depth_val))
        transfer = np.cosh(k * depth_val) / (RHO * G)
        seta = (transfer ** 2) * Spp[:, col]
        Seta[:, col] = seta

        def spectral_moment(mask: np.ndarray) -> float:
            if not np.any(mask):
                return 0.0
            return float(np.trapz(seta[mask], freqs[mask]))

        m0_total = spectral_moment(mask_total)
        m0_sea_swell = spectral_moment(mask_ss)
        m0_ig = spectral_moment(mask_ig)

        hs_tot.append(4.0 * np.sqrt(max(m0_total, 0.0)))
        hs_ss.append(4.0 * np.sqrt(max(m0_sea_swell, 0.0)))
        hs_ig.append(4.0 * np.sqrt(max(m0_ig, 0.0)))

        if np.any(mask_ss):
            ss_slice = seta[mask_ss]
            if np.all(np.isfinite(ss_slice)) and np.nanmax(ss_slice) > 0.0:
                fp = freqs[mask_ss][np.nanargmax(ss_slice)]
                tp_ss.append(1.0 / fp if fp > 0.0 else np.nan)
            else:
                tp_ss.append(np.nan)
        else:
            tp_ss.append(np.nan)

    hs_tot = np.asarray(hs_tot)
    hs_ss = np.asarray(hs_ss)
    hs_ig = np.asarray(hs_ig)
    tp_ss = np.asarray(tp_ss)

    if dataviz:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(time_centers, hs_tot, label="Hs (total)")
        ax.plot(time_centers, hs_ig, label="Hs (IG): 0.004–0.04 Hz")
        ax.plot(time_centers, hs_ss, label="Hs (sea/swell): 0.05–0.33 Hz")
        ax.set_ylabel("Significant wave height [m]")
        ax.set_xlabel("Date (local time)")
        ax.set_title("Spectrogram-derived significant wave height")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        fig.autofmt_xdate()
        plt.show()

    return hs_tot, hs_ig, hs_ss, tp_ss, time_centers







# def compute_spectrogram_summary(
#     series: SensorSeries,
#     *,
#     nperseg: int = 8192,
#     noverlap: int | None = None,
#     window: str = "hann",
#     detrend: str | None = "constant",
# ) -> SpectralSummary:
#     """
#     Time-evolving autospectra using scipy.signal.spectrogram.
#     Longer segment lengths (``nperseg``) push the frequency resolution to lower
#     bands while retaining a rolling estimate of the sea/swell energy.
#     """
#     p_hp = nan_interp(series.pressure_hp)
#     depth_lp = nan_interp(series.depth_lp)
#     n_samples = p_hp.size
#     if nperseg > n_samples:
#         raise ValueError("nperseg exceeds available samples")
#     if noverlap is None:
#         noverlap = nperseg // 2
#     if not 0 <= noverlap < nperseg:
#         raise ValueError("noverlap must satisfy 0 <= noverlap < nperseg")


#     omega = 2.0 * np.pi * freqs
#     seconds = np.asarray(series.seconds, dtype=np.float64)
#     time_int = series.time.astype("datetime64[ns]").astype("int64")

#     depth_interp = np.interp(times_sec, seconds, depth_lp)
#     time_ns = np.interp(times_sec, seconds, time_int)
#     time_centers = np.round(time_ns).astype("int64").astype("datetime64[ns]")

#     Setas = np.empty_like(Spp)
#     mask_total = freqs > 0.0
#     mask_ss = (freqs >= 0.05) & (freqs <= 0.33)
#     mask_ig = (freqs >= 0.004) & (freqs <= 0.04)

#     hs_total = []
#     hs_ss = []
#     hs_ig = []
#     tp_ss = []

#     for j, h_eff in enumerate(depth_interp):
#         k = wavenumber(omega, float(h_eff))
#         transfer = np.cosh(k * h_eff) / (RHO_SEAWATER * G)
#         Setas[:, j] = (transfer ** 2) * Spp[:, j]

#         spec = Setas[:, j]

#         def m0(mask: np.ndarray) -> float:
#             if not np.any(mask):
#                 return 0.0
#             return float(np.trapz(spec[mask], freqs[mask]))

#         m0_total = m0(mask_total)
#         m0_ss = m0(mask_ss)
#         m0_ig = m0(mask_ig)

#         hs_total.append(4.0 * math.sqrt(max(m0_total, 0.0)))
#         hs_ss.append(4.0 * math.sqrt(max(m0_ss, 0.0)))
#         hs_ig.append(4.0 * math.sqrt(max(m0_ig, 0.0)))

#         if np.any(mask_ss):
#             ss_spec = spec[mask_ss]
#             if np.all(np.isfinite(ss_spec)) and np.nanmax(ss_spec) > 0.0:
#                 fp = freqs[mask_ss][np.nanargmax(ss_spec)]
#                 tp_ss.append(1.0 / fp if fp > 0.0 else np.nan)
#             else:
#                 tp_ss.append(np.nan)
#         else:
#             tp_ss.append(np.nan)

#     return SpectralSummary(
#         freqs=freqs,
#         surface_spectra=Setas.T,
#         time_centers=time_centers,
#         hs_total=np.array(hs_total),
#         hs_sea_swell=np.array(hs_ss),
#         hs_infragravity=np.array(hs_ig),
#         tp_sea_swell=np.array(tp_ss),
#     )
3
