# Calculating wave properties

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd


G = 9.81
RHO_SEAWATER = 1025.0


def _tspec_to_seconds_and_centers(t_spec, t1):
    """Normalize t_spec input.

    Accepts either:
    - numeric seconds offsets (array-like), or
    - datetime-like array of absolute times (datetime64 or pandas.Timestamp).

    Returns (t_spec_seconds, time_centers, t0) where:
    - t_spec_seconds: numpy float array of seconds since t0
    - time_centers: numpy datetime64[ns] array of absolute times for each window
    - t0: numpy datetime64[ns] scalar representing the reference start time
    """
    # make numpy array without forcing numeric dtype yet
    arr = np.asarray(t_spec)
    t0 = np.array(pd.to_datetime(t1[0]), dtype="datetime64[ns]")
    if np.issubdtype(arr.dtype, np.datetime64):
        # t_spec already absolute datetimes -> compute seconds since t0
        t_spec_dt = np.asarray(pd.to_datetime(arr), dtype="datetime64[ns]")
        time_centers = t_spec_dt
        t_spec_seconds = ((t_spec_dt - t0) / np.timedelta64(1, "s")).astype(float)
    else:
        # assume numeric seconds offsets
        t_spec_seconds = np.asarray(arr, dtype=float)
        time_offsets = (t_spec_seconds * 1e9).astype("timedelta64[ns]")
        time_centers = t0 + time_offsets
    return t_spec_seconds, time_centers, t0

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


def Spp_to_Seta(
    Spp: np.ndarray,
    freqs: np.ndarray,
    t_spec: np.ndarray,
    t1: np.ndarray,
    h1: np.ndarray,
    *,
    depth_interp: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert a pressure spectrogram into surface-elevation spectra.

    Parameters
    ----------
    Spp : np.ndarray
        Pressure power spectral density (Pa^2/Hz), shape (n_freqs, n_windows).
    freqs : np.ndarray
        Frequency vector [Hz] corresponding to the rows of ``Spp``.
    t_spec : np.ndarray
        Spectrogram time centers in seconds since the start of the record.
    t1 : np.ndarray
        Native time array (datetime64 or compatible) for the 1 Hz pressure record.
    h1 : np.ndarray
        Hydrostatic depth series (may include NaNs).
    depth_interp : np.ndarray, optional
        Pre-interpolated depth at 1 Hz. If omitted, a linear interpolation
        (limit 20 samples) is performed internally.

    Returns
    -------
    Seta : np.ndarray
        Surface-elevation spectra (m^2/Hz), same shape as ``Spp``.
    time_centers : np.ndarray
        Datetime64[ns] array for the center of each spectrogram window.
    depth_at_centers : np.ndarray
        Depth used for each window (meters).
    """
    Spp = np.asarray(Spp, dtype=np.float64)
    freqs = np.asarray(freqs, dtype=np.float64)
    # accept either seconds offsets or datetime-like t_spec
    t_spec_seconds, time_centers, t0 = _tspec_to_seconds_and_centers(t_spec, t1)
    if Spp.ndim != 2:
        raise ValueError("Spp must be 2-D (n_freqs, n_windows)")
    if freqs.ndim != 1 or freqs.size != Spp.shape[0]:
        raise ValueError("freqs must be 1-D and match the first dimension of Spp")

    time_index = pd.to_datetime(t1)
    if depth_interp is None:
        depth_interp = (
            pd.Series(h1, index=time_index)
            .interpolate(method="linear", limit=20, limit_direction="both")
            .to_numpy()
        )
    else:
        depth_interp = np.asarray(depth_interp, dtype=np.float64)
        if depth_interp.shape != h1.shape:
            raise ValueError("depth_interp must match the shape of h1")

    seconds_full = (time_index.to_numpy() - t0) / np.timedelta64(1, "s")
    depth_at_centers = np.interp(t_spec_seconds, seconds_full, depth_interp)

    omega = 2.0 * np.pi * freqs
    Seta = np.empty_like(Spp, dtype=np.float64)
    for col, depth_val in enumerate(depth_at_centers):
        k = wavenumber(omega, float(depth_val))
        transfer = np.cosh(k * depth_val) / (RHO_SEAWATER * G)
        Seta[:, col] = (transfer ** 2) * Spp[:, col]

    return Seta, time_centers, depth_at_centers


def sig_wave_height(
    Seta: np.ndarray,
    freqs: np.ndarray,
    t_spec: np.ndarray,
    t1: np.ndarray,
    h1: np.ndarray,
    bands,
    *,
    depth_interp: np.ndarray | None = None,
    dataviz: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute bulk wave statistics from a surface elevation spectrogram.

    Parameters
    ----------
    Seta : np.ndarray
        Surface elevation spectra (m^2/Hz) with shape (n_freqs, n_windows).
    freqs : np.ndarray
        Frequency vector corresponding to rows of ``Spp`` [Hz].
    t_spec : np.ndarray
        Spectrogram time centers in seconds since the start of the record.
    t1 : np.ndarray
        Native datetime64 array for the 1 Hz series.
    h1 : np.ndarray
        Raw hydrostatic depth series (may contain NaNs).
    bands : list
        desired inputs of frequency bands to be integrated     
    depth_interp : np.ndarray, optional
        Pre-interpolated depth series matching ``h1`` (1 Hz). If not supplied, it is built via linear interpolation.
    dataviz : bool, optional
        When True, generate a quick-look plot of the Hs time series.

    Returns
    -------
    hs_tot, hs_ig, hs_ss, hs_input, tp_ss, time_centers : tuple of np.ndarray
        Bulk metrics for each spectrogram window and the corresponding datetime64 centers.
    """
    RHO = 1025.0

    if depth_interp is None:
        depth_interp = (
            pd.Series(h1, index=pd.to_datetime(t1))
            .interpolate(method="linear", limit=20, limit_direction="both")
            .to_numpy()
        )

    # spectrogram centers as datetime64 and seconds offsets for interpolation
    t_spec_seconds, time_centers, t0 = _tspec_to_seconds_and_centers(t_spec, t1)
    seconds_full = (t1 - t1[0]) / np.timedelta64(1, "s")
    depth_at_centers = np.interp(t_spec_seconds, seconds_full, depth_interp)

    mask_total = freqs > 0.0
    mask_ss = (freqs >= 0.05) & (freqs <= 0.33)
    mask_ig = (freqs >= 0.004) & (freqs <= 0.04)
    mask_input = (freqs >= bands[0]) & (freqs <= bands[1])

    omega = 2.0 * np.pi * freqs

    hs_tot: list[float] = []
    hs_ss: list[float] = []
    hs_ig: list[float] = []
    hs_bands: list[float] = []
    tp_ss: list[float] = []


    for col, depth_val in enumerate(depth_at_centers):
        seta = Seta[:, col]

        def spectral_moment(mask: np.ndarray) -> float:
            if not np.any(mask):
                return 0.0
            return float(np.trapz(seta[mask], freqs[mask]))

        m0_total = spectral_moment(mask_total)
        m0_sea_swell = spectral_moment(mask_ss)
        m0_ig = spectral_moment(mask_ig)
        m0_input = spectral_moment(mask_input)

        hs_tot.append(4.0 * np.sqrt(max(m0_total, 0.0)))
        hs_ss.append(4.0 * np.sqrt(max(m0_sea_swell, 0.0)))
        hs_ig.append(4.0 * np.sqrt(max(m0_ig, 0.0)))
        hs_bands.append(4.0 * np.sqrt(max(m0_input, 0.0)))

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
    hs_bands = np.asarray(hs_bands)
    tp_ss = np.asarray(tp_ss)

    if dataviz:
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(time_centers, hs_tot, label="Hs (total)")
        ax.plot(time_centers, hs_ig, label="Hs (IG): 0.004–0.04 Hz")
        ax.plot(time_centers, hs_ss, label="Hs (sea/swell): 0.05–0.33 Hz")
        ax.set_ylabel("Significant wave height [m]")
        ax.set_xlabel("Date (local time)")
        ax.set_title("Zero-th Moment Derived significant wave height")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        fig.autofmt_xdate()
        plt.show()

    return hs_tot, hs_ig, hs_ss, hs_bands, tp_ss, time_centers


from utide import solve, reconstruct



def detide(
    p: np.ndarray,
    t: np.ndarray,
    LAT, 
    avg = '15min'
    ):

    ''' Removes the tidal signal given a pressure/surface elevation/current array

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

    #  5min means to help with computation of tides
    p_30 = p_interp.resample(avg).mean().dropna()
    idx_30s = p_30.index

    # tideal computation on 30s grid
    coef = solve(
        idx_30s.to_numpy(),
        p_30.to_numpy(),
        lat=LAT,
        method="ols",
        conf_int="linear", # 'none' skips uncertainty calc, 'linear' uses quick lin. estimate, 'MC' run monte-carlo (mad spendy)
        nodal=True,
        trend=True,
        verbose=False)

    # # Reconstruct on the 1 Hz grid (chunking to avoid crashes)
    chunks = np.array_split(idx_1hz, 8)  # adjust slice count for memory
    tide_parts = []
    t0 = idx_1hz[0]
    for chunk in chunks:
        # chunk is a DatetimeIndex
        recon_chunk = reconstruct(chunk.to_numpy(), coef)
        tide_parts.append(pd.Series(recon_chunk.h.astype("float32"), index=chunk))
    tide_full = pd.concat(tide_parts).sort_index()
    p_wave = (p_interp - tide_full) # water level chane due to wave effects only
    p_prime = p_wave.resample('1s').mean().dropna()
    p_full = p_interp

    # # keep it on the averaged grid
    # chunks = np.array_split(idx_30s, 8)  # adjust slice count for memory
    # tide_parts = []
    # t0 = idx_30s[0]
    # for chunk in chunks:
    #     # chunk is a DatetimeIndex
    #     recon_chunk = reconstruct(chunk.to_numpy(), coef)
    #     tide_parts.append(pd.Series(recon_chunk.h.astype("float32"), index=chunk))
    # tide_full = pd.concat(tide_parts).sort_index()
    # tide_full = tide_full.reindex(p_interp.index).interpolate("time")

    # p_prime = (p_interp - tide_full) # water level chane due to wave effects only
    return tide_full, p_prime, p_interp



