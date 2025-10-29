from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat

# --- constants ---
PSI_TO_PA   = 6894.757293168     # 1 psi  -> Pa
DBAR_TO_PA  = 1e4                # 1 dbar -> Pa
MATLAB_EPOCH_OFFSET = 719529.0   # days between 0000-01-01 and 1970-01-01
SECONDS_PER_DAY     = 86400.0

def matlab_datenum_to_datetime64(tnums):
    """
    Convert MATLAB datenum(s) to numpy datetime64[ns].
    """
    tnums = np.asarray(tnums, dtype=np.float64)
    delta_days = tnums - MATLAB_EPOCH_OFFSET
    ns_since_unix = np.round(delta_days * SECONDS_PER_DAY * 1e9).astype("int64")
    return np.datetime64("1970-01-01") + ns_since_unix.astype("timedelta64[ns]")

def _normalize_burst_starts(tclip, nsamp, nburst):
    """
    Return a 1D array of length nburst with MATLAB datenum start times per burst.
    Accepts:
      - scalar (same start for all bursts)
      - (nburst,), (1, nburst), (nburst, 1)
      - (nsamp, nburst)  -> uses the first sample of each burst as its start
    """
    tc = np.asarray(tclip, dtype=float)

    # If nothing usable, bail early with a clear error
    if tc.size == 0:
        raise ValueError("tclip is empty; cannot construct time index.")

    # Scalar -> broadcast to all bursts
    if tc.size == 1:
        return np.full(nburst, float(np.squeeze(tc)))

    # 1-D vector length nburst
    if tc.ndim == 1 and tc.shape[0] == nburst:
        return tc.astype(float)

    # 2-D forms with a single row/col that flatten to length nburst
    if tc.ndim == 2 and 1 in tc.shape and max(tc.shape) == nburst:
        return tc.reshape(nburst).astype(float)

    # Full per-sample times: (nsamp, nburst) -> take first sample as burst start
    if tc.ndim == 2 and tc.shape == (nsamp, nburst):
        return tc[0, :].astype(float)

    raise ValueError(
        f"Unrecognized tclip shape {tc.shape}; expected scalar, "
        f"(nburst,), (1,nburst), (nburst,1), or (nsamp,nburst)."
    )

def loadPressureData(
    mat_path: Path,
    *,
    sample_rate_hz: float = 1.0,      # samples per second within each burst
    gap_seconds: float = 20.0,        # seconds of gap to insert between bursts
    units: str = "psi",               # "psi", "dbar", or "pa"
    is_gauge: bool = False,           # True if pclip already has atmosphere removed
    patm_psi: float = 14.6959,        # local atmospheric pressure [psi] for absolute psi
    patm_dbar: float = 10.1325,       # local atmospheric pressure [dbar] for absolute dbar
    patm_pa: float = 101325.0,        # local atmospheric pressure [Pa]  for absolute Pa
    rho: float = 1025.0,              # seawater density [kg/m^3]
    gravity: float = 9.81,            # gravitational acceleration [m/s^2]
    clip_negative_gauge: bool = True  # clip tiny negative gauge pressures to 0 Pa
) -> pd.DataFrame:
    """
    Load bursty pressure data from a MATLAB file and return a continuous time series.

    Expects .mat to contain:
      - pclip : (nsamp, nburst) pressure in 'units'
      - tclip : scalar, (nburst,), (1,nburst), (nburst,1), or (nsamp,nburst)
                representing burst start time(s) or full per-sample times (MATLAB datenum)

    Returns a DataFrame with datetime index and columns:
      - pressure_raw         : raw pressure in original 'units'
      - pressure_gauge_pa    : gauge pressure in Pa (>=0 if clipping enabled)
      - h                    : hydrostatic depth (m) = pressure_gauge_pa / (rho*g)
    """
    # --- load from .mat ---
    data = loadmat(mat_path)
    if "pclip" not in data or "tclip" not in data:
        raise KeyError("Input .mat must contain 'pclip' and 'tclip' variables.")

    pclip = np.asarray(data["pclip"], dtype=float)   # (nsamp, nburst)
    tclip = np.asarray(data["tclip"], dtype=float)   # various shapes
    if pclip.ndim != 2:
        raise ValueError(f"pclip must be 2-D (nsamp, nburst); got shape {pclip.shape}")

    nsamp, nburst = pclip.shape

    # --- choose atmospheric pressure in incoming units and convert to Pa ---
    u = units.lower()
    if u == "psi":
        patm_in = 0.0 if is_gauge else float(patm_psi)
        p_gauge_pa = (pclip - patm_in) * PSI_TO_PA
    elif u == "dbar":
        patm_in = 0.0 if is_gauge else float(patm_dbar)
        p_gauge_pa = (pclip - patm_in) * DBAR_TO_PA
    elif u == "pa":
        patm_in = 0.0 if is_gauge else float(patm_pa)
        p_gauge_pa = (pclip - patm_in)
    else:
        raise ValueError("units must be one of {'psi','dbar','pa'}")

    if clip_negative_gauge:
        p_gauge_pa = np.maximum(p_gauge_pa, 0.0)

    # --- hydrostatic depth (water column above the transducer) ---
    depth_m_bursty = p_gauge_pa / (rho * gravity)  # (nsamp, nburst)

    # --- stitch bursts with NaN gaps (column-major, to preserve within-burst order) ---
    gap_samples = int(round(gap_seconds * sample_rate_hz))
    block_len = nsamp + (gap_samples if gap_samples > 0 else 0)

    # Prepare containers (raw in original units; gauge in Pa; depth in m)
    stacked_raw = np.full((block_len, nburst), np.nan)
    stacked_gpa = np.full((block_len, nburst), np.nan)
    stacked_h   = np.full((block_len, nburst), np.nan)

    stacked_raw[:nsamp, :] = pclip
    stacked_gpa[:nsamp, :] = p_gauge_pa
    stacked_h[:nsamp,   :] = depth_m_bursty

    # Flatten in Fortran order to match burst sequencing
    pressure_raw      = stacked_raw.reshape(-1, order="F")
    pressure_gauge_pa = stacked_gpa.reshape(-1, order="F")
    depth_m           = stacked_h.reshape(-1, order="F")

    # Drop the trailing gap after the last burst
    if gap_samples > 0:
        pressure_raw      = pressure_raw[:-gap_samples]
        pressure_gauge_pa = pressure_gauge_pa[:-gap_samples]
        depth_m           = depth_m[:-gap_samples]

    # --- build the time index robustly (per-burst, respecting gaps) ---
    burst_starts = _normalize_burst_starts(tclip, nsamp, nburst)  # length nburst, MATLAB datenums
    # Time within a burst as MATLAB datenum increments
    within_burst = (np.arange(nsamp, dtype=float) / SECONDS_PER_DAY)

    # Preallocate stitched MATLAB datenums with NaN gaps to mirror data stitching
    tstack = np.full((block_len, nburst), np.nan, dtype=float)
    for j in range(nburst):
        tstart = float(burst_starts[j])
        tstack[:nsamp, j] = tstart + within_burst

    # Flatten time the same way; drop trailing gap
    matlab_time_series = tstack.reshape(-1, order="F")
    if gap_samples > 0:
        matlab_time_series = matlab_time_series[:-gap_samples]

    # Convert MATLAB datenums -> datetime64[ns]
    dt = matlab_datenum_to_datetime64(matlab_time_series)

    # --- assemble tidy DataFrame ---
    df = pd.DataFrame(
        {
            "p_raw": pressure_raw,            # in 'units'
            "p": pressure_gauge_pa,  # Pa
            "h": depth_m,                            # m
        },
        index=pd.to_datetime(dt)
    )

    # --- attach metadata for provenance ---
    df.attrs["sample_rate_hz"] = float(sample_rate_hz)
    df.attrs["rho"] = float(rho)
    df.attrs["gravity"] = float(gravity)
    df.attrs["units"] = u
    df.attrs["is_gauge"] = bool(is_gauge)
    if u == "psi":
        df.attrs["patm_used_psi"] = 0.0 if is_gauge else float(patm_psi)
    elif u == "dbar":
        df.attrs["patm_used_dbar"] = 0.0 if is_gauge else float(patm_dbar)
    else:
        df.attrs["patm_used_pa"] = 0.0 if is_gauge else float(patm_pa)

    return df
