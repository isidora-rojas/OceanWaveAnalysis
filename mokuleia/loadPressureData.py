from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import loadmat

# constants
PSI_TO_PA = 6894.757293168
MATLAB_EPOCH_OFFSET = 719529.0
SECONDS_PER_DAY = 86400.0

def matlab_datenum_to_datetime64(tnums):
    tnums = np.asarray(tnums, dtype=np.float64)
    delta_days = tnums - MATLAB_EPOCH_OFFSET
    ns_since_unix = np.round(delta_days * SECONDS_PER_DAY * 1e9).astype("int64")
    return np.datetime64("1970-01-01") + ns_since_unix.astype("timedelta64[ns]")

def loadPressureData(
    mat_path: Path,
    *,
    sample_rate_hz: float = 1.0,
    gap_seconds: float = 20.0,
    patm_psi: float = 14.7,
    rho: float = 1025.0,
    gravity: float = 9.81,
) -> pd.DataFrame:
    """
    Return DataFrame with index=time, and columns:
        pressure_pa, depth_m
    """
    
    data = loadmat(mat_path)
    pclip = np.asarray(data["pclip"], float)
    tclip = np.asarray(data["tclip"], float)

    nsamp, nburst = pclip.shape
    gap_samples = int(round(gap_seconds * sample_rate_hz))
    block_len = nsamp + gap_samples if gap_samples > 0 else nsamp

    pgauge_pa = (pclip - patm_psi) * PSI_TO_PA

    # stitch bursts with gaps as NaN
    stacked = np.full((block_len, nburst), np.nan)
    stacked[:nsamp] = pgauge_pa
    pressure_pa = stacked.reshape(-1, order="F")
    if gap_samples > 0:
        pressure_pa = pressure_pa[:-gap_samples]

    depth_stack = np.full((block_len, nburst), np.nan)
    depth_stack[:nsamp] = pgauge_pa / (rho * gravity)
    depth_m = depth_stack.reshape(-1, order="F")
    if gap_samples > 0:
        depth_m = depth_m[:-gap_samples]

    samples = pressure_pa.size
    seconds = np.arange(samples) / sample_rate_hz
    matlab_time = tclip[0, 0] + seconds / SECONDS_PER_DAY
    dt = matlab_datenum_to_datetime64(matlab_time)

    # Return a tidy DataFrame
    df = pd.DataFrame(
        {"pressure_pa": pressure_pa,
         "h": depth_m},
        index=pd.to_datetime(dt)
    )

    # metadata preserved on df
    df.attrs["sample_rate_hz"] = sample_rate_hz
    df.attrs["rho"] = rho
    df.attrs["gravity"] = gravity

    return df
