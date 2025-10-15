#!/usr/bin/env python3
"""
Utilities for turning Mokuleia pressure sensor burst MAT files into continuous
Python series.

The notebook ``sensor1_wave_analysis.ipynb`` includes an exploratory version of
this workflow under the "Load the raw pressure bursts" section.  This module
collects that logic so you can reuse it for any Mokuleia sensor by pointing at
the appropriate ``MoA*.mat`` file.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict

import numpy as np
from scipy.io import loadmat

# constants
PSI_TO_PA = 6894.757293168  # [Pa/psi]
MATLAB_EPOCH_OFFSET = 719529.0
SECONDS_PER_DAY = 86400.0


@dataclass
class PressureBursts:
    """Flattened bursts with gaps reinserted and converted to SI units."""

    time: np.ndarray  # datetime64[ns]
    seconds: np.ndarray  # elapsed seconds from start
    pressure_pa: np.ndarray  # gauge pressure [Pa]
    depth_m: np.ndarray  # hydrostatic depth [m]
    sample_rate_hz: float

    def to_dict(self) -> Dict[str, np.ndarray]:
        """Compact representation for saving with numpy.savez."""
        data = asdict(self)
        data.pop("sample_rate_hz", None)
        return data


def matlab_datenum_to_datetime64(tnums: np.ndarray) -> np.ndarray:
    """Convert MATLAB datenums to numpy datetime64[ns]."""
    tnums = np.asarray(tnums, dtype=np.float64)
    delta_days = tnums - MATLAB_EPOCH_OFFSET
    ns_since_unix = np.round(delta_days * SECONDS_PER_DAY * 1e9).astype("int64")
    return np.datetime64("1970-01-01T00:00:00") + ns_since_unix.astype("timedelta64[ns]")


def load_pressure_bursts(
    mat_path: Path,
    *,
    sample_rate_hz: float = 1.0,
    gap_seconds: float = 20.0,
    patm_psi: float = 14.7,
    rho: float = 1025.0,
    gravity: float = 9.81,
) -> PressureBursts:
    """
    Load and stitch bursts from a Mokuleia sensor MAT file.

    Parameters
    ----------
    mat_path : Path
        Location of the ``MoA*.mat`` file with ``pclip`` and ``tclip`` variables.
    sample_rate_hz : float, optional
        Native burst sample rate in Hz. Defaults to 1 Hz.
    gap_seconds : float, optional
        Length of the download gap between bursts (seconds). Defaults to 20 s.
    patm_psi : float, optional
        Atmospheric pressure used to convert absolute to gauge (psi).
    rho : float, optional
        Water density [kg/m^3], defaults to seawater ~1025.
    gravity : float, optional
        Gravitational acceleration [m/s^2].

    Returns
    -------
    PressureBursts
        Flattened time series with gaps reinserted as NaNs and converted to SI units.
    """
    data = loadmat(mat_path)
    for key in ("pclip", "tclip"):
        if key not in data:
            raise KeyError(f"MAT file must contain '{key}'")

    pclip = np.asarray(data["pclip"], dtype=np.float64)
    tclip = np.asarray(data["tclip"], dtype=np.float64)
    if pclip.shape != tclip.shape:
        raise ValueError("pclip and tclip must share the same shape")

    nsamp, nburst = pclip.shape
    gap_samples = int(round(gap_seconds * sample_rate_hz))
    block_len = nsamp + gap_samples if gap_samples > 0 else nsamp

    # convert absolute pressure to gauge pressure in Pascals
    pgauge_pa = (pclip - patm_psi) * PSI_TO_PA

    stacked = np.full((block_len, nburst), np.nan, dtype=np.float64)
    stacked[:nsamp, :] = pgauge_pa
    pressure_pa = stacked.reshape(-1, order="F")
    if gap_samples > 0:
        pressure_pa = pressure_pa[:-gap_samples]

    depth_stack = np.full((block_len, nburst), np.nan, dtype=np.float64)
    depth_stack[:nsamp, :] = pgauge_pa / (rho * gravity)
    depth_m = depth_stack.reshape(-1, order="F")
    if gap_samples > 0:
        depth_m = depth_m[:-gap_samples]

    samples = pressure_pa.size
    seconds = np.arange(samples, dtype=np.float64) / sample_rate_hz
    matlab_time = tclip[0, 0] + seconds / SECONDS_PER_DAY
    time = matlab_datenum_to_datetime64(matlab_time)

    return PressureBursts(
        time=time,
        seconds=seconds,
        pressure_pa=pressure_pa,
        depth_m=depth_m,
        sample_rate_hz=sample_rate_hz,
    )


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Load Mokuleia pressure bursts into a continuous Python series."
    )
    parser.add_argument("mat_path", type=Path, help="Path to the MoA*.mat sensor file.")
    parser.add_argument("--sample-rate", type=float, default=1.0, help="Sensor sampling rate in Hz.")
    parser.add_argument("--gap-seconds", type=float, default=20.0, help="Download gap duration in seconds.")
    parser.add_argument("--patm-psi", type=float, default=14.7, help="Reference atmospheric pressure in psi.")
    parser.add_argument("--rho", type=float, default=1025.0, help="Water density in kg/m^3.")
    parser.add_argument("--gravity", type=float, default=9.81, help="Gravitational acceleration in m/s^2.")
    parser.add_argument(
        "--save-npz",
        type=Path,
        help="Optional .npz output path for saving time, pressure, and depth arrays.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    bursts = load_pressure_bursts(
        args.mat_path,
        sample_rate_hz=args.sample_rate,
        gap_seconds=args.gap_seconds,
        patm_psi=args.patm_psi,
        rho=args.rho,
        gravity=args.gravity,
    )

    duration_hours = bursts.seconds[-1] / 3600.0 if bursts.seconds.size else 0.0
    start = bursts.time[0] if bursts.time.size else "N/A"
    end = bursts.time[-1] if bursts.time.size else "N/A"

    print(f"Loaded {bursts.pressure_pa.size} samples at {bursts.sample_rate_hz:.3f} Hz.")
    print(f"Time span: {start} to {end} (~{duration_hours:.2f} h).")
    print(
        "Depth stats (m): "
        f"mean={np.nanmean(bursts.depth_m):.3f}, "
        f"min={np.nanmin(bursts.depth_m):.3f}, "
        f"max={np.nanmax(bursts.depth_m):.3f}"
    )

    if args.save_npz:
        args.save_npz.parent.mkdir(parents=True, exist_ok=True)
        np.savez(args.save_npz, sample_rate_hz=bursts.sample_rate_hz, **bursts.to_dict())
        print(f"Saved arrays to {args.save_npz}")


if __name__ == "__main__":
    main()
