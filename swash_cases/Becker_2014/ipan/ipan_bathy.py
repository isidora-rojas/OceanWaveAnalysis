# Create SWASH .bot bathymetry files in /mnt/data for the user's 1D transect.
# - Reef flat: 0–400 m, depth = 0.5 m
# - Reef face: 400–580 m, slope = 0.06 (depth increases linearly to 11.3 m at 580 m)
# - Offshore buffer: 580–650 m, constant 11.3 m
# Also create tide-shifted variants for a mean 0.5 m and range 0.7 m (±0.35 m).

import numpy as np
from pathlib import Path

out_dir = Path("Becker_2014/")
out_dir.mkdir(parents=True, exist_ok=True)

dx = 1.0
x_end = 650.0   # domain 0..650 m
nx = int(x_end / dx) + 1  # include endpoint
x = np.linspace(0.0, x_end, nx)

# Base depth profile
depth = np.zeros_like(x)

# Reef flat: 0–400 m at 0.5 m depth
depth[(x >= 0) & (x <= 400)] = 0.5

# Reef face: 400–580 m, slope 0.06
mask_slope = (x > 400) & (x <= 580)
depth[mask_slope] = 0.5 + 0.06 * (x[mask_slope] - 400.0)

# Offshore buffer: 580–650 m at 11.3 m depth
depth[x > 580] = 11.3

# Write helper to save a .bot file (two columns: x depth)
def write_bot(filepath: Path, xvals, zvals):
    header = "! x [m], depth [m, positive downward]\n"
    # Build 2-col array
    arr = np.column_stack([xvals, zvals])
    # Save with header; SWASH ignores lines starting with '!'
    with open(filepath, "w") as f:
        f.write(header)
        np.savetxt(f, arr, fmt=["%.2f", "%.4f"])

# Base file
base_path = out_dir / "reef_transect.bot"
write_bot(base_path, x, depth)


# Show a quick preview of the first 10 lines of the base file
preview_lines = 12
with open(base_path, "r") as f:
    preview = "".join([next(f) for _ in range(preview_lines)])

preview, str(base_path), [str(p) for p in tide_paths[:4]]
A