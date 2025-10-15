# mkplot_a11.py
# Replicates:
#   load exact.dat
#   load a11stw01.tbl
#   t=a11stw01(:,1); wl=a11stw01(:,2);
#   plot(t,wl,'k'); hold; plot(exact(1:4:601,1),exact(1:4:601,2),'ko');

import numpy as np
import matplotlib.pyplot as plt

# --- Load data (assumes files are in the same folder as this script) ---
a11 = np.loadtxt("a11stw01.tbl")
exact = np.loadtxt("exact.dat")

# Columns: MATLAB is 1-indexed; Python is 0-indexed
t = a11[:, 0]      # a11stw01(:,1)
wl = a11[:, 1]     # a11stw01(:,2)

# Subsample exact like exact(1:4:601,:) in MATLAB:
# -> rows 0..600 step 4 in Python
exact_sub = exact[0:601:4, :]
x_exact = exact_sub[:, 0]
eta_exact = exact_sub[:, 1]

# --- Plot ---
plt.figure(figsize=(8, 4))  # ~2:1 aspect like PlotBoxAspectRatio [2 1 1]
plt.plot(t, wl, '-', linewidth=1.5)              # black line (default color ok)
plt.plot(x_exact, eta_exact, 'o', markersize=4)  # black circles

plt.xlabel("time [s]")
plt.ylabel("water level [m]")
plt.title(r"kh=0.55$\pi$; depth averaged")
plt.tight_layout()
plt.grid(True, alpha=0.3)

# Save like MATLAB's "print -dpng seiche1.png"
plt.savefig("seiche1.png", dpi=150)
plt.show()
