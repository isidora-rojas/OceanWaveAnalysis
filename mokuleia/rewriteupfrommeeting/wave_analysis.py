import numpy as np
from scipy import signal
from numpy import pi

# --- inputs you must set ---
pp = pp_dtrnd                  # detrended pressure time series (Pa) or (dbar -> convert)
dt = 1.0                       # sampling [s]
rho = 1025.0                   # seawater density [kg/m^3]
g = 9.81
h = 12.0                       # mean water depth [m] (replace with your site)
z = -1.5                       # sensor elevation relative to SWL [m] (negative downward)
Fs = 1.0/dt
fNy = Fs/2

# bands (tune as needed)
f_ig = (0.004, 0.040)
f_ss = (0.050, 0.500)

# Welch settings → one point per window
nperseg = 4096                 # ~68 min @ 1 Hz. Use 2048 (~34 min) if you want finer time sampling
noverlap = nperseg // 2
window = 'hann'
nfft = nperseg                 # can use next pow2

# --- helper: dispersion solver k(ω) (vectorized Newton) ---
def wavenumber(omega, h, g=9.81, max_iter=50, tol=1e-12):
    # initial guess: deep water
    k = (omega**2) / g
    # avoid zero/NaN for ω=0
    k = np.where(omega > 0, k, 0.0)
    for _ in range(max_iter):
        f = g*k*np.tanh(k*h) - omega**2
        df = g*np.tanh(k*h) + g*k*h*(1/np.cosh(k*h))**2
        step = f/np.where(df!=0, df, np.inf)
        k_new = k - step
        if np.nanmax(np.abs(step)) < tol:
            return k_new
        k = np.where(np.isfinite(k_new), k_new, k)
    return k

# --- compute Welch PSD of *pressure* ---
freqs, Spp = signal.welch(
    pp, fs=Fs, window=window, nperseg=nperseg, noverlap=noverlap, nfft=nfft, detrend='constant', scaling='density'
)
# freqs: [0..fNy], one-sided; Spp units: pressure^2 / Hz

# --- convert Spp -> surface-elevation spectrum Sηη using transfer function ---
omega = 2*pi*freqs
k = wavenumber(omega, h, g)

# pressure response function (sensor at z below SWL)
# Avoid division by zero at f=0:
with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
    G = (1.0/(rho*g)) * (np.cosh(k*h)/np.cosh(k*(h+z)))    # amplitude gain p->η
    # optional: cap extreme gains at high f to reduce noise amplification
    G = np.where(np.isfinite(G), G, 0.0)
    G = np.clip(G, -50, 50)  # gentle cap; adjust to your site/instrument

See = (G**2) * Spp   # surface elevation PSD [m^2/Hz]

# --- integrate moments over bands ---
def band_mask(freqs, fmin, fmax):
    return (freqs >= fmin) & (freqs <= min(fmax, fNy))

ms0_ig = np.trapz(See[band_mask(freqs, *f_ig)], freqs[band_mask(freqs, *f_ig)])
ms0_ss = np.trapz(See[band_mask(freqs, *f_ss)], freqs[band_mask(freqs, *f_ss)])
ms0_all = np.trapz(See[(freqs>0)], freqs[(freqs>0)])

Hs_ig  = 4*np.sqrt(ms0_ig)
Hs_ss  = 4*np.sqrt(ms0_ss)
Hs_all = 4*np.sqrt(ms0_all)

# --- peak period (sea/swell band) ---
m = band_mask(freqs, *f_ss)
fp_ss = freqs[m][np.argmax(See[m])] if np.any(m) else np.nan
Tp_ss = 1.0/fp_ss if fp_ss and fp_ss>0 else np.nan
