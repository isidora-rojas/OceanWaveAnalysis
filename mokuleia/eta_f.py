"""
Computation of wave setup at the reef face sensor using equation 4 in Becker(2014) "Water level effects on 
breaking wave setup for Pacific Island fringing reefs".

equation 4: eta_f = - (H_f^2 * k_f) / (8 sinh(2 * k_f * h_f) )

Each variable is seperately computed as:

H_f: rms wave height at reef face
k_f: wavenumber at reef face
h_f: mean water level at reef face
"""
import numpy as np
from BulkWaveStats import sig_wave_height


def H_f(
    Seta: np.ndarray,
    freqs: np.ndarray,
    t_spec: np.ndarray,
    t1: np.ndarray,
    h1: np.ndarray,
    bands: tuple[float, float],
    *,
) -> tuple[np.ndarray, np.ndarray]:
    ''' Return RMS wave height and matching time centers. Feelin a bit lazy so going to compute off of sig wave height function sig_wave_height.
    Since Hs is 4 sqrt(sigma) so Hrms = sqrt(2)/2 * Hs

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
    band: list []
        frequency band interested in integrated under

    Returns
    --------
    )
    H_f: np.ndarray
        root mean square wave height 
    '''
    Hs_tot, Hs_ig, Hs_ss, Hs_input, Tp_ss, t_center = sig_wave_height(Seta, freqs, t_spec, t1, h1)
    Hrms = 0.5 * np.sqrt(2) * Hs_input

    return Hrms, t_center



