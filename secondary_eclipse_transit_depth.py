import numpy as np
from scipy.constants import h, c, k

# Planck function (spectral radiance) in W·m⁻²·sr⁻¹·m⁻¹
def planck_lambda(wavelength, T):
    wl = wavelength  # meters
    return (2*h*c**2 / wl**5) / (np.exp(h*c / (wl*k*T)) - 1)

# Band-integrated surface brightness ratio
def band_integrated_ratio(T_p, T_s, wavelengths, response):
    """
    Compute band-averaged planet/star surface brightness ratio.

    Parameters:
    T_p : float, planet brightness temperature [K]
    T_s : float, stellar effective temperature [K]
    wavelengths : array [m]
        Wavelength grid matching the response curve
    response : array
        Normalized filter transmission curve (same length as wavelengths)
    """
    I_p = planck_lambda(wavelengths, T_p)
    I_s = planck_lambda(wavelengths, T_s)

    # Weighted integrals over the bandpass
    num = np.trapezoid(I_p * response, wavelengths)
    den = np.trapezoid(I_s * response, wavelengths)
    return num / den

# Eclipse depth calculator with bandpass option
def secondary_eclipse_depth(Rp, Rs, T_p, T_s, a, Ag=0.0,
                             wavelengths=None, response=None,
                             wavelength=4.5e-6):
    """
    Calculate secondary eclipse depth.

    Parameters:
    Rp : float
        Planet radius [m]
    Rs : float
        Stellar radius [m]
    T_p : float
        Planet dayside brightness temperature [K]
    T_s : float
        Stellar effective temperature [K]
    a : float
        Semi-major axis [m]
    Ag : float, optional
        Geometric albedo (default=0)
    wavelengths : array, optional
        Wavelength grid [m] for band integration
    response : array, optional
        Filter transmission curve (normalized)
    wavelength : float, optional
        Single central wavelength [m] (if no bandpass given)

    Returns:
    depth : float
        Eclipse depth (fraction)
    """

    if wavelengths is not None and response is not None:
        ratio = band_integrated_ratio(T_p, T_s, wavelengths, response)
    else:
        I_p = planck_lambda(wavelength, T_p)
        I_s = planck_lambda(wavelength, T_s)
        ratio = I_p / I_s

    thermal = (Rp / Rs)**2 * ratio
    reflected = Ag * (Rp / a)**2  # full phase (Phi = 1)
    depth = thermal + reflected
    return depth



def eclipse_snr(depth_ppm, duration_min, sigma30_ppm=400, N_tel=1, N_ecl=1):
    """
    Estimate SNR of a secondary eclipse with NGTS.

    Parameters
    ----------
    depth_ppm : float
        Eclipse depth in ppm (parts per million).
    duration_min : float
        Eclipse duration in minutes.
    sigma30_ppm : float, optional
        RMS noise in ppm per 30 minutes (default = 400 ppm, typical single NGTS telescope).
        Use ~152 ppm for combined NGTS telescopes, or ~1000 ppm for conservative spec.
    N_tel : int, optional
        Number of NGTS telescopes combined simultaneously (default = 1).
    N_ecl : int, optional
        Number of eclipses observed/stacked (default = 1).

    Returns
    -------
    snr : float
        Signal-to-noise ratio of the eclipse detection.
    sigma_eclipse : float
        Effective noise level for one eclipse (ppm).
    """
    # Noise per eclipse (binned to eclipse duration)
    sigma_eclipse = sigma30_ppm * np.sqrt(30.0 / duration_min)
    
    # SNR for a single eclipse, single telescope
    snr_single = depth_ppm / sigma_eclipse
    
    # Scale by multiple eclipses and telescopes
    snr_total = snr_single * np.sqrt(N_tel * N_ecl)
    
    return snr_total, sigma_eclipse


print(f"Noise per eclipse: {sigma_eclipse:.1f} ppm")
print(f"SNR (N_tel={N_tel}, N_ecl={N_ecl}): {snr:.2f}")


