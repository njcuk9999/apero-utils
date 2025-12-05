from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from typing import Union
import warnings
from scipy.optimize import curve_fit
# import univariate spline as uis
from scipy.interpolate import InterpolatedUnivariateSpline as uis
from scipy.constants import c


# =========================
# Gaussian + Slope Function
# =========================
def gauss_fit_s(x: Union[float, np.ndarray], a: float, x0: float, sigma: float,
                zp: float, slope: float) -> Union[float, np.ndarray]:
    """
    Gaussian fit with a slope

    :param x: numpy array (1D), the x values for the gauss fit
    :param a: float, the amplitude
    :param x0: float, the mean position
    :param sigma: float, the FWHM
    :param zp: float, the dc level
    :param slope: the float (x-x0) * slope

    :return: np.ndarray - the gaussian value with slope correction
    """
    # calculate gaussian
    with warnings.catch_warnings(record=True) as _:
        gauss = a * np.exp(-0.5 * (x - x0) ** 2 / (sigma ** 2)) + zp
    correction = (x - x0) * slope
    return gauss + correction


# =========================
# Peak Finder
# =========================
def is_peak(tmp_flux):
    """
    Function to find peaks in the flux data.
    A peak is defined as a point that is higher than its immediate neighbors.
    """
    tmp_flux = np.array(tmp_flux, dtype=float)
    # which point is higher than either of its neighbors (+-2 steps)

    roll_n2 = np.roll(tmp_flux, 2)  # flux shifted by -2
    roll_n1 = np.roll(tmp_flux, 1)  # flux shifted by -1
    roll_p1 = np.roll(tmp_flux, -1)  # flux shifted by +1
    roll_p2 = np.roll(tmp_flux, -2)  # flux shifted by +2

    # Find indices where the value is greater than all four neighbors
    peak = np.where((tmp_flux > roll_n2) * (tmp_flux > roll_n1) * (tmp_flux > roll_p1) * (tmp_flux > roll_p2))[0]

    return peak


# =========================
# Data Loading
# =========================
# Load wavelength solution (2D array: order x pixel)
wave = np.array(fits.getdata('2F3798BAE7a_pp_e2dsff_AB_wavesol_ref_AB.fits'), dtype=float)
# Load flux data (2D array: order x pixel)
flux = np.array(fits.getdata('F883D6B407c_pp_e2dsff_AB.fits'), dtype=float)
hc_box_size = 5  # half-width of the box around each peak for fitting

# =========================
# Output Arrays
# =========================
nsig_all = []  # S/N of each detected line
fwhm_all = []  # FWHM (km/s) of each detected line
wave_line_all = []  # Wavelength of each detected line
orders_all = []  # Order index for each detected line

index0 = np.arange(len(wave[0]))  # pixel indices for one order

# =========================
# Main Loop Over Orders
# =========================
for iord in range(len(wave)):
    print(iord)
    peaks = is_peak(flux[iord])  # Find peaks in this order

    # Only keep peaks that are not too close to the edges
    keep = np.where((peaks > hc_box_size) & (peaks < len(flux[iord]) - hc_box_size))[0]
    peaks = peaks[keep]

    # Spline to interpolate wavelength as a function of pixel index
    spl_wave = uis(index0, wave[iord], k=3)

    # Loop over all peaks in this order
    for ipeak in range(len(peaks)):
        # Extract a segment of the flux around the peak
        segment = flux[iord][peaks[ipeak] - hc_box_size:peaks[ipeak] + hc_box_size]
        if False in np.isfinite(segment):
            print(f'Order {iord + 1}, Peak {ipeak + 1}: Segment contains NaN or Inf values, skipping...')
            continue

        index = index0[peaks[ipeak] - hc_box_size:peaks[ipeak] + hc_box_size]

        ymax = np.max(segment)
        ymin = np.min(segment)
        posmax = index[np.argmax(segment)]  # pixel index of max

        guess_hc_ewid = 1  # initial guess for width
        guess = [ymax - ymin, posmax, guess_hc_ewid, ymin, 0]  # [ampl, center, width, offset, slope]

        try:
            fit, cov = curve_fit(gauss_fit_s, index, segment, p0=guess)
        except RuntimeError as e:
            # Fit failed, skip this peak
            continue
        recon = gauss_fit_s(index, *fit)  # reconstructed fit

        nsig = fit[0] / np.nanstd(segment - recon)  # S/N of the line

        if nsig < 3:
            # S/N too low, skip
            continue

        wave_line = spl_wave(fit[1])  # wavelength at fitted center

        # FWHM in km/s, using the local wavelength solution
        fwhm = (spl_wave(fit[1] + fit[2]) / wave_line - 1) * c * 2.3548 / 1000  # fwhm in km/s

        # Store results
        nsig_all.append(nsig)
        fwhm_all.append(fwhm)
        wave_line_all.append(wave_line)
        orders_all.append(iord)

        rms = np.nanstd(segment - recon)
        nsig = fit[0] / rms
        # print(f'Order {iord + 1}, Peak {ipeak + 1}: FWHM = {fit[2]:.3f}, NSIG = {nsig:.3f}')

    plt.plot(wave[iord], flux[iord], label=f'Order {iord + 1}')
    # plt.scatter(wave[iord][peaks], flux[iord][peaks], color='red')

plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux (e-/s)')

# =========================
# Convert Lists to Arrays
# =========================
nsig_all = np.array(nsig_all)
fwhm_all = np.array(fwhm_all)
wave_line_all = np.array(wave_line_all)
orders_all = np.array(orders_all)

# =========================
# Quality Cuts
# =========================
keep = (fwhm_all > 3) * (fwhm_all < 7) * (nsig_all > 20)
nsig_all = nsig_all[keep]
fwhm_all = fwhm_all[keep]
wave_line_all = wave_line_all[keep]
orders_all = orders_all[keep]

# =========================
# Merge Duplicates Across Orders
# =========================
keep = np.ones_like(wave_line_all, dtype=bool)
# loop such that lines overlapping in orders are merged together
for i in range(len(wave_line_all)):
    # find doublets where 2 lines are within 1km/s
    g = np.abs(wave_line_all / wave_line_all[i] - 1) * c < 1000
    if np.sum(g) > 1:
        print(f'Merging {np.sum(g)} lines for order {orders_all[i] + 1} at {wave_line_all[i]:.3f} nm')
        # merge lines
        wave_line_all[i] = np.mean(wave_line_all[g])
        nsig_all[i] = np.mean(nsig_all[g])
        fwhm_all[i] = np.mean(fwhm_all[g])
        orders_all[i] = np.min(orders_all[g])
        tmp = keep[g]
        tmp[1:] = False  # keep only the first occurrence
        keep[g] = tmp

# Keep only merged lines
nsig_all = nsig_all[keep]
fwhm_all = fwhm_all[keep]
wave_line_all = wave_line_all[keep]
orders_all = orders_all[keep]

# =========================
# Final Plot
# =========================
plt.plot(wave_line_all, np.zeros_like(wave_line_all), 'ko', markersize=2, label='Lines')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux (e-/s)')
plt.title('Flux vs Wavelength')
plt.grid()
plt.show()