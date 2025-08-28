"""
===============================================================================
Code Steps (matches main code steps):
===============================================================================
    1. Change to working directory.
    2. Load and filter the Uranium-Neon line catalog.
    3. Deduplicate lines, keeping only the brightest within 20 km/s.
    4. Define function to get approximate wavelength for an order.
    5. Load hollow cathode (HC) and Fabry-Perot (FP) spectra.
    6. Prepare order list (start from middle order and alternate outwards).
    7. Initialize cavity fit and bookkeeping arrays.
    8. If enough previous solutions exist, robustly fit the cavity using all orders.
    9. Check for integer offset in cavity numbering and correct if needed.
    10. Fill missing values and robustly filter orders for cavity fit.
    11. Final robust fit to the cavity using all FP data and user validation.
    12. Main loop over orders to build wavelength solution:
        a. Extract 1D spectrum for this order.
        b. Find the NLINES brightest HC lines in this order.
        c. Find Fabry-Perot peaks for this order.
        d. Estimate the step between FP peaks and fit a polynomial to it.
        e. Normalize the step size to integer multiples of the fit.
        f. Assign a running index to each FP peak (cavity order).
        g. Remove outliers in the FP peak sequence.
        h. Create a synthetic spectrum with spikes at the HC line positions.
        i. Get the approximate wavelength range for this order.
        j. Estimate the range of possible FP cavity orders for this order.
        k. Try all possible FP cavity order guesses and find the best alignment.
        l. Normalize the nvalid2 array for plotting.
        m. Plot the results and ask the user for validation.
        n. If user accepts, save the wavelength solution and pickle.
    13. Build the final 2D wavelength solution for all orders:
        a. For each order, fit a Chebyshev polynomial to the wavelength solution.
        b. Smooth the Chebyshev coefficients across orders to remove outliers.
        c. For each order, compute the final wavelength solution and store coefficients in header.
    14. Save the final wavelength solution to a FITS file.
    15. Save the cavity fit coefficients to a file.
    16. Plot the final wavelength solution for all orders.
===============================================================================
"""
from scipy.optimize import curve_fit
import scipy.optimize
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from astropy.table import Table
import os
import scipy.optimize
from tqdm import tqdm
from typing import List, Optional, Union, Any, Tuple
import scipy
from apero.core.math import lowpassfilter
import pickle
import glob
from scipy.interpolate import InterpolatedUnivariateSpline as ius

# TODO: Loop 3 times
#       - first time will be bad
#       - second time should be good
#       - third time to make sure
# --- CONSTANTS AND PARAMETERS -------------------------------------------------
CAVITY0 = 2.4004e7  # Initial guess for cavity length (in nm)   # TODO: ASK USER
WAVE_DOMAIN = [955, 2450]  # First to last orders (nm)   # TODO: ASK USER
N_ORDERS = 49  # Number of spectral orders
WAVE_APPROX = 0.05  # Fractional range for approximate wavelength
NLINES = 50  # Number of lines to use per order
WAVEDEGN = 5  # Degree of polynomial for wavelength solution
FP_PEAK_STEP_POLY_DEG = 1  # Degree for robust_polyfit of FP peak step
FP_PEAK_POLY_DEG = 3  # Degree for robust_polyfit of FP peak count
FP_PEAK_KEEP_PERCENTILE = 90  # Percentile for FP peak amplitude threshold
FP_PEAK_KEEP_FACTOR = 0.3  # Factor for FP peak amplitude threshold
FP_PEAK_STEP_WINDOW = 5  # Step for polyfit in wave_guess/fp_pix
FP_POLY_DEG = 5  # Degree for polyfit in wave_guess/fp_pix
NSIG_ACCEPT_FP = 15   # Number of sigma to blindly accept an order as valid
FP_STEP_OFFS_RANGE = 20  # Range for integer offset search
FP_STEP_OFFS_STEP = 1  # Step for integer offset search
FP_STEP_MAD_BINS = 40  # Number of bins for histogram in mini
FP_STEP_MAD_RANGE = (-10, 10)  # Range for histogram in mini
FP_STEP_MAD_THRESHOLD = 2.0  # Threshold for mini selection
FP_STEP_VALID_MIN = 5  # Minimum number of valid points for fit
CHEBY_FIT_DEG = 5  # Degree for Chebyshev fit
ROBUST_POLYFIT_DEG = 7  # Degree for robust_polyfit on Chebyshev coeffs
ROBUST_POLYFIT_SIGMA = 3  # Sigma cut for robust_polyfit on Chebyshev coeffs
hc_window_kms = 50  # Window size for HC lines (in km/s)

path = '/data/spip/misc/wavesol/2025-05-09/'

# hc_spectrum_file = '/Users/eartigau/spip/test2/180881c_pp_e2ds_AB.fits'
# fp_spectrum_file = '/Users/eartigau/spip/test2/180879a_pp_e2ds_AB.fits'
hc_spectrum_file = '/scratch2/spip/misc/download/2025-05-23/1F5B7E10A5c_pp_e2dsff_AB.fits'
fp_spectrum_file = '/scratch2/spip/misc/download/2025-05-23/36ECD2422Aa_pp_e2dsff_AB.fits'

working_dir = '/scratch2/spip/misc/statics/2025-05-23'


# --- Utility Functions --------------------------------------------------------
def save_pickle(data: object, filename: str) -> None:
    """Save data to a pickle file."""
    with open(filename, 'wb') as f:
        pickle.dump(data, f)


def load_pickle(filename: str) -> object:
    """Load data from a pickle file."""
    with open(filename, 'rb') as f:
        return pickle.load(f)


def mad(v: np.ndarray) -> float:
    """Compute the median absolute deviation (MAD) of an array."""
    return np.nanmedian(np.abs(v - np.nanmedian(v)))


def gauss(x: float, mu: float, amp: float, sigma: float, zp: float) -> float:
    """Gaussian function."""
    return amp * np.exp(-0.5 * ((x - mu) / sigma) ** 2) + zp


def get_lines_pix(sp_tmp, fp=False):
    """
    Detect spectral line peaks and estimate their positions, FWHM, and amplitudes.
    """
    sp = np.array(sp_tmp.copy(), dtype=float)

    # --- Remove low-frequency background for HC spectra ---
    if not fp:
        sp -= lowpassfilter(sp, 15)

    # --- Clean up spectrum: set non-finite and edge values to zero ---
    sp[~np.isfinite(sp)] = 0
    sp[0:5] = 0
    sp[-5:] = 0

    # --- Find local maxima (peaks) in the spectrum ---
    peaks = np.where((sp > np.roll(sp, 1)) & (sp > np.roll(sp, -1)))[0]

    peaks_pixels = []
    peaks_max = []
    fwhm_pixels = []

    for ipeak in range(len(peaks)):
        if not fp:
            pixidbit = [-1, 0, 1]
            pixbit = sp[peaks[ipeak] - 1:peaks[ipeak] + 2] / sp[peaks[ipeak]]
            if np.min(pixbit) < 0:
                continue
            fit = np.polyfit(pixidbit, pixbit, 2)
            fwhm = 2 * np.sqrt(-0.5 / fit[0])
            if (fwhm < 1) or (fwhm > 5):
                continue
            peak_pos = -0.5 * fit[1] / fit[0]
            peaks_pixels.append(peak_pos + peaks[ipeak])
            fwhm_pixels.append(fwhm)
            peaks_max.append(sp[peaks[ipeak]])
        else:
            if ipeak == 0 or ipeak == len(peaks) - 1:
                continue
            w = (peaks[ipeak + 1] - peaks[ipeak - 1]) / 4.0
            w = int(np.round(w))
            pixidbit = np.arange(-w, w + 1)
            pixbit = sp[peaks[ipeak] - w:peaks[ipeak] + w + 1]
            try:
                fit, _ = scipy.optimize.curve_fit(gauss, pixidbit, pixbit, p0=[0, sp[peaks[ipeak]], 1, 0])
                fwhm = 2.35482 * fit[2]
                peaks_pixels.append(fit[0] + peaks[ipeak])
                fwhm_pixels.append(fwhm)
                peaks_max.append(fit[1])
            except:
                continue

    mus = np.array(peaks_pixels)
    sigmas = np.array(fwhm_pixels)
    flux = np.array(peaks_max)

    keep = np.isfinite(mus) & np.isfinite(sigmas) & np.isfinite(flux) & (flux > 0)
    mus = mus[keep]
    sigmas = sigmas[keep]
    flux = flux[keep]

    ord = np.argsort(mus)
    mus = mus[ord]
    sigmas = sigmas[ord]
    flux = flux[ord]

    if fp:
        # --- 12d. Estimate the step between FP peaks and fit a polynomial to it ---
        peak_step = np.diff(mus)
        mus = mus[1:]
        flux = flux[1:]

        for ite in range(2):
            fit, keep = robust_polyfit(mus, peak_step, FP_PEAK_STEP_POLY_DEG, 3)
            peak_step = peak_step[keep]
            mus = mus[keep]
            flux = flux[keep]

        # --- 12e. Normalize the step size to integer multiples of the fit ---
        nn = np.round(peak_step / np.polyval(fit, mus))
        peak_step /= nn

        # --- 12f. Assign a running index to each FP peak (cavity order) ---
        peak_count = np.zeros(len(mus), dtype=int)
        for i in range(1, len(mus)):
            step_from_last = (mus[i] - mus[i - 1]) / np.polyval(fit, mus[i])
            peak_count[i] = peak_count[i - 1] + np.round(step_from_last).astype(int)

        # --- 12g. Remove outliers in the FP peak sequence ---
        fit, keep = robust_polyfit(peak_count, mus, FP_PEAK_POLY_DEG, 7)

        residual = mus - np.polyval(fit, peak_count)
        # plt.plot(peak_count, residual, 'k.')
        # plt.plot(peak_count[keep], residual[keep], 'r.')
        # plt.xlabel('Peak count')
        # plt.ylabel('Residual')
        # plt.title('Residuals of FP peaks')
        # plt.show()

        mus = mus[keep]
        # mu_pix = mu_pix[keep]
        flux = flux[keep]
        peak_count = peak_count[keep]

        return mus, sigmas, flux, peak_count
    else:
        return mus, sigmas, flux


def fit_cheby(xvector: np.ndarray, yvector: np.ndarray, deg: int,
              domain: List[float], weight: Optional[np.ndarray] = None
              ) -> Union[np.ndarray, Any]:
    """
    Fit a chebyshev polynomial in form y(x) = T0(x) + T1(x) + ... Tn(x)
    returns the chebyshev polynomial coefficients
    """
    domain_cheby = 2 * (xvector - domain[0]) / (domain[1] - domain[0]) - 1
    coeffs = np.polynomial.chebyshev.chebfit(domain_cheby, yvector, deg, w=weight)
    return coeffs


def val_cheby(coeffs: np.ndarray, xvector: Union[np.ndarray, int, float],
              domain: List[float]) -> Union[np.ndarray, int, float]:
    """
    Using the output of fit_cheby calculate the fit to x  (i.e. y(x))
    where y(x) = T0(x) + T1(x) + ... Tn(x)
    """
    domain_cheby = 2 * (xvector - domain[0]) / (domain[1] - domain[0]) - 1
    yvector = np.polynomial.chebyshev.chebval(domain_cheby, coeffs)
    return yvector


def get_peak0_guess(iord):
    """
    Get the initial guess for the peak0 value.
    This function is a placeholder and should be implemented based on specific requirements.
    """

    pkls = glob.glob('wave_order_*.pkl')

    if len(pkls) < 5:
        return np.nan

    orders = []
    peak0_guesses = []
    for pkl in pkls:
        dict_fp = load_pickle(pkl)
        orders.append(int(pkl.split('_')[-1].split('.')[0]))
        peak0_guesses.append(dict_fp['peak0_guess'])

    peak0_guesses = np.array(peak0_guesses)
    peak0_guesses = peak0_guesses[~np.isnan(peak0_guesses)]

    ord = np.argsort(orders)
    peak0_guesses = peak0_guesses[ord]
    orders = np.array(orders)[ord]

    fit = robust_polyfit(orders, peak0_guesses, 1, 3)[0]

    return np.polyval(fit, iord)


def robust_polyfit(x: np.ndarray, y: np.ndarray, degree: int, nsigcut: float, accept_width: Optional[float] = None) -> \
Tuple[np.ndarray, np.ndarray]:
    """
    Perform a robust polynomial fit to the data, iteratively rejecting outliers.
    """
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    degree = np.array(degree, dtype=int)
    if accept_width is not None:
        nsigcut = 1.0
    keep = np.isfinite(y)
    nsigmax = np.inf
    fit = None
    while nsigmax > nsigcut:
        fit = np.polyfit(x[keep], y[keep], degree)
        res = y - np.polyval(fit, x)
        if accept_width is None:
            sig = np.nanmedian(np.abs(res))
        else:
            nsigcut = 1
            sig = accept_width
        if sig == 0:
            nsig = np.zeros_like(res)
            nsig[res != 0] = np.inf
        else:
            nsig = np.abs(res) / sig
        nsigmax = np.max(nsig[keep])
        keep = nsig < nsigcut
    return fit, keep


class YesNoButtonGraph:
    def __init__(self, fig):
        self.fig = fig
        self.response = 'n'
        self.question = ''
        self.yes_button = None
        self.no_button = None

    def add(self, question):
        self.question = str(question)
        # Add question text centered above the buttons
        self.fig.text(0.5, 0.20, question, ha='center', va='center',
                      fontsize=12)
        # Move subplots up to leave room for buttons/text
        plt.subplots_adjust(bottom=0.4)
        # Define button positions (x, y, width, height)
        yes_ax = self.fig.add_axes([0.42, 0.08, 0.08, 0.05])
        no_ax = self.fig.add_axes([0.52, 0.08, 0.08, 0.05])
        # Create buttons
        self.yes_button = Button(yes_ax, 'Yes', color='lightgreen',  hovercolor='green')
        self.no_button = Button(no_ax, 'No', color='lightcoral', hovercolor='red')
        self.yes_button.on_clicked(self.on_yes)
        self.no_button.on_clicked(self.on_no)

    def on_yes(self, event):
        print(f'{self.question}: yes')
        self.response = 'y'
        plt.close(self.fig)

    def on_no(self, event):
        print(f'{self.question}: no')
        self.response = 'n'
        plt.close(self.fig)


# --- MAIN PIPELINE ------------------------------------------------------------

if __name__ == "__main__":

    # --- 1. Change to working directory ---
    os.chdir(working_dir)

    # --- 2. Load and filter the Uranium-Neon line catalog ---
    #   Our defaut comes from Redman et al 2011: https://iopscience.iop.org/article/10.1088/0067-0049/195/2/24
    #   via Vizier: https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/ApJS/195/24/table2
    #   Note if downloading from Vizier - select FITS (binary) Table and set max rows to unlimited!

    tbl_hc_lines_file = 'J_ApJS_195_24_table2.dat.fits'
    tbl = Table.read(tbl_hc_lines_file)
    ion = tbl['Ion']
    keep = np.array([True if 'U' in ion[i] else False for i in range(len(ion))])
    tbl = tbl[keep]
    keep = (tbl['lambda'].data > WAVE_DOMAIN[0] * (1 - WAVE_APPROX)) * (
                tbl['lambda'].data < (WAVE_DOMAIN[1] * (1 + WAVE_APPROX)))
    tbl = tbl[keep]
    wave_ref0 = tbl['lambda'].data
    flux_ref0 = tbl['RFlux'].data

    # --- 3. Deduplicate lines, keeping only the brightest within 20 km/s ---
    keep = np.zeros_like(wave_ref0, dtype=bool)
    for i in range(len(wave_ref0)):
        dv = (1 - wave_ref0[i] / wave_ref0) * 3e5
        g = np.abs(dv) < hc_window_kms
        if flux_ref0[i] == np.max(flux_ref0[g]):
            keep[i] = True
    tbl = tbl[keep]
    wave_ref0 = tbl['lambda'].data

    # --- 4. Define function to get approximate wavelength for an order ---
    def get_approx_wave(ord):
        wave0 = 1 / np.polyval(np.polyfit([0., N_ORDERS], [1 / WAVE_DOMAIN[0], 1 / WAVE_DOMAIN[1]], 1), ord)
        return wave0, wave0 * (1 - WAVE_APPROX), wave0 * (1 + WAVE_APPROX)


    # --- 5. Load hollow cathode (HC) and Fabry-Perot (FP) spectra ---
    sp1 = fits.getdata(hc_spectrum_file)
    fp1 = fits.getdata(fp_spectrum_file)
    norders = sp1.shape[0]

    # --- 6. Prepare order list (start from middle order and alternate outwards) ---
    orders = np.arange(norders)
    orders = orders[np.argsort(np.abs(orders - norders / 2))]

    # --- 7. Initialize cavity fit and bookkeeping arrays ---
    flag_known_cavity = False
    n_pickles = len(glob.glob('wave_order_*.pkl'))
    all_fp_wave = np.zeros(0)
    all_int_fp = np.zeros(0)
    orders_known = []
    orders_mid = []

    fit_cavity = [CAVITY0]

    # --- 8. If enough previous solutions exist, robustly fit the cavity using all orders ---
    if n_pickles > 5:
        flag_known_cavity = True
        iord_wave_center = np.zeros(len(orders)) + np.nan
        iord_cavity_center = np.zeros(len(orders)) + np.nan
        fit_cavity = None

        for ii, iord in enumerate(orders):
            wave_order_file = f'wave_order_{iord}.csv'
            if os.path.exists(wave_order_file):
                tbl = Table.read(wave_order_file)
                orders_known.append(iord)
                orders_mid.append(np.mean(tbl['wavelength']))
                dict_fp = load_pickle(wave_order_file.replace('.csv', '.pkl'))

                if fit_cavity is None:
                    # --- Initial fit for the first order ---
                    # fit_cavity, _ = robust_polyfit(
                    #    dict_fp['fp_wave'],
                    #    dict_fp['fp_wave'] * dict_fp['int_fp'],
                    #    1, 3
                    # )
                    fit_cavity = [np.nanmedian(dict_fp['fp_wave'] * dict_fp['int_fp'])]
                    all_fp_wave = np.array(dict_fp['fp_wave'])
                    all_int_fp = np.round(np.polyval(fit_cavity, all_fp_wave) / all_fp_wave).astype(int)
                    all_ord = np.zeros(len(dict_fp['fp_wave'])) + iord
                else:

                    # --- Update fit for subsequent orders ---
                    cavity = np.polyval(fit_cavity, dict_fp['fp_wave'])
                    prev_int_fp = dict_fp['int_fp']
                    new_int_fp = np.round(cavity / dict_fp['fp_wave']).astype(int)
                    if False in (new_int_fp == prev_int_fp):
                        print(f"Warning: Intensity values have changed for order {iord}.")
                        dict_fp['int_fp'] = new_int_fp
                    all_fp_wave = np.append(all_fp_wave, dict_fp['fp_wave'])
                    all_int_fp = np.append(all_int_fp, dict_fp['int_fp'])
                    all_ord = np.append(all_ord, np.zeros(len(dict_fp['fp_wave'])) + iord)
                    save_pickle(dict_fp, wave_order_file.replace('.csv', '.pkl'))
                    fit_cavity, _ = robust_polyfit(
                        all_fp_wave,
                        all_fp_wave * all_int_fp,
                        3, 8
                    )
                    iord_wave_center[ii] = np.median(dict_fp['fp_wave'])
                    iord_cavity_center[ii] = np.median(dict_fp['fp_wave'] * dict_fp['int_fp'])

                # TODO: Plot afterwards?
                plt.plot(dict_fp['fp_wave'], dict_fp['fp_wave'] * dict_fp['int_fp'], '.', alpha=0.3)
                med_wave = np.nanmedian(dict_fp['fp_wave'])
                med_cavity = np.nanmedian(dict_fp['fp_wave'] * dict_fp['int_fp'])

                plt.text(med_wave, med_cavity, '\n\n\n' + str(iord), fontsize=8)

        # --- 10. Fill missing values and robustly filter orders for cavity fit ---
        orders_known = np.array(orders_known)
        orders_mid = np.array(orders_mid)
        all_ord = np.array(all_ord)
        invalid = ~np.isfinite(iord_wave_center)
        iord_wave_center[invalid] = np.polyval(np.polyfit(orders_known, orders_mid, 2), orders[invalid])
        iord_cavity_center[invalid] = np.polyval(fit_cavity, iord_wave_center[invalid])
        fit, keep_orders = robust_polyfit(orders_known, 1 / orders_mid, 1, 3)
        keep_iord = np.zeros_like(iord_wave_center, dtype=bool)
        for ii in range(len(iord_wave_center)):
            keep_iord[ii] = all_ord[ii] in orders_known[keep_orders]
        iord_wave_center = iord_wave_center[keep_iord]
        iord_cavity_center = iord_cavity_center[keep_iord]

        # --- 11. Final robust fit to the cavity using all FP data and user validation ---
        fit_cavity_tmp, _ = robust_polyfit(
            all_fp_wave,
            all_int_fp * all_fp_wave,
            7, 3
        )
        oo = np.argsort(all_fp_wave)

        fig, frame = plt.subplots(ncols=1, nrows=1)
        frame.plot(all_fp_wave[oo], np.polyval(fit_cavity_tmp, all_fp_wave[oo]), 'k-', label='Cavity fit')
        plt.legend()

        # yes / no button on graph
        yninst = YesNoButtonGraph(fig)
        yninst.add('Is this cavity valid and to be used as input?')
        plt.show(block=True)

        input_user = yninst.response

        if input_user == 'y':
            fit_cavity = fit_cavity_tmp

    for ite in range(2):
        # --- 12. Main loop over orders to build wavelength solution ---
        for iord in orders:
            wave_order_file = f'wave_order_{iord}.csv'
            if not os.path.exists(wave_order_file):
                print(f'Processing order {iord}...')

                # --- 12a. Extract 1D spectrum for this order ---
                sp = fits.getdata(hc_spectrum_file)[iord]

                # --- 12b. Find the NLINES brightest HC lines in this order ---
                linepix, _, flux = get_lines_pix(sp)
                if len(linepix) > NLINES:
                    oo = np.argsort(-flux)
                    linepix = linepix[oo[0:NLINES]]
                    flux = flux[oo[0:NLINES]]

                hdr_wavesol = fits.getheader(hc_spectrum_file)

                # --- 12c. Find Fabry-Perot peaks for this order ---
                fp_pix, mu_pix, amp, peak_count = get_lines_pix(fp1[iord], fp=True)

                # --- 12h. Create a synthetic spectrum with spikes at the HC line positions ---
                sp_test = np.zeros_like(sp)
                sp_test[linepix.astype(int)] = 1
                gg = gauss(np.arange(-5, 6), 0, 1, 1.5, 0)
                sp_test += np.convolve(sp_test, gg, mode='same')

                # --- 12i. Get the approximate wavelength range for this order ---
                wave0, wave_start, wave_end = get_approx_wave(iord)
                g = (wave_ref0 > wave_start) & (wave_ref0 < wave_end)
                wave_ref = wave_ref0[g]

                # --- 12j. Estimate the range of possible FP cavity orders for this order ---

                peak0_guess = get_peak0_guess(iord)

                if np.isnan(peak0_guess):
                    i_start = int(CAVITY0 / wave_end)
                    i_end = int(CAVITY0 / wave_start)
                    step = 0.2
                else:
                    i_start = peak0_guess - 200
                    i_end = peak0_guess + 200
                    step = 0.05

                peak0_guesses = np.arange(i_start, i_end, step)
                nvalid = np.zeros(len(peak0_guesses)) + np.nan
                nvalid2 = np.zeros(len(peak0_guesses)) + np.nan
                best_nvalid = 0

                # --- 12k. Try all possible FP cavity order guesses and find the best alignment ---
                for ii in tqdm(range(len(peak0_guesses))):

                    # --- Guess the peak0 value for this iteration ---
                    peak0_guess = peak0_guesses[ii]

                    # --- Compute the guessed wavelength solution for this peak0 ---
                    wave_guess = np.polyval(fit_cavity, wave0) / (peak0_guess - peak_count)

                    # --- Fit a polynomial between guessed wavelengths and FP pixel positions ---
                    fit = np.polyfit(
                        wave_guess[::FP_PEAK_STEP_WINDOW],
                        fp_pix[::FP_PEAK_STEP_WINDOW],
                        FP_POLY_DEG
                    )

                    # --- Map reference wavelengths to pixel positions using the fit ---
                    pix_ref = np.polyval(fit, wave_ref)
                    g = (pix_ref > 0) & (pix_ref < len(sp))

                    # --- If enough valid points, count how many HC lines match FP peaks ---
                    if np.sum(g) > FP_STEP_VALID_MIN:
                        pix_ref2 = pix_ref[g].astype(int)
                        nvalid[ii] = np.sum(sp_test[pix_ref2])

                    # --- Compute a normalized metric for plotting and selection ---
                    nvalid2[ii] = nvalid[ii] - np.nanmedian(nvalid[ii - 11:ii])

                    # --- If this guess is promising, refine the alignment ---
                    if nvalid2[ii] > 0:  # Only keep best

                        # --- Compute pixel offsets between HC and FP lines ---
                        mini = np.zeros(len(linepix))
                        mini_wave = np.zeros(len(linepix))

                        for i in range(len(linepix)):
                            imin = np.argmin(np.abs(pix_ref2 - linepix[i]))
                            mini[i] = pix_ref2[imin] - linepix[i]
                            mini_wave[i] = wave_ref[g][imin]

                        # --- Histogram the offsets to find the best alignment cluster ---
                        n, vals = np.histogram(
                            mini,
                            bins=FP_STEP_MAD_BINS,
                            range=FP_STEP_MAD_RANGE
                        )

                        g = np.abs(mini - vals[np.argmax(n)]) < FP_STEP_MAD_THRESHOLD

                        # --- If enough lines are well-aligned, fit a polynomial to them ---
                        if np.sum(g) > FP_STEP_VALID_MIN:
                            hc_pix2 = linepix[g]
                            mini2 = mini[g]
                            mini_wave2 = mini_wave[g]

                            if nvalid2[ii] > best_nvalid:
                                best_nvalid = nvalid2[ii]  # Update best_nvalid
                                # --- Robust fit to the best-matched lines ---
                                _, keep = robust_polyfit(hc_pix2, mini_wave2, 2, 3)
                                hc_pix2 = hc_pix2[keep]
                                mini_wave2 = mini_wave2[keep]

                                # --- Final polynomial fit for wavelength solution ---
                                # using warnings to remove pri
                                best_fit = np.polyfit(hc_pix2, mini_wave2, 3)
                                fp_wave = np.polyval(best_fit, fp_pix)

                                # --- Compute the cavity length from the polynomial fit ---
                                cavity = np.polyval(fit_cavity, fp_wave)

                                int_fp = np.array(np.round(cavity / fp_wave), dtype=int)

                                # --- Search for the best integer offset for FP orders ---
                                offs = np.arange(-FP_STEP_OFFS_RANGE, FP_STEP_OFFS_RANGE + 1, FP_STEP_OFFS_STEP)
                                mads = np.zeros(len(offs))

                                for ioff, off in enumerate(offs):
                                    mads[ioff] = mad((int_fp + off) * fp_wave)

                                # --- Compute the best-fit wavelength solution for all pixels ---
                                best_wave = np.polyval(best_fit, np.arange(len(sp)))
                                int_fp = int_fp + offs[np.argmin(mads)]

                                # --- Select valid FP peaks within the fitted wavelength range ---
                                valid = (fp_wave > np.min(mini_wave)) * (fp_wave < np.max(mini_wave))

                                # --- Only proceed if there are valid FP peaks ---
                                if len(fp_wave) > 0 and len(fp_pix) > 0:
                                    pass

                                # --- Store the FP solution for this guess ---
                                dict_fp = {
                                    'int_fp': int_fp[valid],
                                    'fp_wave': fp_wave[valid],
                                    'fp_pix': fp_pix[valid].astype(float),
                                    'peak0_guess': peak0_guess
                                }

                # --- 12l. Normalize the nvalid2 array for plotting ---
                nvalid2 /= mad(nvalid2)
                print(np.nanargmax(nvalid2), np.nanmax(nvalid2))

                if (ite == 0) * (np.nanmax(nvalid2) < NSIG_ACCEPT_FP):
                    print('No valid solution found for this order')
                    continue

                # --- 12m. Plot the results and ask the user for validation ---
                if (ite == 0) * (np.nanmax(nvalid2) > NSIG_ACCEPT_FP):
                    input_user = 'y'
                else:
                    continue

            # --- Create figure for diagnostics ---
            fig, ax = plt.subplots(2, 1, figsize=(10, 5))
            
            # --- Plot number of valid lines as a function of peak0 guess ---
            ax[0].plot(peak0_guesses, nvalid2, 'g-')
            ax[0].plot(peak0_guesses, nvalid2, 'r.')
            ax[0].set_xlabel('Peak0 guess')
            ax[0].set_ylabel('Number of valid lines')
            ax[0].set_title('Number of valid lines as a function of peak0 guess')
            
            # --- Plot spectrum and mark reference wavelengths ---
            ax[1].plot(best_wave, sp)
            keep = (wave_ref > np.min(best_wave)) & (wave_ref < np.max(best_wave))
            wave_ref2 = wave_ref[keep]
            for i in range(len(wave_ref2)):
                ax[1].axvline(wave_ref2[i], color='0.5', alpha=0.5)
            ax[1].set_yscale('log')
            floor_val = np.nanmedian(sp[sp != 0]) * 0.1
            ax[1].set_ylim(floor_val, np.nanmax(sp[sp != 0]))
            ax[1].set_xlabel('Wavelength')
            ax[1].set_ylabel('Flux')
            
            # yes / no button on graph
            yninst = YesNoButtonGraph(fig)
            yninst.add('Is this valid?')
            plt.show(block=True)

            input_user = yninst.response

            # --- 12n. If user accepts, save the wavelength solution and pickle ---
            if 'y' in input_user.lower():
                tbl = Table((np.arange(len(best_wave)), best_wave),
                            names=('pixel', 'wavelength'))
                tbl.write(wave_order_file, format='csv', overwrite=True)
                print(f"Saved to {wave_order_file}")
                pickle_file = wave_order_file.replace('.csv', '.pkl')
                save_pickle(dict_fp, pickle_file)

    # --- 13. Build the final 2D wavelength solution for all orders ---
    final_wave_sol = np.zeros((N_ORDERS, sp1.shape[1])) + np.nan
    hdr_wavesol = fits.Header()
    all_cheby = np.zeros((N_ORDERS, WAVEDEGN + 1)) + np.nan
    ith_coeff = 0

    middle_pixel_fp_peak = np.zeros(sp1.shape[0]) + np.nan
    cavity_middle = np.zeros(sp1.shape[0]) + np.nan

    # --- 13a. For each order, fit a Chebyshev polynomial to the wavelength solution ---
    for iord in range(sp1.shape[0]):
        pkl_file = f'wave_order_{iord}.pkl'
        if not os.path.exists(pkl_file):
            continue
        dict_fp = load_pickle(pkl_file)
        all_peaks = np.array(dict_fp['fp_wave'])
        cavity_residual = all_peaks - np.polyval(fit_cavity, all_peaks) / np.round(
            np.polyval(fit_cavity, all_peaks) / all_peaks)
        wave_from_cavity = np.polyval(fit_cavity, dict_fp['fp_wave']) / dict_fp['int_fp']

        fit_wave = np.polyfit(dict_fp['fp_pix'], wave_from_cavity, WAVEDEGN)
        final_wave = np.polyval(fit_wave, np.arange(sp1.shape[1]))
        fit_cheby_wave = fit_cheby(np.arange(sp1.shape[1]), final_wave, CHEBY_FIT_DEG, [0, sp1.shape[1] - 1])
        all_cheby[iord] = fit_cheby_wave

        wave_middle = final_wave[len(final_wave) // 2]
        cavity_middle[iord] = np.polyval(fit_cavity, wave_middle)
        middle_pixel_fp_peak[iord] = cavity_middle[iord] / wave_middle

    fit_middle_fp, keep = robust_polyfit(np.arange(sp1.shape[0]), middle_pixel_fp_peak, 5, 8)

    diff = middle_pixel_fp_peak - np.polyval(fit_middle_fp, np.arange(sp1.shape[0]))

    fig, ax = plt.subplots(2, 1, figsize=(10, 5))
    ax[0].plot(np.arange(sp1.shape[0]), middle_pixel_fp_peak, 'k.')
    ax[1].plot(np.arange(sp1.shape[0]), diff, 'r.')
    ax[1].plot(np.arange(sp1.shape[0])[keep], diff[keep], 'g.')
    ax[0].set_xlabel('Order')
    ax[0].set_ylabel('Cavity length')
    ax[1].set_xlabel('Order')
    ax[1].set_ylabel('Cavity length - fit')
    ax[1].grid()
    plt.show()

    ith_coeff = 0
    # --- 13c. For each order, compute the final wavelength solution and store coefficients in header ---
    for iord in (range(sp1.shape[0])):
        fp_pix, mu_pix, amp, peak_count = get_lines_pix(fp1[iord], fp=True)
        # because the Nth peak is decreasing
        peak_count = np.max(peak_count) - peak_count

        if np.isfinite(all_cheby[iord][0]):
            wave = val_cheby(all_cheby[iord], fp_pix, domain=[0, fp1.shape[1] - 1])
            cavity = np.polyval(fit_cavity, wave)
            int_fp = cavity / wave
        else:
            valid = np.isfinite(cavity_middle)
            ord = np.arange(sp1.shape[0])[valid]
            fit_cavity_to_order = np.polyfit(ord, cavity_middle[valid], 7)
            cavity = np.polyval(fit_cavity_to_order, iord)
            int_fp = peak_count

        middle_fp = ius(fp_pix, int_fp, k=1)(2044)

        err_fp_pos = middle_fp - np.polyval(fit_middle_fp, iord)
        off = np.round(err_fp_pos)

        int_fp = int_fp - off
        print(iord, err_fp_pos, off)

        int_fp = np.round(int_fp).astype(int)

        for ite in range(3):
            wave = cavity / int_fp
            cavity = np.polyval(fit_cavity, wave)

        all_cheby[iord] = fit_cheby(fp_pix, wave, CHEBY_FIT_DEG, [0, sp1.shape[1] - 1])
        for icoeff in range(WAVEDEGN + 1):
            hdr_wavesol[f'WAVE0{str(ith_coeff).zfill(3)}'] = fit_cheby_wave[icoeff], \
                f'Wavelength coefficients order={iord} coeffs={icoeff}'
            ith_coeff += 1
        final_wave_sol[iord] = val_cheby(all_cheby[iord], np.arange(sp1.shape[1]), domain=[0, sp1.shape[1] - 1])
    # plt.show()

    hdr_wavesol['WAVEORDN'] = N_ORDERS
    hdr_wavesol['WAVEDEGN'] = WAVEDEGN

    # --- 14. Save the final wavelength solution to a FITS file ---
    fits.writeto('spip_wave_sol.fits', final_wave_sol, header=hdr_wavesol, overwrite=True)

    # --- 15. Save the cavity fit coefficients to a file ---
    outname_cavity_fit = 'cavity_length_ll_fit.dat'
    with open(outname_cavity_fit, 'w') as f:
        for i in range(len(fit_cavity)):
            f.write(f"{fit_cavity[i] / 2}\n")

    # --- 16. Plot the final wavelength solution for all orders ---
    fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True)
    for iord in range(sp1.shape[0]):
        ax[0].plot(final_wave_sol[iord], sp1[iord], '-', alpha=0.5)
        ax[0].text(np.mean(final_wave_sol[iord]), 0, str(iord), fontsize=8)
        ax[1].plot(final_wave_sol[iord], fp1[iord], '-', alpha=0.5)
        ax[0].set_yscale('log')
        ax[1].set_yscale('log')
        meanwave_order = np.mean(final_wave_sol[iord])
        mean_log_flux = np.log10(np.nanmedian(sp1[iord][sp1[iord] > 0]))
        ax[0].text(meanwave_order, mean_log_flux, str(iord), fontsize=8)
    plt.show(block=True)