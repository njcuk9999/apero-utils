import glob
import os
# Import Union and List from typing
from typing import Union, List

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.constants import c
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from tqdm import tqdm

from etienne_tools import mjd_to_matplotlib_date, odd_ratio_mean
import getpass

whoami = getpass.getuser()

if whoami == 'eartigau':
    calib_dir = '/Users/eartigau/glitch_fp/data_NIRPS_HE/'
    output_slinky = '/Users/eartigau/glitch_fp/data_NIRPS_HE_output/'
    if not os.path.isdir(output_slinky):
        os.mkdir(output_slinky)

    doplot = True
    sync = False

if whoami == 'spirou':
    calib_dir = '/space/spirou/SLINKY/data_NIRPS_HE/'
    output_slinky = '/space/spirou/SLINKY/data_NIRPS_HE_output/'
    if not os.path.isdir(output_slinky):
        os.mkdir(output_slinky)

    doplot = False
    sync = True


def val_cheby(coeffs: np.ndarray, xvector: Union[np.ndarray, int, float],
              domain: List[float]) -> Union[np.ndarray, int, float]:
    """
    Using the output of fit_cheby calculate the fit to x  (i.e. y(x))
    where y(x) = T0(x) + T1(x) + ... Tn(x)

    :param coeffs: output from fit_cheby
    :param xvector: x value for the y values with fit
    :param domain: domain to be transformed to -1 -- 1. This is important to
    keep the components orthogonal. For SPIRou orders, the default is 0--4088.
    You *must* use the same domain when getting values with fit_cheby
    :return: corresponding y values to the x inputs
    """
    # transform to a -1 to 1 domain
    domain_cheby = 2 * (xvector - domain[0]) / (domain[1] - domain[0]) - 1
    # fit values using the domain and coefficients
    yvector = np.polynomial.chebyshev.chebval(domain_cheby, coeffs)
    # return y vector
    return yvector


def gp_project(x, y, yerr, wslinky=1e-1, xmin=980, xmax=1850, npts=100000):
    """
    Project data points onto a regular grid using a Gaussian weight.

    :param x: The x values for which we have data and errors
    :param y: The y values for which we have data and errors
    :param yerr: The error on y
    :param wslinky: The e-width of the Gaussian kernel
    :param xmin: The starting point of the grid
    :param xmax: The end point of the grid
    :param npts: The number of points in the grid
    :return: The x and y values of the grid at which we have projected the data
    """
    # Create a grid of x values
    xv = np.linspace(xmin, xmax, npts)
    # Initialize weights and y values for the grid
    weights = np.full(npts, 1e-12)
    yv = np.zeros(npts)

    xvbis = xv / wslinky
    xbis = x / wslinky
    # Loop over each data point
    for i in tqdm(range(len(x)), leave = False):
        # Calculate the distance between the grid points and the data point
        dd = xvbis - xbis[i]
        g = np.abs(dd) < 10

        dd2 = dd[g]

        # Calculate the weight of the data point
        w2 = np.exp(-0.5 * dd2 ** 2) / yerr[i] ** 2
        # Add the weight to the grid weights
        weights[g] += w2
        # Add the weighted y value to the grid y values
        yv[g] += w2 * y[i]
    # Normalize the y values by the weights
    yv /= weights

    return xv, yv


def odd_ratio_linfit(x, y, yerr):
    """
    Fit a linear model to the data using an iterative weighted least squares method.

    :param x: Abscissa
    :param y: Ordinate
    :param yerr: Error on the ordinate
    :return: Linear fit and error on the fit
    """
    # Remove NaN values
    g = np.isfinite(y + yerr + x)
    x = x[g]
    y = y[g]
    yerr = yerr[g]
    # Initialize weights
    w = np.ones(len(x))

    # Iterate until weights converge
    sum = 1.0
    sum0 = 0.0
    while np.abs(sum0 - sum) > 1e-6:
        sum0 = np.sum(w)
        # Fit the data with current weights
        fit, sig = np.polyfit(x, y, 1, w=w / yerr, cov=True)
        errfit = np.sqrt(np.diag(sig))
        # Compute residuals and update weights
        res = (y - np.polyval(fit, x)) / yerr
        p1 = np.exp(-0.5 * res ** 2)
        p2 = 1e-6
        w = p1 / (p1 + p2)
        sum = np.sum(w)

    return fit, errfit


def sigma(v):
    """
    Compute the standard deviation of the data using a robust estimator.

    :param v: Vector
    :return: Half-width of the 68\% confidence interval
    """
    n1, p1 = np.nanpercentile(v, [16, 84])
    return 0.5 * (p1 - n1)


def search_fits_with_mjd(search_string, mjdkey='MJDMID'):
    """
    Search for FITS files with a given search string and sort them by MJD.

    :param search_string: Search string
    :param mjdkey: MJD key in the header
    :return: Files sorted by MJD, MJDs sorted
    """
    files = glob.glob(search_string)
    mjds = np.zeros(len(files))
    for ifile, file in enumerate(files):
        h = fits.getheader(file)
        mjds[ifile] = h[mjdkey]
    order = np.argsort(mjds)
    files = np.array(files)[order]
    mjds = mjds[order]

    return files, mjds

# Leverage point of the linear fit, typically the RV content barycenter of the domain
wave_leverage = 1600
# Autocorrelation length of the (semi)GP kernel used to model the slinky effect
wslinky = 1e-1

# which instrument are we using
#
instrument = 'NIRPS_HE'
#instrument = 'SPIROU'


if not os.path.isdir(calib_dir):
    os.mkdir(calib_dir)

if instrument == 'SPIROU':
    maestria_path = 'spirou_offline'
    inst_short = 'spirou'
    fiber = 'AB'
    # FOR NIRPS_HE only, probably OK for NIRPS_HA, **not** of for SPIRou
    inst_wavestart = 965
    inst_waveend = 2500
elif instrument == 'NIRPS_HE':
    maestria_path = 'nirps_he_offline'
    inst_short = 'nirps'
    fiber = 'A'
    # FOR NIRPS_HE only, probably OK for NIRPS_HA, **not** of for SPIRou
    inst_wavestart = 965
    inst_waveend = 1950
elif instrument == 'NIRPS_HA':
    maestria_path = 'nirps_ha_offline'
    inst_short = 'nirps'
    fiber = 'A'
    # FOR NIRPS_HE only, probably OK for NIRPS_HA, **not** of for SPIRou
    inst_wavestart = 965
    inst_waveend = 1950
else:
    raise ValueError(f'Instrument {instrument} not recognized')

if sync:
    if whoami == 'artigau':
        # Sync files from the server to the local calibration directory
        cmd = f'rsync -avz spirou@maestria:/cosmos99/{inst_short}/apero-data/{maestria_path}/red/*/*lines_{fiber}.fits' \
              f' {calib_dir}'
        os.system(cmd)
        cmd = f'rsync -avz spirou@maestria:/cosmos99/{inst_short}/apero-data/{maestria_path}/calib/*_wavesol_ref_{fiber}.fits {calib_dir}'
        os.system(cmd)
        cmd = f'rsync -avz spirou@maestria:/cosmos99/{inst_short}/apero-data/{maestria_path}/calib/*pp_e2dsff_{fiber}_wave_night_{fiber}.fits {calib_dir}'
        os.system(cmd)

    elif whoami == 'spirou':
        # Sync files from the server to the local calibration directory
        cmd = f'rsync -avz maestria:/cosmos99/{inst_short}/apero-data/{maestria_path}/red/*/*lines_{fiber}.fits' \
              f' {calib_dir}'
        os.system(cmd)
        cmd = f'rsync -avz /cosmos99/{inst_short}/apero-data/{maestria_path}/calib/*_wavesol_ref_{fiber}.fits {calib_dir}'
        os.system(cmd)
        cmd = f'rsync -avz /cosmos99/{inst_short}/apero-data/{maestria_path}/calib/*pp_e2dsff_{fiber}_wave_night_{fiber}.fits {calib_dir}'
        os.system(cmd)

# Load the reference wave solution
wavesol = glob.glob(f'{calib_dir}/*_wavesol_ref_{fiber}.fits')[0]
ref_wave_sol = fits.getdata(wavesol)

# Find all FP and HC files and their corresponding MJDs
files_fp, mjds_fp = search_fits_with_mjd(f'{calib_dir}/*wave_fplines_{fiber}.fits')
files_hc, mjds_hc = search_fits_with_mjd(f'{calib_dir}/*wave_hclines_{fiber}.fits')

#files_hc = files_hc[::10]
#mjds_hc = mjds_hc[::10]

dt = np.zeros(len(files_fp))
for i_fp in range(len(files_fp)):
    dt[i_fp] = np.min(np.abs(mjds_fp[i_fp] - mjds_hc))
g = dt<0.1
files_fp = files_fp[g]
mjds_fp = mjds_fp[g]

# files_hc = np.array(files_hc)
# order_files = np.argsort(np.random.rand(len(files_hc)))

# files_hc = files_hc[order_files]
# mjds_hc = mjds_hc[order_files]


# Process each HC file
for i_hc in range(len(files_hc)):
    print('Processing file {}/{}'.format(i_hc + 1, len(files_hc)))
    hdr = fits.getheader(files_hc[i_hc])
    tbl_hc = Table.read(files_hc[i_hc], 'WAVE_HCLIST')

    if 'CAVITY' in tbl_hc.colnames:
        print('Cavity already computed, skipping')
        continue
    else:
        print('Computing cavity for', files_hc[i_hc])

    # Find the corresponding FP file
    i_fp = np.argmin(np.abs(mjds_fp - mjds_hc[i_hc]))
    delta_t = mjds_fp[i_fp] - mjds_hc[i_hc]

    print('Delta t =', delta_t)
    if np.abs(delta_t) > 0.5:
        print('No corresponding FP file found, skipping')
        continue

    tbl_fp = Table.read(files_fp[i_fp], 'WAVE_FPLIST')

    tbl_hc['CAVITY'] = np.nan
    current_order = -1

    mask = tbl_hc['PIXEL_MEAS'].mask

    all_frac_peak = np.zeros_like(tbl_hc['PIXEL_MEAS'])
    all_cavity = np.zeros_like(tbl_hc['PIXEL_MEAS'])
    wave_ref = np.array(tbl_hc['WAVE_REF'])
    pixel_meas = np.array(tbl_hc['PIXEL_MEAS'])

    for i in tqdm(range(len(tbl_hc)), leave=False):
        if mask[i]:
            continue

        order = tbl_hc['ORDER'][i]
        if order != current_order:
            tbl_fp_order = tbl_fp[tbl_fp['ORDER'] == order]
            tbl_fp_order = tbl_fp_order[~tbl_fp_order['PIXEL_MEAS'].mask]
            current_order = order

            spl = ius(tbl_fp_order['PIXEL_MEAS'], tbl_fp_order['PEAK_NUMBER'], k=1, ext=1)

        all_frac_peak[i] = spl(pixel_meas[i])  # np.polyval(np.polyfit(pixel_meas2, peak_number2, 2), pixel_meas)
        all_cavity[i] = wave_ref[i] * all_frac_peak[i]

    bad = all_frac_peak == 0
    all_frac_peak[bad] = np.nan
    all_cavity[bad] = np.nan

    tbl_hc['PEAK_NUMBER'] = all_frac_peak
    tbl_hc['CAVITY'] = all_cavity

    # Update the FITS file with the new table
    for col in tbl_hc.colnames:
        try:
            tbl_hc[col][tbl_hc[col].mask] = np.nan
        except:
            pass

    with fits.open(files_hc[i_hc]) as hdul:
        tbl_hc = tbl_hc.as_array()
        hdul[1].data = tbl_hc
        hdul.writeto(files_hc[i_hc], overwrite=True)

# Load a reference table and header
ref_file_hc = search_fits_with_mjd(f'{calib_dir}/*waveref_hclines_{fiber}.fits')[0][0]
tbl1 = Table.read(ref_file_hc, 'WAVE_HCLIST_REF')
href = fits.getheader(ref_file_hc)
cavity_polynomial = np.array([href[key] for key in href['WCAV0*'].keys()])
WCAV_PED = href['WCAV_PED']

# Find all wave HC files
order = np.argsort(mjds_hc)
mjds_hc = mjds_hc[order]
files_hc = files_hc[order]

# Initialize an array to store cavity values
all_cavity = np.zeros([len(files_hc), len(tbl1)], dtype=float)
for ifile, file in tqdm(enumerate(files_hc), leave=False):
    all_cavity[ifile] = Table.read(file, 'WAVE_HCLIST')['CAVITY'].data.data

# Replace zero values with NaN
all_cavity[all_cavity == 0] = np.nan
# Compute the median cavity per line
med_per_line = np.nanmedian(all_cavity, axis=0)

meds = np.zeros(len(files_hc))
# Normalize the cavity values
for iepoch in range(len(files_hc)):
    meds[iepoch] = np.nanmedian(all_cavity[iepoch] - med_per_line)
    all_cavity[iepoch] -= meds[iepoch]
plt.plot(mjds_hc, meds, '.')
plt.show()

# Update the median cavity per line
med_per_line = np.nanmedian(all_cavity, axis=0)

cavity_ref = val_cheby(cavity_polynomial, tbl_hc['WAVE_REF'], domain=[inst_wavestart, inst_waveend])+WCAV_PED
dv_ref = c*(med_per_line/cavity_ref-1)
bad = np.abs(dv_ref)>1000
dv_ref[bad] = np.nan
# we offset the reference cavity to match the median cavity
tbl1['WAVE_REF'] = tbl1['WAVE_REF']*(1-dv_ref/c)

med_per_line[bad] = np.nan


n1, p1 = np.nanpercentile(all_cavity, [16, 84], axis=0)
sig_per_line = (p1 - n1) / 2

bad = np.sum(np.isfinite(all_cavity), axis=0) < all_cavity.shape[0] // 5
sig_per_line[sig_per_line == 0] = np.nan
sig_per_line[bad] = np.nan

sig_per_line_ms = c * (sig_per_line / med_per_line)

# Initialize arrays to store fit parameters
all_slopes = np.zeros_like(files_hc, dtype=float)
all_errslopes = np.zeros_like(files_hc, dtype=float)
all_pedestals = np.zeros_like(files_hc, dtype=float)
all_errpedestals = np.zeros_like(files_hc, dtype=float)

if doplot:
    # Create a figure for plotting with 2 rows and 1 column, sharing the x-axis
    fig, ax = plt.subplots(ncols=1, nrows=2, sharex=True)
    # Initialize the plots with NaN values to set up the date format
    ax[0].plot_date([np.nan], [np.nan])
    ax[1].plot_date([np.nan], [np.nan])

# Loop through each wave HC file
for ifile, file in enumerate(files_hc):
    # Read the table from the current file
    tbl2 = Table.read(file, 'WAVE_HCLIST')

    # Skip files that do not have the 'CAVITY' column
    if 'CAVITY' not in tbl2.colnames:
        continue
    # Get the header of the current file
    h = fits.getheader(file)

    # Calculate the cavity deviation from the median per line, normalized by WCAV_PED and converted to m/s
    dcavity = c * (tbl2['CAVITY'].data.data / med_per_line - 1)  # / WCAV_PED
    sdcavity = c * sig_per_line / med_per_line
    # Get the NSIG values from the table
    nsig = tbl2['NSIG'].data.data
    # Normalize the WAVE_REF values and shift by the leverage point
    wave2 = tbl2['WAVE_REF'] / 1e3 - wave_leverage / 1e3
    # Initialize the fit parameters
    fit = [1, 0]
    sig_fit = [1, 1]

    # fit[0],sig_fit[0] = odd_ratio_mean(dcavity, sdcavity)
    moy, err = odd_ratio_mean(dcavity, sdcavity)
    fit[1] = moy
    sig_fit[1] = err
    slope, err_slope = odd_ratio_mean((dcavity - moy) / wave2, sdcavity / wave2)
    fit[0] = slope
    sig_fit[0] = err_slope

    # Iteratively fit the data to refine the fit parameters
    # Calculate the ratio of the standard deviation of the residuals
    # ratio = sigma((dcavity - np.polyval(fit, wave2)) * nsig)
    # Calculate the scaled standard deviation of the cavity
    # sdcavity = ratio / nsig
    # Remove outliers by setting large deviations to NaN
    # dcavity[np.abs(dcavity-np.nanmedian(dcavity)) > 1000] = np.nan
    # Perform a weighted linear fit to the data
    fit, sig_fit = odd_ratio_linfit(wave2, dcavity - moy, sdcavity)
    fit[1] += moy

    # Print the fit results for the current file
    args = (ifile + 1, len(files_hc), fit[1], sig_fit[1], fit[0], sig_fit[0])
    print('\t{}/{} zp {:5.2f}+-{:5.2f} [m/s], slope {:5.2f}+-{:5.2f} [m/s/µm]'.format(*args))

    if doplot:
        # Plot the zero-point and slope with error bars
        ax[0].errorbar(mjd_to_matplotlib_date(h['MJDMID']), fit[1], yerr=sig_fit[1], fmt='.g')

        ax[1].errorbar(mjd_to_matplotlib_date(h['MJDMID']), fit[0], yerr=sig_fit[0], fmt='.g')

    # rotate at 45 deg
    # ax[0].text(mjd_to_matplotlib_date(h['MJDMID']), fit[1], file.split('/')[-1].split('_')[0], fontsize = 8,
    # rotation = 45)

    # Store the fit parameters and their errors
    all_slopes[ifile] = fit[0]
    all_errslopes[ifile] = sig_fit[0]
    all_pedestals[ifile] = fit[1]
    all_errpedestals[ifile] = sig_fit[1]
    # Store the MJD of the current file

if doplot:
    ax[0].grid(color='grey', linestyle='--', linewidth=0.5)
    ax[1].grid(color='grey', linestyle='--', linewidth=0.5)
    ax[0].set(ylabel='Zero-point [m/s]')
    ax[1].set(ylabel='Slope [m/s/µm]', xlabel='MJD')
    ax[1].set(xlabel='Date')
    # Display the plot
    plt.show()

# Print the standard deviation of the differences between consecutive pedestals
print(sigma(all_pedestals - np.roll(all_pedestals, 1)) / np.sqrt(2))


recovered_pedestal = np.zeros_like(all_pedestals)
recovered_slope = np.zeros_like(all_slopes)
recovered_errslope = np.zeros_like(all_errslopes)
recovered_errpedestal = np.zeros_like(all_errpedestals)

# Loop through each FP file
for i_fp in range(len(files_fp)):
    print('Processing wavelength solution of night {}/{}'.format(i_fp + 1, len(files_fp)))
    file_fp = files_fp[i_fp]
    hdr = fits.getheader(file_fp)
    wavefile = f'{calib_dir}' + hdr['WAVEFILE']

    # Check if the wavefile exists
    #if not os.path.exists(wavefile):
    #    print('No wavefile found, skipping')
    #    continue

    hdr = fits.getheader(wavefile)

    # Skip if the wavelength solution is already corrected for the cavity effect
    #if False:
    #    # if 'SLINKY' in hdr.keys():
    #    print('Wavelength solution already corrected for cavity effect, skipping')
    #    continue

    tbl_fp = Table.read(file_fp, 'WAVE_FPLIST')
    hdr_fp = fits.getheader(file_fp)

    # Find the corresponding HC (Hollow Cathode) file based on the closest MJD
    i_hc = np.argmin(np.abs(mjds_hc - hdr_fp['MJDMID']))
    slope_hc = all_slopes[i_hc]
    err_slope_hc = all_errslopes[i_hc]
    pedestal_hc = all_pedestals[i_hc]
    err_pedestal_hc = all_errpedestals[i_hc]

    # Apply the polynomial cavity correction to the wavelength model
    wavelength_model = np.array(tbl_fp['WAVE_REF'].data)
    doppler = (wavelength_model / 1e3 - wave_leverage / 1e3) * slope_hc + pedestal_hc
    #doppler = 0

    for ite in range(4):
        peak_number = tbl_fp['PEAK_NUMBER'].data
        cavity = val_cheby(cavity_polynomial, wavelength_model, domain=[inst_wavestart, inst_waveend]) + WCAV_PED


        # This is wrong and leads to a big error in RV
        # wavelength_model = cavity / peak_number * (1 - doppler / c)
        wavelength_model = cavity / peak_number * (1 + doppler / c)

    # diff = c*(wavelength_model / tbl_fp['WAVE_REF'].data-1)
    # fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = True)
    # ax[0].plot(tbl_fp['WAVE_REF'], tbl_fp['WAVE_REF'] * tbl_fp['PEAK_NUMBER'])
    tbl_fp['WAVE_REF'] = wavelength_model
    # ax[0].plot(tbl_fp['WAVE_REF'], tbl_fp['WAVE_REF'] * tbl_fp['PEAK_NUMBER'],'r.')
    # ax[1].plot(tbl_fp['WAVE_REF'], diff)
    # plt.show()

    # Initialize the wavelength map with the reference wave solution
    wavemap = np.array(ref_wave_sol)

    dv_residuals = []
    dv_residuals_err = []
    wave_ref_residuals = []

    # Loop through each order in the table
    for order in np.unique(tbl_fp['ORDER']):
        g = tbl_fp['ORDER'] == order
        tbl_order = tbl_fp[g]
        pixel_meas = np.array(tbl_order['PIXEL_MEAS'])
        wave_ref = np.array(tbl_order['WAVE_REF'])
        valid = np.isfinite(pixel_meas + wave_ref)
        pixel_meas = pixel_meas[valid]
        wave_ref = wave_ref[valid]

        # Fit a polynomial to the pixel measurements and wave references
        fit = np.polyfit(pixel_meas, wave_ref, 5)
        residual = (wave_ref / np.polyval(fit, pixel_meas) - 1) * c
        wave_order = np.polyval(fit, np.arange(4088))

        # Update the wavelength map if the wave order is valid
        if np.all(np.isfinite(wave_order)):
            wavemap[order] = np.polyval(fit, np.arange(4088))

        # Calculate the error on the residuals
        err = sigma(residual - np.roll(residual, 1)) / np.sqrt(2)
        err = np.ones_like(residual) * err

        dv_residuals.append(residual)
        dv_residuals_err.append(err)
        wave_ref_residuals.append(wave_ref)

    # Concatenate the residuals and errors
    dv_residuals = np.concatenate(dv_residuals)
    dv_residuals_err = np.concatenate(dv_residuals_err)
    wave_ref_residuals = np.concatenate(wave_ref_residuals)

    # We have 5 points per scale length of slinky
    npts = int((np.nanmax(wavemap) - np.nanmin(wavemap)) / wslinky * 5)
    xmin = np.nanmin(wavemap)
    xmax = np.nanmax(wavemap)
    # Project the residuals onto a regular grid using a Gaussian weight
    xv, yv = gp_project(wave_ref_residuals, dv_residuals, dv_residuals_err, wslinky=wslinky, xmin=xmin, xmax=xmax,
                        npts=npts)
    spl = ius(xv, yv, k=2)

    # Correct the wavelength map for the cavity effect
    wavemap = wavemap * (1 + spl(wavemap) / c)

    # this is a sanity check that the wavelength solution is now consistent. We backproject the HC line position
    # onto the wavelength grid and compare with the nominal HC line position

    i_hc = np.argmin(np.abs(hdr_fp['MJDMID'] - mjds_hc))
    tbl_hc = Table.read(files_hc[i_hc])
    tbl_hc['WAVE_REF'] = tbl1['WAVE_REF']

    v_hc = []
    v_fp = []
    hc_meas = np.array(tbl_hc['PIXEL_MEAS'])

    if doplot:
        fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)

    for iord in np.unique(tbl_hc['ORDER']):
        g = tbl_hc['ORDER'] == iord
        wave_ord = wavemap[iord]
        pix = np.arange(4088)
        spl = ius(pix, wave_ord, k=2)
        # tbl_hc['WAVE_MEAS'][g] = spl(tbl_hc['PIXEL_MEAS'][g])

        hc_meas[g] = spl(tbl_hc['PIXEL_MEAS'][g])

        v = np.array(tbl_hc['WAVE_REF'][g] / hc_meas[g] - 1) * c
        v_hc.append(v)

        # print(iord, np.nanmedian(v))
        if doplot:
            ax[0].plot(tbl_hc['WAVE_REF'][g], v, '.')

        g = tbl_fp['ORDER'] == iord

        v = np.array(tbl_fp['WAVE_REF'][g] / spl(tbl_fp['PIXEL_MEAS'][g]) - 1) * c

        if doplot:
            ax[1].plot(tbl_fp['WAVE_REF'][g], v, '.')
        v_fp.append(v)

    if doplot:
        # ax[0].plot(tbl_hc['WAVE_REF'],c*(tbl_hc['WAVE_MEAS']/tbl_hc['WAVE_REF']-1),'.')
        ax[0].set(xlabel='Nominal HC line position [nm]')
        ax[1].set(xlabel='Nominal FP line position [nm]')
        ax[0].set(ylabel='RV [m/s]')
        ax[1].set(ylabel='RV [m/s]')
        ax[0].set(ylim=[-100, 100])
        ax[1].set(ylim=[-100, 100])
        plt.show()

    dv_hc =  np.array( (hc_meas / tbl_hc['WAVE_REF']-1)*c)

    valid = np.isfinite(dv_hc)
    fit, sig = odd_ratio_linfit((tbl_hc['WAVE_REF'] - wave_leverage)[valid]/1e3, dv_hc[valid], sig_per_line_ms[valid])

    recovered_pedestal[i_fp] = fit[1]
    recovered_slope[i_fp] = fit[0]
    recovered_errslope[i_fp] = sig[0]
    recovered_errpedestal[i_fp] = sig[1]

    print('{} -> zp {:5.2f}+-{:5.2f} [m/s], slope {:5.2f}+-{:5.2f} [m/s/µm]'.format(i_fp,fit[1], sig[1],fit[0],
                                                                                    sig[0]))

    v_hc = np.concatenate(v_hc)
    v_fp = np.concatenate(v_fp)

    print('HC', sigma(v_hc))
    print('FP', sigma(v_fp))

    # Update the wavefile with the corrected wavelength map
    with fits.open(wavefile) as hdul:
        hdul[1].data = wavemap
        hdul[0].header['SLINKY'] = True, 'Wavelength solution corrected for cavity effect'
        hdul.writeto(wavefile, overwrite=True)

print('RMS : {:.3f} m/s pedestal'.format(np.std(recovered_pedestal)))
print('RMS : {:.3f} m/s/µm slope'.format(np.std(recovered_slope)))

if doplot:
    fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
    ax[0].errorbar(mjds_hc, recovered_pedestal, yerr=recovered_errpedestal, fmt='g.')
    ax[1].errorbar(mjds_hc, recovered_slope, yerr=recovered_errslope, fmt='g.')
    ax[0].set(ylabel='Zero-point [m/s]')
    ax[1].set(ylabel='Slope [m/s/µm]', xlabel='MJD')
    ax[1].set(xlabel='Date')
    plt.show()