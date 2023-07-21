"""
Script to plot the extracted and telluric corrected apero reduced spectra of a e2dsff file
Optionally, it can also plot the geneva reduced spectra
Currently it only shows the plots and doesn't save it
This script creates folders named apero_data and geneva_data where it stores the downloaded data

How to use:
	Find a the path to a s2d apero spectra (should end with _pp_e2dsff_A.fits)
	Put this path in the ext_data_path variable
	Select the (apero) order that you want to plot and put it in order_num
    Specify the standard deviation of the gaussian filter applied to the data (gaussian_filter_sigma), put 0 to ignore the gaussian filter
	Specify the plot title in the variable plot_title
	Tell whether or not to also plot the geneva reduced spectra in the variable include_geneva
	Run python plot_ext_tcorr_s2d.py

@author: Luc Bazinet
"""

import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d


########## Only modify these parameters #############
ext_data_path = '/cosmos99/nirps/apero-data/nirps_he_07276_online/red/2023-04-02/NIRPS_2023-04-03T04_04_01_090_pp_e2dsff_A.fits'
order_num = 10 # apero order
gaussian_filter_sigma = 2
include_geneva = True
plot_title = f'WASP-178\n{ext_data_path.split("/")[-1]}\norder {order_num}'
#####################################################


def read_NIRPS_wave(header, order):
    """
    Extracts the waveloengths of a specific order given the header of an apero fits file
    """
    nord = header['WAVEDEGN']
    npix = 4088
    domain = [0, npix]
    coeffs = []
    for i in range(order * (nord + 1), (order + 1) * (nord + 1)):
        coeffs.append(header[f'WAVE{i:04d}'])
    coeffs = np.array(coeffs)
    domain_cheby = 2 * (np.arange(npix) - domain[0]) / (domain[1] - domain[0]) - 1
    yvector = np.polynomial.chebyshev.chebval(domain_cheby, coeffs)
    return yvector


def rsync_nirpsclient(location, data_folder):
    """
    Downloads file from nirps-client@maestria.astro.umontreal.ca:location where location is the argument given in this method
    The downloaded file will be in data_folder
    """
    subprocess.run(['rsync', '-avu', 'nirps-client@maestria.astro.umontreal.ca:' + location, data_folder])


if not ext_data_path.endswith('pp_e2dsff_A.fits'):
    raise Exception("File path doesn't end with 'pp_e2dsff_A.fits'. Make sure that the file is a e2dsff APERO file.")

tcorr_data_path = ext_data_path.replace('e2dsff_A', 'e2dsff_tcorr_A')
obs_date, ext_filename = ext_data_path.split('/')[-2:]
tcorr_filename = tcorr_data_path.split('/')[-1]

if not os.path.exists('apero_data/'):
    os.makedirs('apero_data/')

if not os.path.exists('apero_data/ext/' + ext_filename):
    rsync_nirpsclient(ext_data_path, 'apero_data/ext/')

if not os.path.exists('apero_data/tcorr/' + tcorr_filename):
    rsync_nirpsclient(tcorr_data_path, 'apero_data/tcorr/')

ext_fits = fits.open('apero_data/ext/' + ext_filename)
tcorr_fits = fits.open('apero_data/tcorr/' + tcorr_filename)

wave = read_NIRPS_wave(ext_fits[0].header, order_num)
ext_flux = ext_fits[1].data[order_num]
tcorr_flux = tcorr_fits[1].data[order_num]

if gaussian_filter_sigma > 0:
    ext_flux = gaussian_filter1d(ext_flux, gaussian_filter_sigma)
    tcorr_flux = gaussian_filter1d(tcorr_flux, gaussian_filter_sigma)

plt.plot(wave, ext_flux, linewidth=0.7, label='Apero extracted spectra')
plt.plot(wave, tcorr_flux, linewidth=0.7, label='Apero telluric corrected spectra')

if include_geneva:

    # The geneva drs remove heavily telluric contaminated orders, this next array is to correspond the apero order to the geneva order
    # This is not robust and a better way to equate apero to geneva order is welcomed
    apero_to_geneva_orders = np.array([-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                                       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
                                       41, 42, -1, -1, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 56, 57,
                                       58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, -1])
    geneva_order = apero_to_geneva_orders[order_num]
    if geneva_order == -1:
        print(f'Order {order_num} is not in Geneva drs')
    else:
        if not os.path.exists('geneva_data/'):
            os.makedirs('geneva_data/')
        filedate = ':'.join(ext_filename.split('_')[1:4]) + '.' + ext_filename.split('_')[4]
        geneva_ext_filename = f'r.NIRPS.{filedate}_S2D_BLAZE_A.fits'
        geneva_tcorr_filename = geneva_ext_filename.replace('S2D_BLAZE_A', 'S2D_BLAZE_TELL_CORR_A')
        if not os.path.exists('geneva_data/ext/' + geneva_ext_filename):
            rsync_nirpsclient(f'/cosmos99/nirps/geneva-data/DRS-3.0.0/reduced/{obs_date}/{geneva_ext_filename}', 'geneva_data/ext/')
        if not os.path.exists('geneva_data/tcorr/' + geneva_tcorr_filename):
            rsync_nirpsclient(f'/cosmos99/nirps/geneva-data/DRS-3.0.0/reduced/{obs_date}/{geneva_tcorr_filename}', 'geneva_data/tcorr/')

        geneva_ext_fits = fits.open('geneva_data/ext/' + geneva_ext_filename)
        geneva_tcorr_fits = fits.open('geneva_data/tcorr/' + geneva_tcorr_filename)

        geneva_ext_data_cube = np.array(geneva_ext_fits[1].data.tolist())
        geneva_wave = geneva_ext_fits[4].data[geneva_order] / 10 # convert Angstrom to nm
        
        c_kms = 299792.458
        berv = geneva_ext_fits[0].header['HIERARCH ESO QC BERV']
        berv_doppler_factor = np.sqrt((1+berv/c_kms)/(1-berv/c_kms))
        geneva_wave = geneva_wave / berv_doppler_factor # revert the BERV correction that the geneva pipeline does
        
        geneva_ext_flux = geneva_ext_fits[1].data[geneva_order]
        geneva_tcorr_flux = geneva_tcorr_fits[1].data[geneva_order]
        geneva_tcorr_flux[np.where(geneva_tcorr_flux == 0)] = np.nan

        if gaussian_filter_sigma > 0:
            geneva_ext_flux = gaussian_filter1d(geneva_ext_flux, gaussian_filter_sigma)
            geneva_tcorr_flux = gaussian_filter1d(geneva_tcorr_flux, gaussian_filter_sigma)

        plt.plot(geneva_wave, geneva_ext_flux, linewidth=0.7, label='Geneva extracted spectra')
        plt.plot(geneva_wave, geneva_tcorr_flux, linewidth=0.7, label='Geneva telluric corrected spectra')

plt.ylabel('Flux')
plt.xlabel('Wavelength [nm]')
plt.title(plot_title)
plt.legend()
plt.show()
