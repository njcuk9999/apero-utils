"""
Script to plot the extracted and telluric corrected apero reduced spectra of a s1d file
Optionally, it can also plot the geneva reduced spectra
Currently it only shows the plots and doesn't save it
This script creates folders named apero_data and geneva_data where it stores the downloaded data

How to use:
	Find a the path to a s1d apero spectra (should end with _pp_s1d_v_A.fits)
	Put this path in the ext_data_path variable
	Specify the x limits of the plotting (xlim), use None to plot the whole data in that direction (e.g [None, 1360] will plot every lower than 1360 nm)
	Specify the standard deviation of the gaussian filter applied to the data (gaussian_filter_sigma), put 0 to ignore the gaussian filter
	Tell whether or not to also plot the geneva reduced spectra in the variable include_geneva
	Run python plot_ext_tcorr.py

@author: Luc Bazinet
"""

import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d
import pdb


########## Only modify these parameters #############
ext_data_path = '/cosmos99/nirps/apero-data/nirps_he_07276_online/red/2023-04-02/NIRPS_2023-04-03T04_04_01_090_pp_s1d_v_A.fits'
xlim = [None, None] # in nm
gaussian_filter_sigma = 0
include_geneva = True
plot_title = f'WASP-178\n{ext_data_path.split("/")[-1]}'
#####################################################


def rsync_nirpsclient(location, data_folder):
	"""
	Downloads file from nirps-client@maestria.astro.umontreal.ca:location where location is the argument given in this method
	The downloaded file will be in data_folder
	"""
	subprocess.run(['rsync', '-avu', 'nirps-client@maestria.astro.umontreal.ca:' + location,  data_folder])


if not ext_data_path.endswith('pp_s1d_v_A.fits'):
	raise Exception("File path doesn't end with 'pp_s1d_v_A.fits'. Make sure that the file is a s1d APERO file.")

tcorr_data_path = ext_data_path.replace('s1d_v_A', 's1d_v_tcorr_A')
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

ext_data_cube = np.array(ext_fits[1].data.tolist())
wave = ext_data_cube[:, 0]
ext_flux = ext_data_cube[:, 1]
tcorr_flux = np.array(tcorr_fits[1].data.tolist())[:, 1]

if xlim[0] is None:
	if xlim[1] is None:
		where = np.arange(len(wave))
	else:
		where = np.where(wave < xlim[1])
elif xlim[1] is None:
	where = np.where(wave > xlim[0])
else:
	where = np.where(np.logical_and(wave > xlim[0], wave < xlim[1]))
ext_flux_cut = ext_flux[where]
tcorr_flux_cut = tcorr_flux[where]
wave_cut = wave[where]

if gaussian_filter_sigma > 0:
	ext_flux_cut = gaussian_filter1d(ext_flux_cut, gaussian_filter_sigma)
	tcorr_flux_cut = gaussian_filter1d(tcorr_flux_cut, gaussian_filter_sigma)

factor = np.nanpercentile(ext_flux_cut, 50)

plt.plot(wave_cut, ext_flux_cut, linewidth=0.7, label='Apero extracted spectra')
plt.plot(wave_cut, tcorr_flux_cut, linewidth=0.7, label='Apero telluric corrected spectra')
plt.xlim(xlim)

if include_geneva:
	if not os.path.exists('geneva_data/'):
		os.makedirs('geneva_data/')

	filedate = ':'.join(ext_filename.split('_')[1:4]) + '.' + ext_filename.split('_')[4]
	geneva_ext_filename = f'r.NIRPS.{filedate}_S1D_A.fits'
	geneva_tcorr_filename = geneva_ext_filename.replace('S1D_A', 'S1D_SKYSUB_A')
	if not os.path.exists('geneva_data/ext/' + geneva_ext_filename):
		rsync_nirpsclient(f'/cosmos99/nirps/geneva-data/DRS-3.0.0/reduced/{obs_date}/{geneva_ext_filename}', 'geneva_data/ext/')
	if not os.path.exists('geneva_data/tcorr/' + geneva_tcorr_filename):
		rsync_nirpsclient(f'/cosmos99/nirps/geneva-data/DRS-3.0.0/reduced/{obs_date}/{geneva_tcorr_filename}', 'geneva_data/tcorr/')

	geneva_ext_fits = fits.open('geneva_data/ext/' + geneva_ext_filename)
	geneva_tcorr_fits = fits.open('geneva_data/tcorr/' + geneva_tcorr_filename)

	geneva_ext_data_cube = np.array(geneva_ext_fits[1].data.tolist())
	geneva_wave = geneva_ext_data_cube[:, 0]/10  # Angstrom to nm conversion
	geneva_ext_flux = geneva_ext_data_cube[:, 2]
	geneva_tcorr_flux = np.array(geneva_tcorr_fits[1].data.tolist())[:, 2]

	if xlim[0] is None:
		if xlim[1] is None:
			where = np.arange(len(geneva_wave))
		else:
			where = np.where(geneva_wave < xlim[1])
	elif xlim[1] is None:
		where = np.where(geneva_wave > xlim[0])
	else:
		where = np.where(np.logical_and(geneva_wave > xlim[0], geneva_wave < xlim[1]))
	geneva_ext_flux_cut = geneva_ext_flux[where]
	geneva_tcorr_flux_cut = geneva_tcorr_flux[where]
	geneva_wave_cut = geneva_wave[where]
	
	if gaussian_filter_sigma > 0:
		geneva_ext_flux_cut = gaussian_filter1d(geneva_ext_flux_cut, gaussian_filter_sigma)
		geneva_tcorr_flux_cut = gaussian_filter1d(geneva_tcorr_flux_cut, gaussian_filter_sigma)

	geneva_factor = np.nanpercentile(geneva_ext_flux_cut, 50)

	plt.plot(geneva_wave_cut, geneva_ext_flux_cut*(factor/geneva_factor), linewidth=0.7, label='Geneva extracted spectra')
	plt.plot(geneva_wave_cut, geneva_tcorr_flux_cut*(factor/geneva_factor), linewidth=0.7, label='Geneva sky sub spectra')

plt.ylabel('Flux')
plt.xlabel('Wavelength [nm]')
plt.title(plot_title)
plt.legend()
plt.show()

#pdb.set_trace()