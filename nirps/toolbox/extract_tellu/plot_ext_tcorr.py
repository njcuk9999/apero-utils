'''
Script to plot the extracted and telluric corrected apero reduced spectra of a given target, date and exposure
Optionally, it can also plot the geneva reduced spectra
Currently it only shows the plots and doesn't save it
This script creates folders apero_data and geneva_data where it stores the downloaded data

Usage:
	Find a the path to a s1d apero spectra (should end with _pp_s1d_v_A.fits)
	Put this path in the ext_data_path variable
	Tell whether or not to also plot the geneva reduced spectra in the variable include_geneva
	Run python plot_ext_tcorr.py

@author: Luc Bazinet
'''

import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pdb


def rsync_nirpsclient(location, data_folder):
	'''
	Downloads file from nirps-client@maestria.astro.umontreal.ca:location where location is the argument given in this method
	The downloaded file will be in data_folder
	'''
	subprocess.run(['rsync', 'nirps-client@maestria.astro.umontreal.ca:' + location ,  data_folder])



ext_data_path = '/cosmos99/nirps/apero-data/nirps_he_07276_online/red/2023-04-30/NIRPS_2023-05-01T03_05_50_459_pp_s1d_v_A.fits'
plot_title = 'Tau Boo\nNIRPS_2023-05-01T03_05_50_459'
include_geneva = True



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


factor = np.nanpercentile(ext_flux, 50)
plt.plot(wave, ext_flux, linewidth=0.7, label='Apero extracted spectra')
plt.plot(wave, tcorr_flux, linewidth=0.7, label='Apero telluric corrected spectra')

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
	geneva_wave = geneva_ext_data_cube[:, 0]  # in Angstrom
	geneva_ext_flux = geneva_ext_data_cube[:, 2]
	geneva_tcorr_flux = np.array(geneva_tcorr_fits[1].data.tolist())[:, 2]

	geneva_factor = np.nanpercentile(geneva_ext_flux, 50)
	plt.plot(geneva_wave/10, geneva_ext_flux*(factor/geneva_factor), linewidth=0.7, label='Geneva extracted spectra')
	plt.plot(geneva_wave/10, geneva_tcorr_flux*(factor/geneva_factor), linewidth=0.7, label='Geneva extracted spectra')

	#plt.plot(geneva_wave/10, geneva_ext_flux, linewidth=0.7, label='Geneva extracted spectra')
	#plt.plot(geneva_wave/10, geneva_tcorr_flux, linewidth=0.7, label='Geneva telluric corrected spectra')

plt.ylabel('Flux')
plt.xlabel('Wavelength [nm]')
plt.title(plot_title)
plt.legend()
plt.show()

#pdb.set_trace()