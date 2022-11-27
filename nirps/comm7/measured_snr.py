"""
Code to tell us the maximum flux in a given set of files
and related this to the ND

auther: cook
date: 2022-11-23
"""
import glob
import os
from typing import List

import apero.core.math as mp
import numpy as np
from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

# ==================================================================
# VARIABLES
# ==================================================================
# choose which mode to use (either HA or HE)
modes = ['HA', 'HE']
# ------------------------------------------------------------------
# file string
FILE_STR = '*e2dsff_A.fits'
# set up parameters for HE and HA
params = dict(HE=dict(), HA=dict())
params['HA']['PATH'] = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202211/red/'
params['HE']['PATH'] = '/nirps_raw/nirps/apero-data/drs-data/nirps_he_202211/red/'

# header keys
KW_OBJNAME = 'OBJECT'
KW_SPT = 'HIERARCH ESO OCS TARG SPTYPE'
KW_DATE = 'DATE'
KW_EXPTIME = 'EXPTIME'
KW_AIRMASS = 'HIERARCH ESO TEL AIRM START'
KW_SEEING = 'HIERARCH ESO TEL AMBI FWHM START'
KW_MODE = 'HIERARCH ESO INS MODE'

# define orders to measure SNR in
orders = [19, 29, 39, 49, 59]
echelle_orders = [130, 120, 110, 100, 90]
# output tables
OUTTABLE_CSV = '/nirps_raw/nirps/misc/comm7/snr_measurement/all_results.csv'
OUTTABLE_FITS = '/nirps_raw/nirps/misc/comm7/snr_measurement/all_results.fits'
OUTTABLE_TXT = '/nirps_raw/nirps/misc/comm7/snr_measurement/all_results.txt'


# ==================================================================
# CLASSES
# ==================================================================
class NIRPSException(Exception):
	pass


# ==================================================================
# FUNCTIONS
# ==================================================================
def consecutive_spectra(files: np.ndarray) -> tuple[List[str], List[str]]:
	# storage for output
	first_file = []
	second_file = []
	# print progress
	print('Finding consecutive spectra of the same object')
	# loop around all files except the last one
	for it, filename1 in tqdm(enumerate(files[:-1])):
		# get the second filename
		filename2 = files[it + 1]
		# load header filename 1
		hdr1 = fits.getheader(filename1)
		objname1 = hdr1[KW_OBJNAME]
		# load header filename 2
		hdr2 = fits.getheader(filename2)
		objname2 = hdr2[KW_OBJNAME]
		# do not compare different modes
		if hdr1[KW_MODE] != hdr2[KW_MODE]:
			continue
		# if the next file has the same object name we say they are consecutive
		if objname1 == objname2:
			first_file.append(filename1)
			second_file.append(filename2)
	# return two list of files
	return first_file, second_file

# ==================================================================
# START OF CODE
# ==================================================================
if __name__ == '__main__':

	# get all files in path
	files = []
	for mode in modes:
		# choose parameters for instrument
		iparams = params[mode]
		# get files
		files += glob.glob(os.path.join(iparams['PATH'], '*', FILE_STR))

	# sort files by name
	files = np.sort(files)

	# get consecutive frames
	files1, files2 = consecutive_spectra(files)

	# column storage
	sdict = dict()
	sdict['target'] = []
	sdict['file1'], sdict['file2'] = [], []
	sdict['spt'], sdict['date1'], sdict['date2'] = [], [], []
	sdict['mode'], sdict['Texp1'], sdict['Texp1'] = [], [], []
	sdict['AM1'], sdict['seeing1'] = [], []
	sdict['AM2'], sdict['seeing2'] = [], []
	for order_num in echelle_orders:
		sdict[f'measureSNR1[{order_num}]'] = []
		sdict[f'measureSNR2[{order_num}]'] = []

	# print process
	print('Calculating SNR')
	# work out the snr
	for filename1, filename2 in tqdm(zip(files1, files2)):
		# ----------------------------------------------------------
		# load files
		data1 = fits.getdata(filename1)
		hdr1 = fits.getheader(filename1)
		data2 = fits.getdata(filename2)
		hdr2 = fits.getheader(filename2)
		# ----------------------------------------------------------
		# add keys
		sdict['target'].append(hdr1[KW_OBJNAME])
		sdict['file1'].append(os.path.basename(filename1))
		sdict['file2'].append(os.path.basename(filename2))
		sdict['spt'].append(hdr1[KW_SPT])
		sdict['date1'].append(hdr1[KW_DATE])
		sdict['date2'].append(hdr2[KW_DATE])
		sdict['mode'].append(hdr1[KW_MODE])
		sdict['Texp1'].append(hdr1[KW_EXPTIME])
		sdict['Texp2'].append(hdr2[KW_EXPTIME])
		sdict['AM1'].append(hdr1[KW_AIRMASS])
		sdict['seeing1'].append(hdr1[KW_SEEING])
		sdict['AM2'].append(hdr2[KW_AIRMASS])
		sdict['seeing2'].append(hdr2[KW_SEEING])
		# ----------------------------------------------------------
		# difference image
		diff = data1 - data2
		# ----------------------------------------------------------
		# loop around choosen orders
		for jt, order_num in enumerate(orders):
			# get rms
			rms = mp.robust_sigma(diff[order_num])
			# get median of data 1
			med1 = mp.nanmedian(data1)
			med2 = mp.nanmedian(data2)
			# get snr
			snr1 = med1 / rms
			snr2 = med2 / rms
			# add to storage
			sdict[f'measureSNR1[{echelle_orders[jt]}]'].append(snr1)
			sdict[f'measureSNR2[{echelle_orders[jt]}]'].append(snr1)
	# ----------------------------------------------------------
	# convert dict to table
	table = Table(sdict)
	# save to disk
	print(f'Saving to {OUTTABLE_CSV}')
	table.write(OUTTABLE_CSV, format='csv')
	print(f'Saving to {OUTTABLE_FITS}')
	table.write(OUTTABLE_FITS, format='fits')
	print(f'Saving to {OUTTABLE_TXT}')
	table.write(OUTTABLE_TXT, format='ascii.fixed_width')

# ==================================================================
# END OF CODE
# ==================================================================