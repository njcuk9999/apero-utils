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
from apero.core import constants
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from tqdm import tqdm
import matplotlib.pyplot as plt

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
KW_DPRTYPE = 'HIERARCH ESO DPR TYPE'


# define orders to measure SNR in
orders = [29, 59]
echelle_orders = [120, 90]
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
	for it in tqdm(range(len(files[:-1]))):
		# get the second filename
		filename1 = files[it]
		filename2 = files[it + 1]
		# load header filename 1
		hdr1 = fits.getheader(filename1)
		objname1 = hdr1[KW_OBJNAME]
		# load header filename 2
		hdr2 = fits.getheader(filename2)
		objname2 = hdr2[KW_OBJNAME]
		# do not compare non objects
		if 'OBJECT' not in hdr1[KW_DPRTYPE] or 'OBJECT' not in hdr2[KW_DPRTYPE]:
			continue
		# do not compare different modes
		if hdr1[KW_MODE] != hdr2[KW_MODE]:
			continue
		# if exptime is very different don't use
		delta_exptime = abs(hdr1[KW_EXPTIME]/hdr2[KW_EXPTIME])
		if delta_exptime < 0.8 or delta_exptime > 1.2:
			continue

		# if the next file has the same object name we say they are consecutive
		if objname1 != objname2:
			continue
		# if the date is not within 2 exptimes
		delta_time = Time(hdr1[KW_DATE], format='fits') - Time(hdr2[KW_DATE], format='fits')
		if abs(delta_time) > 2 * hdr1[KW_EXPTIME]:
			continue
		# if we reached this point files are considered consecutive
		first_file.append(filename1)
		second_file.append(filename2)
	# return two list of files
	return first_file, second_file


# ==================================================================
# START OF CODE
# ==================================================================
if __name__ == '__main__':
	# get apero parameters
	aparams = constants.load()
	blaze_size = aparams['FF_BLAZE_HALF_WINDOW']
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
	sdict['target1'], sdict['target2'] = [], []
	sdict['file1'], sdict['file2'] = [], []
	sdict['spt1'], sdict['spt2'] = [], []
	sdict['date1'], sdict['date2'] = [], []
	sdict['mode1'], sdict['mode2'] = [], []
	sdict['Texp1'], sdict['Texp2'] = [], []
	sdict['AM1'], sdict['AM2'] = [], []
	sdict['seeing1'], sdict['seeing2'] = [], []
	for order_num in echelle_orders:
		sdict[f'measure[{order_num}]SNR1'] = []
		sdict[f'measure[{order_num}]SNR2'] = []
		sdict[f'hdr[{order_num}]SNR1'] = []
		sdict[f'hdr[{order_num}]SNR2'] = []
	# print process
	print('Calculating SNR')
	# work out the snr
	for it in tqdm(range(len(files1))):
		# get filenames
		filename1, filename2 = files1[it], files2[it]
		# ----------------------------------------------------------
		# load files
		data1 = fits.getdata(filename1)
		hdr1 = fits.getheader(filename1)
		data2 = fits.getdata(filename2)
		hdr2 = fits.getheader(filename2)
		# ----------------------------------------------------------
		# add keys
		sdict['target1'].append(hdr1[KW_OBJNAME])
		sdict['target2'].append(hdr2[KW_OBJNAME])
		sdict['file1'].append(os.path.basename(filename1))
		sdict['file2'].append(os.path.basename(filename2))
		sdict['spt1'].append(hdr1[KW_SPT])
		sdict['spt2'].append(hdr2[KW_SPT])
		sdict['date1'].append(hdr1[KW_DATE])
		sdict['date2'].append(hdr2[KW_DATE])
		sdict['mode1'].append(hdr1[KW_MODE])
		sdict['mode2'].append(hdr2[KW_MODE])
		sdict['Texp1'].append(hdr1[KW_EXPTIME])
		sdict['Texp2'].append(hdr2[KW_EXPTIME])
		sdict['AM1'].append(hdr1[KW_AIRMASS])
		sdict['AM2'].append(hdr2[KW_AIRMASS])
		sdict['seeing1'].append(hdr1[KW_SEEING])
		sdict['seeing2'].append(hdr2[KW_SEEING])
		# ----------------------------------------------------------
		# ----------------------------------------------------------
		# loop around choosen orders
		for jt, order_num in enumerate(orders):
			# rms = mp.estimate_sigma(diff[order_num])
			# get central pixel
			cent_pos = len(data1[order_num]) // 2
			blaze_lower = cent_pos - blaze_size
			blaze_upper = cent_pos + blaze_size
			# cut the data
			datacut1 = data1[order_num][blaze_lower:blaze_upper]
			datacut2 = data2[order_num][blaze_lower:blaze_upper]
			# get median of data 1
			med1 = mp.nanmedian(datacut1)
			med2 = mp.nanmedian(datacut2)
			# get rms
			rms = mp.estimate_sigma((datacut1/med1)/(datacut2/med2))
			# work out the measured noise
			noise = np.sqrt(2) / rms
			# get snr
			snr1 = med1 / noise
			snr2 = med2 / noise
			# get hdr snr
			hsnr1 = hdr1[f'EXTSN{order_num:03d}']
			hsnr2 = hdr2[f'EXTSN{order_num:03d}']

			# add to storage
			sdict[f'measure[{echelle_orders[jt]}]SNR1'].append(snr1)
			sdict[f'measure[{echelle_orders[jt]}]SNR2'].append(snr2)
			sdict[f'hdr[{echelle_orders[jt]}]SNR1'].append(hsnr1)
			sdict[f'hdr[{echelle_orders[jt]}]SNR2'].append(hsnr2)


		for jt, order_num in enumerate(orders):
			# rms = mp.estimate_sigma(diff[order_num])
			# get central pixel
			cent_pos = len(data1[order_num]) // 2
			blaze_lower = cent_pos - blaze_size
			blaze_upper = cent_pos + blaze_size
			# cut the data
			datacut1 = data1[order_num][blaze_lower:blaze_upper]
			datacut2 = data2[order_num][blaze_lower:blaze_upper]
			# get median of data 1
			med1 = mp.nanmedian(datacut1)
			med2 = mp.nanmedian(datacut2)
			# get snr
			hsnr1 = sdict[f'hdr[{echelle_orders[jt]}]SNR1'][-1]
			hsnr2 = sdict[f'hdr[{echelle_orders[jt]}]SNR2'][-1]
			snr1 = sdict[f'measure[{echelle_orders[jt]}]SNR1'][-1]
			snr2 = sdict[f'measure[{echelle_orders[jt]}]SNR2'][-1]

			if abs(hsnr1/snr1) < 0.8 or abs(hsnr1/snr1) > 1.2:

				keys, values, colors, store = [], [], [], []
				for kt, key in enumerate(sdict):

					if kt % 2 == 0:
						store = [sdict[key][-1]]
					else:
						store += [sdict[key][-1]]
						values.append(store)
						colors.append(['blue', 'red'])
						store = []
						keys.append(key[:-1])


				plt.close()
				fig, frames = plt.subplots(ncols=1, nrows=3)
				frames[0].plot(data1[order_num]/med1, color='b', alpha=0.7, label='obj1')
				frames[0].plot(data2[order_num]/med2, color='r', alpha=0.7, label='obj2')
				frames[0].set(title='data/med')
				frames[0].legend(loc=0)

				frames[1].plot((data1[order_num]/med1) / (data2[order_num]/med2), color='g')
				frames[1].set(title='ratio')

				tbl = frames[2].table(cellText=values, rowLabels=keys,
				                      colLabels=['obj1', 'obj2'], loc=0,
				                      cellColours=colors)
				for cell in tbl._cells:
					tbl._cells[cell].set_alpha(.3)
				tbl.auto_set_font_size(False)
				tbl.set_fontsize(8)
				tbl.scale(0.75, 0.75)
				frames[2].axis('off')
				plt.subplots_adjust(left=0.05, right=0.95, bottom=0.15, top=0.95)
				plt.show()
				plt.close()


	# ----------------------------------------------------------
	# convert dict to table
	table = Table(sdict)
	# save to disk
	print(f'Saving to {OUTTABLE_CSV}')
	table.write(OUTTABLE_CSV, format='csv', overwrite=True)
	print(f'Saving to {OUTTABLE_FITS}')
	table.write(OUTTABLE_FITS, format='fits', overwrite=True)
	print(f'Saving to {OUTTABLE_TXT}')
	table.write(OUTTABLE_TXT, format='ascii.fixed_width', overwrite=True)


# ==================================================================
# END OF CODE
# ==================================================================
