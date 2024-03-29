#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Code for measuring SNR using Rene's magic SNR function

Created on 2022-11-20 at 19:19

@author: cook
"""
import os
import warnings
import sys

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.special import erf


# =============================================================================
# Define variables
# =============================================================================
# set up path
if len(sys.argv) == 2:
	PATH = sys.argv[1]
else:
	raise ValueError('Must provide a path as an argument')

FILE_SUFFIX = 'e2dsff_A.fits'
# define number of pixels away from the center to use
WIDTH = 100
# debug plot
DEBUG = False
# list of debug plot orders (echelle order number)
DEBUG_ECHELLE_ORDERS = [90]
# None for all otheriwse a list
PRINT_ECHELLE_ORDERS = [90]
# keywords
KW_OBJNAME = 'OBJECT'
KW_SPT = 'HIERARCH ESO OCS TARG SPTYPE'
KW_DATE = 'DATE'
KW_EXPTIME = 'EXPTIME'
KW_AIRMASS = 'HIERARCH ESO TEL AIRM START'
KW_SEEING = 'HIERARCH ESO TEL AMBI FWHM START'
KW_MODE = 'HIERARCH ESO INS MODE'
KW_DPRTYPE = 'HIERARCH ESO DPR TYPE'
KW_EXTSN059 = 'EXTSN059'
# output tables
OUTTABLE_CSV = '/nirps_raw/nirps/misc/comm7/rene_snr/all_results.csv'
OUTTABLE_FITS = '/nirps_raw/nirps/misc/comm7/rene_snr/all_results.fits'
OUTTABLE_TXT = '/nirps_raw/nirps/misc/comm7/rene_snr/all_results.txt'


# =============================================================================
# Define functions
# =============================================================================
def get_files(path):
	all_files = []
	# loop around all sub directories
	for root, dirs, files in os.walk(path):
		for filename in files:
			if FILE_SUFFIX is not None:
				if not filename.endswith(FILE_SUFFIX):
					continue
			all_files.append(os.path.join(root, filename))
	return all_files


def normal_fraction(sigma=1.0):
	"""
	Return the expected fraction of population inside a range
	(Assuming data is normally distributed)

	:param sigma: the number of sigma away from the median to be
	:return:
	"""
	# set function name
	# _ = display_func('normal_fraction', __NAME__)
	# return error function
	return erf(sigma / np.sqrt(2.0))


def estimate_sigma(tmp, sigma=1.0):
	"""
	Return a robust estimate of N sigma away from the mean

	:param tmp: np.array (1D) - the data to estimate N sigma of
	:param sigma: int, number of sigma away from mean (default is 1)

	:return: the sigma value
	"""
	# get formal definition of N sigma
	sig1 = normal_fraction(sigma)
	# get the 1 sigma as a percentile
	p1 = (1 - (1 - sig1) / 2) * 100
	# work out the lower and upper percentiles for 1 sigma
	upper = np.nanpercentile(tmp, p1)
	lower = np.nanpercentile(tmp, 100 - p1)
	# return the mean of these two bounds
	return (upper - lower) / 2.0


def measure_snr(vector, order_num, echelle_order):
	"""
	Measure the SNR in the order using Rene's magic fit
	(fit a quadratic of form ax^2 + bx + c)

	to points (x1, y1), (x2, y2), (x3, y3), (x4, y4)
	where x1 = -2, x2 = -1, x3 = 0, x4 = +1
	      y3 = c

	y1 = 4a - 2b + c
	y2 = a - b + c
	y4 = a + b + c

	Thus:
		c = y4/3 + y2 - y1/3

	:param vector:
	:param order_num:
	:param echelle_order:
	:return:
	"""
	# calculate all dy values
	dy_arr = []
	fit_arr = []
	real_arr = []
	# loop around valid pixels
	for it in range(2, len(vector) - 1):
		# get the four points for rene diff
		points = [vector[-2 + it], vector[-1 + it], vector[it], vector[it + 1]]
		# check for nans
		if np.sum(np.isfinite(points)) != 4:
			continue
		# rene fit
		yfit, dy = rene_fit(*points)
		# add value to dy array
		dy_arr.append(dy)
		fit_arr.append(yfit)
		real_arr.append(vector[it])
	# robust sigma of dy_arr
	rms = estimate_sigma(dy_arr)
	# measure noise (from Rene magic)
	#  rms^2 = (20/9) * noise^2
	noise = np.sqrt(9 / 20) * rms
	# measure median
	med = np.nanmedian(vector)
	# calculate snr
	snr = med / noise
	# debug plot
	if DEBUG and echelle_order in DEBUG_ECHELLE_ORDERS:
		fig, frames = plt.subplots(ncols=1, nrows=2)
		frames[0].plot(real_arr, color='g', marker='o', label='data')
		frames[0].plot(fit_arr, color='r', marker='x', label='fit')
		frames[0].legend(loc=0)
		frames[1].plot(dy_arr, color='k', marker='o', ls='None')
		plt.suptitle(f'EXT[{order_num}]  ECHELLE[{echelle_order}]: \n'
		             f'rms={rms:.4f}    noise={noise:.4f}    med={med:.4f}    snr={snr:.4f}')
		plt.show()

	# property dict
	pdict = dict()
	pdict['SNR'] = snr
	pdict['MED'] = med
	pdict['NOISE'] = noise
	pdict['RMS'] = rms

	return pdict


def rene_fit(y1, y2, y3, y4):
	"""
	Rene's magic fit function

	c = y4/3 + y2 - y1/3

	:param y1: a point 2 pixels before fit pixel
	:param y2: a point 1 pixel before fit pixel
	:param y3: the point to fit
	:param y4: the point 1 pixel after fit pixel

	:return: fit, residual
	"""
	fity3 = y4 / 3 + y2 - y1 / 3
	resdiual = fity3 - y3
	# return outputs
	return fity3, resdiual


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
	# ----------------------------------------------------------------------
	if os.path.isdir(PATH):
		individual = False
		files = get_files(PATH)
	else:
		individual = True
		files = [PATH]
	# store all data
	sdict = dict()
	sdict['target'] = []
	sdict['file'] = []
	sdict['spt'] = []
	sdict['date'] = []
	sdict['dpr.type'] = []
	sdict['mode'] = []
	sdict['Texp'] = []
	sdict['airmass'] = []
	sdict['seeing'] = []
	sdict['EXTSN059'] = []

	# loop around files
	for abspath in files:
		# load spectrum
		data = fits.getdata(abspath)
		hdr = fits.getheader(abspath)
		# ----------------------------------------------------------------------
		if not individual:
			if 'OBJECT' not in hdr[KW_DPRTYPE]:
				continue
		# ----------------------------------------------------------------------
		# print object name
		print('\n' + '=' * 50 + f'\n{hdr[KW_OBJNAME]}\n' + '=' * 50)
		print('\tDATE: ', hdr[KW_DATE])
		print('\tSPT: ', hdr[KW_SPT])
		print('\tEXPTIME: ', hdr[KW_EXPTIME])
		print('\tAIRMASS: ', hdr[KW_AIRMASS])
		print('\tSEEING: ', hdr[KW_SEEING])
		print('\tMODE: ', hdr[KW_MODE])
		print('\tDPRTYPE: ', hdr[KW_DPRTYPE])
		print('\tEXTSN059: ', hdr[KW_EXTSN059])
		print('=' * 50)
		# ----------------------------------------------------------------------
		# get number of orders and pixel size
		nbo, nbxpix = data.shape
		# storage
		snr, med = dict(), dict()


		# add keys
		sdict['target'].append(hdr[KW_OBJNAME])
		sdict['file'].append(os.path.basename(abspath))
		sdict['spt'].append(hdr[KW_SPT])
		sdict['date'].append(hdr[KW_DATE])
		sdict['dpr.type'].append(hdr[KW_DPRTYPE])
		sdict['mode'].append(hdr[KW_MODE])
		sdict['Texp'].append(hdr[KW_EXPTIME])
		sdict['airmass'].append(hdr[KW_AIRMASS])
		sdict['seeing'].append(hdr[KW_SEEING])
		sdict['EXTSN059'].append(hdr[KW_EXTSN059])

		# loop around each order
		for order_num in range(nbo):
			# ----------------------------------------------------------------------
			# get echelle order
			echelle_order = hdr[f'WAVEEC{order_num:02d}']
			# print value
			if PRINT_ECHELLE_ORDERS is None:
				print_order = True
			elif echelle_order in PRINT_ECHELLE_ORDERS:
				print_order = True
			else:
				print_order = False
			# skip all orders we don't want to print
			if not print_order:
				continue
			# ----------------------------------------------------------------------
			# get central point
			cent = nbxpix // 2
			lower = cent - WIDTH
			upper = cent + WIDTH
			# get cut down range
			measure_data = data[order_num][lower:upper]
			# ----------------------------------------------------------------------
			# measure SNR
			with warnings.catch_warnings(record=True) as _:
				sout = measure_snr(measure_data, order_num, echelle_order)
			# ----------------------------------------------------------------------
			print(f'EXT[{order_num}]  ECHELLE[{echelle_order}] '
			      f'SNR: {sout["SNR"]:.4f}   MED: {sout["MED"]:.4f}\n\n')
			# ----------------------------------------------------------------------
			# add keys
			snr_key = f'SNR[{echelle_order}]'
			med_key = f'MED[{echelle_order}]'
			noise_key = f'NOISE[{echelle_order}]'
			rms_key  = f'RMS[{echelle_order}]'
			# ----------------------------------------------------------------------
			if snr_key not in sdict:
				sdict[snr_key] = [sout['SNR']]
				sdict[med_key] = [sout['MED']]
				sdict[noise_key] = [sout['NOISE']]
				sdict[rms_key] = [sout['RMS']]
			else:
				sdict[snr_key] += [sout['SNR']]
				sdict[med_key] += [sout['MED']]
				sdict[noise_key] += [sout['NOISE']]
				sdict[rms_key] += [sout['RMS']]
	# ------------------------------------------------------------------------------
	# save to disk
	if len(files) > 0:

		if len(np.unique(sdict['mode'])) == 1:
			table_ext = sdict['mode'][0]
			# update filenames
			OUTTABLE_CSV = OUTTABLE_CSV.replace('.csv', table_ext + '.csv')
			OUTTABLE_FITS = OUTTABLE_FITS.replace('.fits', table_ext + '.fits')
			OUTTABLE_TXT = OUTTABLE_TXT.replace('.txt', table_ext + '.txt')
		# convert sdict to table
		table = Table(sdict)
		# save to disk
		print(f'Saving to {OUTTABLE_CSV}')
		table.write(OUTTABLE_CSV, format='csv', overwrite=True)
		print(f'Saving to {OUTTABLE_FITS}')
		table.write(OUTTABLE_FITS, format='fits', overwrite=True)
		print(f'Saving to {OUTTABLE_TXT}')
		table.write(OUTTABLE_TXT, format='ascii.fixed_width', overwrite=True)

# =============================================================================
# End of code
# =============================================================================
