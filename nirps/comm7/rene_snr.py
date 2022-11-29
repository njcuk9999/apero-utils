#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Code for measuring SNR using Rene's magic SNR function

Created on 2022-11-20 at 19:19

@author: cook
"""
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.special import erf


# =============================================================================
# Define variables
# =============================================================================
# set up path
PATH = '/nirps_raw/nirps/apero-data/drs-data/nirps_he_202211/red/'
FILENAME = '2022-11-28/NIRPS_2022-11-29T00_18_16_448_pp_e2dsff_A.fits'
# define number of pixels away from the center to use
WIDTH = 100
# debug plot
DEBUG = True
# list of debug plot orders (echelle order number)
DEBUG_ECHELLE_ORDERS = [90]


# =============================================================================
# Define functions
# =============================================================================
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

	return snr


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
	# get path
	abspath = os.path.join(PATH, FILENAME)
	# load spectrum
	data = fits.getdata(abspath)
	hdr = fits.getheader(abspath)
	# ----------------------------------------------------------------------
	# print object name
	print('=' * 50 + f'\n{hdr["OBJECT"]}\n' + '=' * 50)
	# ----------------------------------------------------------------------
	# get number of orders and pixel size
	nbo, nbxpix = data.shape
	# storage
	snr = dict()
	# loop around each order
	for order_num in range(nbo):
		# ----------------------------------------------------------------------
		# get echelle order
		echelle_order = hdr[f'WAVEEC{order_num:02d}']
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
			snr[order_num] = measure_snr(measure_data, order_num, echelle_order)
		# ----------------------------------------------------------------------
		# print value
		print(f'EXT[{order_num}]  ECHELLE[{echelle_order}]: {snr[order_num]}')


# =============================================================================
# End of code
# =============================================================================
