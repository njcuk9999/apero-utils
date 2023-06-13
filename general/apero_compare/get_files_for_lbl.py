#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-01-10 at 13:36

@author: cook
"""
from astropy.io import fits
import glob
import os


# =============================================================================
# Define variables
# =============================================================================
name1 = 'spirou@rali'
name2 = 'cook@jupiter'
name3 = 'cook@nb19'
name4 = 'newworlds'
name5 = 'lam'
name6 = 'halau'

ref_name = 'spirou@rali'

inpaths = dict()
inpaths[name1] = '/scratch2/spirou/drs-data/minidata2_07284_20230608_rali/red/'
inpaths[name2] = '/scratch2/spirou/drs-data/minidata2_07XXX/reduced/'
inpaths[name3] = '/scratch2/spirou/drs-data/minidata2_07284_20230608_nb19/red/'
# inpaths[name4] = '/scratch2/spirou/drs-data/minidata2_07275_newworld/red/'
# inpaths[name5] = '/scratch2/spirou/drs-data/minidata2_07275_lam/red/'
# inpaths[name6] = '/scratch2/spirou/drs-data/minidata2_07275_halau/red/'

outpaths = dict()
outpaths[name1] = '/scratch3/lbl/data/minidata_comp/minidata2_07284_20230608_rali/science/'
outpaths[name2] = '/scratch3/lbl/data/minidata_comp/minidata2_07284_20230608_jupiter/science/'
outpaths[name3] = '/scratch3/lbl/data/minidata_comp/minidata2_07284_20230608_nb19/science/'
# outpaths[name4] = '/scratch3/lbl/data/minidata_comp/minidata2_07275_newworld/science/'
# outpaths[name5] = '/scratch3/lbl/data/minidata_comp/minidata2_07275_lam/science/'
# outpaths[name6] = '/scratch3/lbl/data/minidata_comp/minidata2_07275_halau/rscienceed/'

objname = 'GL699'


filestr = '*_pp_e2dsff_tcorr_AB.fits'
blazestr = '*_blaze_AB.fits'


# copy Gl699 files

# loop around names
for name in inpaths:
	#
	print('='*50)
	print(f' Processing {name}')
	print('='*50)
	# construct paths
	in_abspath = os.path.join(inpaths[name], '*', filestr)
	out_abspath = os.path.join(outpaths[name], objname)
	# make out_abspath if it doesn't exist
	if not os.path.exists(out_abspath):
		os.makedirs(out_abspath)
	# get all files
	files = glob.glob(in_abspath)
	# loop around files
	for filename in files:

		hdr = fits.getheader(filename)
		if hdr['DRSOBJN'] != 'GL699':
			continue

		out_file = os.path.join(out_abspath, os.path.basename(filename))

		if not os.path.exists(out_file):
			print(f'\tAdding {out_file}')
			os.symlink(filename, out_file)
		else:
			print(f'\tFound {out_file}')


# copy blaze files
for name in inpaths:

	print('='*50)
	print(f' Processing {name}')
	print('='*50)
	# construct paths
	in_abspath = os.path.join(inpaths[name], '*', blazestr)
	out_abspath = os.path.join(outpaths[name].replace('science', 'calib'))
	# make out_abspath if it doesn't exist
	if not os.path.exists(out_abspath):
		os.makedirs(out_abspath)
	# get all files
	files = glob.glob(in_abspath)
	# loop around files
	for filename in files:
		out_file = os.path.join(out_abspath, os.path.basename(filename))
		if not os.path.exists(out_file):
			print(f'\tAdding {out_file}')
			os.symlink(filename, out_file)
		else:
			print(f'\tFound {out_file}')