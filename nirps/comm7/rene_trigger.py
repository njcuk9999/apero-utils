#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Code for measuring SNR using Rene's magic SNR function

Created on 2022-11-20 at 19:19

@author: cook
"""
import argparse
import glob
import os
import sys
from typing import List
import subprocess

from astropy.io import fits

# =============================================================================
# Define variables
# =============================================================================
MODES = ['HA', 'HE']
# set parameters for HA
aparams = dict()
aparams['HA'] = dict()
aparams['HA']['RAWPATH'] = '/nirps_raw/nirps/raw-data/nirps_ha/'
aparams['HA']['SYMPATH'] = '/nirps_raw/nirps/apero-data/common/rawsym202211-HA/'
aparams['HA']['APERO_PROFILE'] = 'setup_07254_nirps_ha'
aparams['HA']['REDPATH'] = '/nirps_raw/nirps/apero-data/drs-data/nirps_ha_202211/red/'
# set parameters for HE
aparams['HE'] = dict()
aparams['HE']['RAWPATH'] = '/nirps_raw/nirps/raw-data/nirps_he/'
aparams['HE']['SYMPATH'] = '/nirps_raw/nirps/apero-data/common/rawsym202211-HE/'
aparams['HE']['APERO_PROFILE'] = 'setup_07254_nirps_ha'
aparams['HE']['REDPATH'] = '/nirps_raw/nirps/apero-data/drs-data/nirps_he_202211/red/'

# dpr.type header keyword
KW_DPR_TYPE = 'HIERARCH ESO DPR TYPE'
# APERO environment
APERO_ENV = 'apero-env-07261'
# default run.ini file
DEFAULT_RUN_INI = 'comm_run.ini'
# test sys.argv
sys.argv = ['', '--mode=HA', '--dprfilter=OBJECT', '--obsdir=2022-11-28']


# ==================================================================
# CLASSES
# ==================================================================
class Color:
	red = '\033[1;91;1m'
	end = '\033[0;0m'


class NIRPSException(Exception):
	def __init__(self, message):
		# Call the base class constructor with the parameters it needs
		super().__init__(Color.red + message + Color.end)



# =============================================================================
# Define functions
# =============================================================================
def get_args():
	parser = argparse.ArgumentParser(description='Rene\'s magic trigger')
	# add mode
	parser.add_argument('--mode', type=str, choices=MODES,
	                    default='None', help='Mode to run trigger for.')
	# add dpr filter
	parser.add_argument('--dprfilter', type=str, default='None',
	                    help='Filter by string in dpr.type (i.e. OBJECT).'
	                         ' For multiple filters use semi-colon ";"')
	# add obs dir
	parser.add_argument('--obsdir', type=str, default='None',
	                    help='Observation directory name (if known)')
	# add obs dir
	parser.add_argument('--profile', type=str, default='None',
	                    help='APERO profile setup command to use')
	# add obs dir
	parser.add_argument('--runfile', type=str, default='None',
	                    help='APERO processing run.ini file to use'
	                         f'(defaults to {DEFAULT_RUN_INI})')
	# load arguments with parser
	args = parser.parse_args()
	# return arguments
	return args


def get_obs_dir(params: dict, args) -> str:
	# get raw path
	rawpath = params['RAWPATH']
	# deal with obs_dir set by user
	if args.obsdir != 'None':
		# get all observation directories
		obs_dirs = get_dirs(rawpath)
		for _obs_dir in obs_dirs:
			if args.obsdir == os.path.basename(_obs_dir):
				return args.obsdir
	# -------------------------------------------------------------------------
	# get all observation directories
	obs_dirs = get_dirs(rawpath)
	# get an observation directory
	user_obs_dir = None
	while user_obs_dir is None:
		print(f'Directories in : {rawpath}')
		for _obs_dir in obs_dirs:
			print(f'\t{_obs_dir}')
		uinput = input('\n\tChoose directory: ')

		if uinput in obs_dirs:
			user_obs_dir = uinput
		else:
			uinput = uinput.replace('\n', '')
			print(f'\t"{uinput}" is invalid. Please try again (Ctrl+C to exit)')
	return user_obs_dir


def get_dirs(path):
	_all = glob.glob(os.path.join(path, '*'))
	directories = []
	for _thing in _all:
		if os.path.isdir(_thing):
			directories.append(os.path.basename(_thing))
	return directories


def find_non_matching_files(rawfiles: List[str],
                            symfiles: List[str]) -> List[str]:
	# get basenames
	sym_basenames = [os.path.basename(filename) for filename in symfiles]
	raw_basenames = [os.path.basename(filename) for filename in rawfiles]
	# storage for non matches
	no_match_files = []
	# loop around files and find non-matches
	for it, filename in enumerate(raw_basenames):
		if filename not in sym_basenames:
			no_match_files.append(rawfiles[it])
	# return non matches
	return no_match_files


def filter_files(args, files: List[str]) -> List[str]:
	# deal with no filter
	if args.dprfilter == 'None':
		return files
	# filter files
	print('\nFiltering files. Please wait...')
	# get dpr filters
	dprfilters = args.dprfilter.split(';')
	# storage for valid files
	valid_files = []
	# loop around files and populate valid_files
	for filename in files:
		# get header
		hdr = fits.getheader(filename)
		# loop around filters
		for dprfilter in dprfilters:
			# if any key is in dpr.type we have a valid file
			if dprfilter in hdr[KW_DPR_TYPE]:
				valid_files.append(filename)
				break
	# return valid files
	return valid_files


def sym_link_files(params: dict, obsdir: str, files: List[str]):
	# get sym-link path
	sympath = params['SYMPATH']
	print('\nLinking files...\n')
	# loop around files
	for filename in files:
		# get base name
		basename = os.path.basename(filename)
		# get new path
		newpath = os.path.join(sympath, obsdir, basename)
		# make symlink
		print(f'Creating link: {newpath}')
		os.symlink(filename, newpath)
	print('\n\n')


def run_apero_processing(params, args, obsdir):
	# apero commands
	commands = ['conda deactivate'] * 5
	commands += [f'conda activate {APERO_ENV}']
	# -------------------------------------------------------------------------
	# start apero
	if args.profile == 'None':
		commands += [f'{params["APERO_PROFILE"]}']
	else:
		commands += [f'{args.profile}']
	# -------------------------------------------------------------------------
	# run apero
	if args.runfile == 'None':
		commands += [f'apero_processing.py {DEFAULT_RUN_INI}']
	else:
		commands += [f'apero_processing.py {args.runfile}']
	# -------------------------------------------------------------------------
	# get reduced dir
	reddir = os.path.join(params['REDPATH'], obsdir)
	# add the cd command
	commands += [f'cd {reddir}']
	# -------------------------------------------------------------------------
	# run commands
	print('\nPlease run the following:\n\n')
	for command in commands:
		print(f'{command}')


def main():
	# get arguments
	uargs = get_args()
	# set up path
	if uargs.mode not in MODES:
		raise NIRPSException(f'Mode: {uargs.mode} not valid')
	else:
		params = aparams[uargs.mode]
	# -------------------------------------------------------------------------
	# get raw path
	rawpath = params['RAWPATH']
	# get sym-link path
	sympath = params['SYMPATH']
	# -------------------------------------------------------------------------
	# get obs_dir
	obs_dir = get_obs_dir(params, uargs)
	# get files in raw obs_dir
	raw_dir = os.path.join(rawpath, obs_dir)
	all_raw_files = glob.glob(os.path.join(raw_dir, '*.fits'))
	# get files in raw sym obs_dur
	sym_dir = os.path.join(sympath, obs_dir)
	if os.path.exists(sym_dir):
		all_sym_files = glob.glob(os.path.join(sym_dir, '*.fits'))
	else:
		# make the directory
		os.mkdir(sym_dir)
		# set the list of rawsym files to None
		all_sym_files = None
	# -------------------------------------------------------------------------
	# create a list of basenames not in rawsym
	if all_sym_files is not None:
		files_to_add = find_non_matching_files(all_raw_files, all_sym_files)
	else:
		files_to_add = all_raw_files
	# -------------------------------------------------------------------------
	# need to filter files
	files_to_add = filter_files(uargs, files_to_add)
	# -------------------------------------------------------------------------
	# deal with no files to add
	if len(files_to_add) == 0:
		print(f'No files to add in {raw_dir}')
		return
	# -------------------------------------------------------------------------
	# sym link files
	sym_link_files(params, obs_dir, files_to_add)
	# -------------------------------------------------------------------------
	# run apero processing
	run_apero_processing(params, uargs, obs_dir)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
	ll = main()





# =============================================================================
# End of code
# =============================================================================
