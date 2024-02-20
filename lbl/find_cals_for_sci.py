#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-02-20 at 12:03

@author: cook
"""
import glob
from tqdm import tqdm
import os
import shutil

# =============================================================================
# Define variables
# =============================================================================
SCI_FILELIST = ['2504190o_pp_e2dsff_tcorr_AB.fits',
                '2502850o_pp_e2dsff_tcorr_AB.fits',
                '2509620o_pp_e2dsff_tcorr_AB.fits',
                '2505260o_pp_e2dsff_tcorr_AB.fits',
                '2513090o_pp_e2dsff_tcorr_AB.fits',
                '2512280o_pp_e2dsff_tcorr_AB.fits',
                '2513610o_pp_e2dsff_tcorr_AB.fits',
                '2516560o_pp_e2dsff_tcorr_AB.fits',
                '2516070o_pp_e2dsff_tcorr_AB.fits',
                '2517810o_pp_e2dsff_tcorr_AB.fits',
                '2516990o_pp_e2dsff_tcorr_AB.fits',
                '2518350o_pp_e2dsff_tcorr_AB.fits',
                '2540640o_pp_e2dsff_tcorr_AB.fits',
                '2518790o_pp_e2dsff_tcorr_AB.fits',
                '2585920o_pp_e2dsff_tcorr_AB.fits',
                '2732630o_pp_e2dsff_tcorr_AB.fits',
                '2758790o_pp_e2dsff_tcorr_AB.fits',
                '2747240o_pp_e2dsff_tcorr_AB.fits',
                '2762570o_pp_e2dsff_tcorr_AB.fits',
                '2762030o_pp_e2dsff_tcorr_AB.fits',
                '2763280o_pp_e2dsff_tcorr_AB.fits',
                '2763100o_pp_e2dsff_tcorr_AB.fits']

SEARCH_PATH = '/cosmos99/spirou/apero-data/spirou_offline/red/'

COPY_PATH = '/cosmos99/spirou/misc/lbl-test/'

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def function1():
    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    cal_files = []
    wave_files = []

    # remove duplicates from sci_filelist
    sci_filelist = list(set(SCI_FILELIST))
    # print progress
    print('\n\nFinding cal files for {0} sci files'.format(len(sci_filelist)))
    # loop around sci files
    for sci_file in tqdm(sci_filelist):
        found_files = glob.glob(os.path.join(SEARCH_PATH, '*', sci_file))
        # assume the first one is the correct one
        found_path = os.path.dirname(found_files[0])

        # get the corresponding cal e2ds file
        cal_file = sci_file.replace('pp_e2dsff_tcorr_AB.fits',
                                    'pp_e2dsff_C.fits')
        cal_files.append(os.path.join(found_path, cal_file))

        # get the corresponding wave file in the obsdir
        wave_files += list(glob.glob(os.path.join(found_path,
                                                  '*wave_night*AB.fits')))
        wave_files += list(glob.glob(os.path.join(found_path,
                                                  '*wave_night*C.fits')))
    # make sure we don't have duplicates
    cal_files = list(set(cal_files))
    wave_files = list(set(wave_files))

    # print progress
    print('\n\nCopying cal files to {0}'.format(COPY_PATH))
    # copy files
    for cal_file in tqdm(cal_files):
        shutil.copy(cal_file, os.path.join(COPY_PATH, 'science', 'FP'))
    # print progress
    print('\n\nCopying wave files to {0}'.format(COPY_PATH))
    # copy files
    for wave_file in tqdm(wave_files):
        shutil.copy(wave_file, os.path.join(COPY_PATH, 'calib'))

# =============================================================================
# End of code
# =============================================================================
