#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Code to get all filetypes for a list of objects for all APERO reductions
inside "workdir"

Created on 2021-08-31

@author: cook
"""
from astropy.io import fits
import glob
import os
import shutil
from tqdm import tqdm


# =============================================================================
# Define variables
# =============================================================================
# define where to get these files from
workdir = '/scratch2/drs-data/'
# workdir = '/spirou2/minidata20210916/'
# define where to copy these files to
output = '/spirou2/minidata/'
# reset
RESET = False
# define the targets (as in DRSOBJN)
# targets = ['GL699']
# targets = None --> all targets
targets = None

# define what type of files to get
outtypes = ['EXT_E2DS_FF', 'TELLU_OBJ', 'TELLU_TEMP', 'TELLU_TEMP_S1D']
fibers = ['AB', 'C']
dprtypes = ['OBJ_FP', 'OBJ_DARK', 'POLAR_FP', 'POLAR_DARK', 'FP_FP']


def get_fiber(header, filename):
    if 'FIBER' in header:
        return str(header['FIBER'])
    if 'AB.fits' in filename:
        return 'AB'
    elif 'A.fits' in filename:
        return 'A'
    elif 'B.fits' in filename:
        return 'B'
    elif 'C.fits' in filename:
        return 'C'
    else:
        return 'NULL'


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # get directories
    directories = glob.glob(workdir + 'minidata*')
    # loop around directories
    for directory in directories:
        # print progress
        print('Processing directory: {0}'.format(directory))
        # get outdir
        outdir = os.path.join(output, os.path.basename(directory))
        # rm outdir if it exists
        if os.path.exists(outdir) and RESET:
            os.system('rm -rfv {0}'.format(outdir))
            # make out dir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        # get pattern
        pattern = directory + '/reduced/'
        # get files
        all_files = glob.glob(pattern + '*/*.fits')
        # loop around files and add to outdir (if correct target)
        for filename in tqdm(all_files):
            # get header for this file
            header = fits.getheader(filename, ext=0)
            if 'DRSOBJN' in header:
                objname = str(header['DRSOBJN'])
            else:
                try:
                    header = fits.getheader(filename, ext=1)
                    objname = str(header['DRSOBJN'])
                # skip if we still don't have object name
                except Exception as _:
                    continue
            # get dprtype fiber and outtype
            dprtype = str(header.get('DPRTYPE', 'NULL'))
            fiber = get_fiber(header, filename)
            outtype = str(header.get('DRSOUTID', 'NULL'))
            # remove header
            del header
            # do not continue if objname not in targets
            if targets is not None:
                if objname not in targets:
                    continue
            if dprtypes is not None:
                if dprtype not in dprtypes:
                    continue
            if fibers is not None:
                if fiber not in fibers:
                    continue
            if outtypes is not None:
                if outtype not in outtypes:
                    continue
            # reate out path
            outpath = os.path.join(outdir, objname)
            # create object name dir
            if not os.path.exists(outpath):
                os.mkdir(outpath)
            # construct outpath
            outfilename = os.path.join(outpath, os.path.basename(filename))
            # copy
            if not os.path.exists(outfilename):
                shutil.copy(filename, outfilename)

# =============================================================================
# End of code
# =============================================================================
