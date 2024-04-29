#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-04-29 at 15:08

@author: cook
"""
import matplotlib.pyplot as plt

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from astropy import constants as cc
from astropy import units as uu

from apero import lang
from apero.base import base
from apero.core import constants
from apero.core.core import drs_database
from apero.core.core import drs_file
from apero.core.core import drs_log
from apero.core.utils import drs_recipe
from apero.core.utils import drs_startup
from apero.io import drs_image
from apero.science import velocity
from apero.science.calib import flat_blaze
from apero.science.calib import gen_calib
from apero.science.calib import wave
from apero.science.extract import other as extractother

from apero.science.velocity.gen_vel import get_coeff_dict

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'apero_wave_night_spirou.py'
__INSTRUMENT__ = 'SPIROU'
__PACKAGE__ = base.__PACKAGE__
__version__ = base.__version__
__author__ = base.__author__
__date__ = base.__date__
__release__ = base.__release__
# Get Logging function
WLOG = drs_log.wlog
# Get Recipe class
DrsRecipe = drs_recipe.DrsRecipe
# Get parameter class
ParamDict = constants.ParamDict
# Get the text types
textentry = lang.textentry
# define extraction code to use
EXTRACT_NAME = 'apero_extract_spirou.py'
# Speed of light
# noinspection PyUnresolvedReferences
speed_of_light_ms = cc.c.to(uu.m / uu.s).value

# ----------------------------------------------------------------------
# Main Code
# ----------------------------------------------------------------------
obs_dir = '2020-05-14'
hcfiles = ['2488251c_pp.fits', '2488410c_pp.fits']
fpfiles = ['2488244a_pp.fits', '2488406a_pp.fits']
kwargs = dict()

# assign function calls (must add positional)
fkwargs = dict(obs_dir=obs_dir, hcfiles=hcfiles, fpfiles=fpfiles,
               **kwargs)
# ----------------------------------------------------------------------
# deal with command line inputs / function call inputs
recipe, params = drs_startup.setup(__NAME__, __INSTRUMENT__, fkwargs)



mainname = __NAME__ + '._main()'
# get files
hcfiles = params['INPUTS']['HCFILES'][1]
# check qc
hcfiles = drs_file.check_input_qc(params, hcfiles, 'hc files')
# deal with (optional fp files)
if len(params['INPUTS']['FPFILES']) == 0:
    fpfiles = None
else:
    fpfiles = params['INPUTS']['FPFILES'][1]
    # check qc
    fpfiles = drs_file.check_input_qc(params, fpfiles, 'fp files')
    # must check fp files pass quality control
    fpfiles = gen_calib.check_fp_files(params, fpfiles)
# get list of filenames (for output)
rawhcfiles, rawfpfiles = [], []
for infile in hcfiles:
    rawhcfiles.append(infile.basename)
# deal with (optional fp files)
if fpfiles is not None:
    for infile in fpfiles:
        rawfpfiles.append(infile.basename)

    # deal with input data from function
    if 'hcfiles' in params['DATA_DICT']:
        hcfiles = params['DATA_DICT']['hcfiles']
        fpfiles = params['DATA_DICT']['fpfiles']
        rawhcfiles = params['DATA_DICT']['rawhcfiles']
        rawfpfiles = params['DATA_DICT']['rawfpfiles']
        combine = params['DATA_DICT']['combine']
    # combine input images if required
    elif params['INPUT_COMBINE_IMAGES']:
        # get combined file
        cond1 = drs_file.combine(params, recipe, hcfiles, math='mean')
        hcfiles = [cond1[0]]
        # get combined file
        if fpfiles is not None:
            cond2 = drs_file.combine(params, recipe, fpfiles, math='mean')
            fpfiles = [cond2[0]]
        combine = True
    else:
        combine = False
    # get the number of infiles
    num_files = len(hcfiles)

    # warn user if lengths differ
    if fpfiles is not None:
        if len(hcfiles) != len(fpfiles):
            wargs = [len(hcfiles), len(fpfiles)]
            WLOG(params, 'error', textentry('10-017-00002', args=wargs))
    # get the number of files
    num_files = len(hcfiles)
    # get the fiber types from a list parameter (or from inputs)
    fiber_types = drs_image.get_fiber_types(params)
    # get wave reference file (controller fiber)
    ref_fiber = params['WAVE_REF_FIBER']
    # load the calibration database
    calibdbm = drs_database.CalibrationDatabase(params)
    calibdbm.load_db()

    # ----------------------------------------------------------------------
    # Loop around input files
    # ----------------------------------------------------------------------
    for it in range(num_files):
        # ------------------------------------------------------------------
        # add level to recipe log
        log1 = recipe.log.add_level(params, 'num', it)
        # ------------------------------------------------------------------
        # set up plotting (no plotting before this)
        recipe.plot.set_location(it)
        # print file iteration progress
        drs_startup.file_processing_update(params, it, num_files)
        # get this iterations files
        hcfile = hcfiles[it]
        if fpfiles is None:
            fpfile = None
        else:
            fpfile = fpfiles[it]

        # -----------------------------------------------------------------
        # load initial wavelength solution (start point) for this fiber
        #    this should only be a reference wavelength solution
        wprops = wave.get_wavesolution(params, recipe, infile=hcfile,
                                        fiber=ref_fiber,
                                        database=calibdbm)

        # set up parameters
        eargs = [params, recipe, EXTRACT_NAME, hcfile, fpfile]
        ekwargs = dict(wavefile=wprops['WAVEFILE'], logger=log1)
        # run extraction
        hc_outputs, fp_outputs = extractother.extract_wave_files(*eargs,
                                                                 **ekwargs)

        # store rvs from ccfs
        rvs_all = dict()
        wprops_all = dict()
        # loop around fibers
        for fiber in fiber_types:
            # get fp e2ds file
            fp_e2ds_file = fp_outputs[fiber]
            # get the wave solution for this fiber
            wprops = wave.get_wavesolution(params, recipe, infile=hcfile,
                                           fiber=fiber,
                                           database=calibdbm)
            # -----------------------------------------------------------------
            # load the blaze file for this fiber
            bout = flat_blaze.get_blaze(params, fp_e2ds_file.header, ref_fiber)
            blaze_file, blaze_time, blaze = bout
            # compute the ccf
            ccfargs = [fp_e2ds_file, fp_e2ds_file.get_data(), blaze,
                       wprops['WAVEMAP'], fiber]
            rvprops = velocity.compute_ccf_fp(params, recipe, *ccfargs)
            # update ccf properties and push into wprops for wave sol outputs
            wprops, rvprops = wave.update_w_rv_props(wprops, rvprops, mainname)
            # -----------------------------------------------------------------
            # add to rv storage
            rvs_all[fiber] = rvprops
            wprops_all[fiber] = wprops

            # add to rv props
            rvs_all[fiber].set('WM_DV', value=np.nan, source=__NAME__)


    # =================================================================
    # Quality control
    # =================================================================
    qc_params = wave.wave_quality_control(params, wprops_all,
                                          rvs_all)



    pab = rvs_all['AB']
    pa = rvs_all['A']
    pb = rvs_all['B']
    pc = rvs_all['C']


    p1 = pab
    p2 = pc

    coeffs1 = get_coeff_dict(p1['CCF_FIT_COEFFS'], p1['CCF_FIT_NAMES'])
    coeffs2 = get_coeff_dict(p2['CCF_FIT_COEFFS'], p2['CCF_FIT_NAMES'])

    diffRV = coeffs1['RV'] - coeffs2['RV']


    differences = dict()
    for name in pab['CCF_FIT_NAMES']:
        differences[name] = coeffs1[name] - coeffs2[name]


    plt.close()
    fig, frames = plt.subplots(ncols=1, nrows=len(differences), sharex='all')

    for it, name in enumerate(list(differences.keys())):
        frames[it].plot(differences[name], marker='o', ls='None')
        frames[it].set_ylabel(f'{name} Difference')

    frames[-1].set_xlabel('Order Number')


    plt.plot(p1['RV_CCF'], p1['CCF_STACK'], 'bx')
    plt.plot(p1['RV_CCF'], p1['CCF_FIT_STACK'], 'b-')

    plt.plot(p2['RV_CCF'], p2['CCF_STACK'], 'rx')
    plt.plot(p2['RV_CCF'], p2['CCF_FIT_STACK'], 'r-')
    plt.show()


# =============================================================================
# End of code
# =============================================================================
