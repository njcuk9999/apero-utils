#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wrap file for LBL recipes

Created on 2024-10-09 14:53:27.453

@author: Neil Cook, Etienne Artigau, Charles Cadieux, Thomas Vandal, Ryan Cloutier, Pierre Larue

user: spirou@maestria
lbl version: 0.64.0
lbl date: 2024-09-16
"""
import sys
from lbl import lbl_wrap
#from pixpca import get_params
#import numpy as np
#import os
#import glob
#yaml = 'yamls/params_tellu05.yaml'
#from time import sleep
#import datetime

# =============================================================================
# Start of code
# =============================================================================
def wrapper_slinky(nth_obj = 0):
    # set up parameters
    rparams = dict()

    rparams['INSTRUMENT'] = 'NIRPS_HE'
    rparams['DATA_DIR'] = '/space/spirou/LBL-PCA/NIRPS_HE'
    rparams['DATA_SOURCE'] = 'CADC'
    inst_folder = 'nirps_he'
    rsync_machine = 'nirps-client@maestria'
    inst_short = 'nirps'



    objs = 'GJ3090','GL406','PROXIMA','TOI756','TOI4552'

    obj = objs[nth_obj]
    objs_to_do = [obj,obj+'mask5.0sig',obj+'mask10.0sig']

    if obj == 'GJ3090':
        objs_templates = ['GJ3090','GJ3090','GJ3090']
    if obj == 'GL406':
        objs_templates = ['GL406','GL406','GL406']
    if obj == 'PROXIMA':
        objs_templates = ['PROXIMA','PROXIMA','PROXIMA']
    if obj == 'TOI756':
        objs_templates = ['GJ3090','GJ3090','GJ3090']
    if obj == 'TOI4552':
        objs_templates = ['PROXIMA','PROXIMA','PROXIMA']


    # The input file string (including wildcards) - if not set will use all
    #   files in the science directory (for this object name)
    # rparams['INPUT_FILE'] = '*'
    # The input science data are blaze corrected
    rparams['BLAZE_CORRECTED'] = False
    # Override the blaze filename
    #      (if not set will use the default for instrument)
    # rparams['BLAZE_FILE'] = 'blaze.fits'
    # -------------------------------------------------------------------------
    # science criteria
    # -------------------------------------------------------------------------
    # The data type (either SCIENCE or FP or LFC)
    rparams['DATA_TYPES'] = ['SCIENCE']*len(objs_to_do)
    # The object name (this is the directory name under the /science/
    #    sub-directory and thus does not have to be the name in the header
    rparams['OBJECT_SCIENCE'] = objs_to_do
    # This is the template that will be used or created (depending on what is
    #   run)
    rparams['OBJECT_TEMPLATE'] = objs_templates
    # This is the object temperature in K - used for getting a stellar model
    #   for the masks it only has to be good to a few 100 K
    rparams['OBJECT_TEFF'] = [3000]*len(objs_to_do)
    # -------------------------------------------------------------------------
    # what to run and skip if already on disk
    # -------------------------------------------------------------------------
    # Whether to run the telluric cleaning process (NOT recommended for data
    #   that has better telluric cleaning i.e. SPIROU using APERO)
    rparams['RUN_LBL_TELLUCLEAN'] = False
    # Whether to create templates from the data in the science directory
    #   If a template has been supplied from elsewhere this set is NOT required
    rparams['RUN_LBL_TEMPLATE'] = True
    # Whether to create a mask using the template created or supplied
    rparams['RUN_LBL_MASK'] = True
    # Whether to run the LBL compute step - which computes the line by line
    #   for each observation
    rparams['RUN_LBL_COMPUTE'] = True
    # Whether to run the LBL compile step - which compiles the rdb file and
    #   deals with outlier rejection
    rparams['RUN_LBL_COMPILE'] = True
    # whether to skip observations if a file is already on disk (useful when
    #   adding a few new files) there is one for each RUN_XXX step
    #   - Note cannot skip tellu clean
    rparams['SKIP_LBL_TEMPLATE'] = True
    rparams['SKIP_LBL_MASK'] = True
    rparams['SKIP_LBL_COMPUTE'] = True
    rparams['SKIP_LBL_COMPILE'] = False
    # -------------------------------------------------------------------------
    # LBL settings
    # -------------------------------------------------------------------------
    # You can change any setting in parameters (or override those changed
    #   by specific instruments) here
    # -------------------------------------------------------------------------
    # Advanced settings
    #   Do not use without contacting the LBL developers
    # -------------------------------------------------------------------------
    # Dictionary of table name for the file used in the projection against the
    #     derivative. Key is to output column name that will propagate into the
    #     final RDB table and the value is the filename of the table. The table
    #     must follow a number of characteristics explained on the LBL website.
    rparams['RESPROJ_TABLES'] = {
            'DTEMP3000': 'temperature_gradient_3000.fits',
            'DTEMP3500': 'temperature_gradient_3500.fits',
            'DTEMP4000': 'temperature_gradient_4000.fits',
            'DTEMP4500': 'temperature_gradient_4500.fits',
            'DTEMP5000': 'temperature_gradient_5000.fits',
            'DTEMP5500': 'temperature_gradient_5500.fits',
            'DTEMP6000': 'temperature_gradient_6000.fits'
            }

    # Rotational velocity parameters, should be a list of two values, one being
    #     the epsilon and the other one being the vsini in km/s as defined in the
    #     PyAstronomy.pyasl.rotBroad function
    # rparams['ROTBROAD'] = []

    # turn on plots
    rparams['PLOT'] = False

    # -------------------------------------------------------------------------
    # Run the wrapper code using the above settings
    # -------------------------------------------------------------------------
    # run main

    #   uncomment this line
    lbl_wrap(rparams)



#wrapper_slinky()

# =============================================================================
# End of code
# =============================================================================
