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
from pixpca import get_params
import numpy as np
import os
import glob
#yaml = 'yamls/params_tellu05.yaml'
from time import sleep
import datetime

# =============================================================================
# Start of code
# =============================================================================
def wrapper_slinky(params):
    objects_science = np.array(params['object_of_interest'])

    # set up parameters
    rparams = dict()


    datadir0 = str(params['output_slinky'])
    if not os.path.exists(datadir0):
        print(f"Data directory {datadir0} does not exist. Please check the path.")
        sys.exit(1)
    rparams['DATA_DIR'] = os.path.join(params['DATA_DIR'],params['now_str'])
    rparams['INSTRUMENT'] = params['instrument']
    
    if not os.path.exists(rparams['DATA_DIR']):
        os.system('mkdir '+rparams['DATA_DIR'])

    for path in ['science'] :
        datadir = os.path.join(rparams['DATA_DIR'], path)
        if not os.path.exists(datadir):
            print('Creating '+datadir)
            os.mkdir(datadir)
        else:
            print('Directory '+datadir+' already exists')
        paths = glob.glob(os.path.join(datadir0, '*'))
        for p in paths:
            if not os.path.exists(os.path.join(datadir, os.path.basename(p))):
                print('Linking '+p+' to '+datadir)
                cmd = 'ln -s '+p+' '+datadir
                os.system(cmd)
                print(cmd)
            else:
                print('File '+os.path.basename(p)+' already exists in '+datadir)
    
    if not os.path.exists(os.path.join(rparams['DATA_DIR'], 'templates')):
        print('Creating templates directory')
        os.mkdir(os.path.join(rparams['DATA_DIR'], 'templates'))
    if not os.path.exists(os.path.join(rparams['DATA_DIR'], 'calib')):
        print('Creating calib directory')
        os.mkdir(os.path.join(rparams['DATA_DIR'], 'calib'))

    rparams['DATA_SOURCE'] = 'CADC'

    # create duplicates with '_slinky' and '_slinky_[batchname]'. In the end there are 3x more objects

    suffixes = ['', '_slinky', '_slinky_'+params['batchname']]
    objs_to_do = []
    for obj in objects_science:
        for suffix in suffixes:
            if suffix == '':
                objs_to_do.append(obj)
            else:
                obj_new = obj+suffix
                if obj_new not in objs_to_do:
                    objs_to_do.append(obj_new)

    print(os.path.join(rparams['DATA_DIR'],'science'))
    print(objects_science)


    for ite in range(2):
        inst = ['SPIROU','NIRPS_HE'][ite]
        print(f'we clean {inst} files')
        f = glob.glob(f'/space/spirou/LBL-PCA/{inst}/science/*/*fits')
        for ff in f:
            try:
                _ = os.stat(ff)
            except:
                os.system('rm '+ff)
                print('err : '+ff)

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
    rparams['OBJECT_TEMPLATE'] = objs_to_do
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
    rparams['SKIP_LBL_TEMPLATE'] = False
    rparams['SKIP_LBL_MASK'] = False
    rparams['SKIP_LBL_COMPUTE'] = False
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

    # TODO uncomment this line
    lbl_wrap(rparams)


if __name__ == "__main__":
    if len(sys.argv) ==2:
        yaml_file = sys.argv[1]
        params = get_params(yaml_file)
        wrapper_slinky(params)
    else:
        print('You must pass a yaml file name')


# =============================================================================
# End of code
# =============================================================================
