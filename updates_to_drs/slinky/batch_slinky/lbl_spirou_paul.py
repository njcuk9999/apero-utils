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
def wrapper_paul():

    dir = '/space/spirou/LBL-PCA/SPIROU/science'

    objs1 = np.array(glob.glob(f'{dir}/*TELLU05'))
    objs2 = np.array(glob.glob(f'{dir}/*SLINKY05'))

    objs = np.append(objs1,objs2)
    
    objs = np.array([ob.split('/')[-1] for ob in objs])

    nn = np.zeros_like(objs,dtype = int)
    for i in range(len(nn)):
        nn[i] = len(glob.glob(dir+'/'+objs[i]+'/*'))

    valid = nn>10
    objs_to_do = list(objs[valid])



    #params = get_params(yaml)
    #objects_science = np.array(params['object_of_interest'])


    now = datetime.datetime.now()
    now_str = (now.isoformat().replace('-','').replace(':','').replace('T','_')).split('.')[0]

    # set up parameters
    rparams = dict()

    print(rparams)

    if True:#params['instrument'].upper() == 'SPIROU':
        rparams['INSTRUMENT'] = 'SPIROU'
        rparams['DATA_DIR'] = '/space/spirou/LBL-PCA/SPIROU'
        inst_folder = 'spirou'
        rsync_machine = 'spirou-client@maestria'
        inst_short = 'spirou'

    #if 'NIRPS' in params['instrument'].upper():
    #    rparams['INSTRUMENT'] = 'NIRPS_HE'
    #    rparams['DATA_DIR'] = '/space/spirou/LBL-PCA/NIRPS_HE'
    #    inst_folder = 'nirps_he'
    #    rsync_machine = 'nirps-client@maestria'
    #    inst_short = 'nirps'


    datadir0 = str(rparams['DATA_DIR'])

    rparams['DATA_DIR'] =  datadir0+'_'+now_str
    os.system('mkdir '+rparams['DATA_DIR'])

    ln_paths = ['calib','models','science']

    for path in ln_paths:
        cmd = 'ln -s '+datadir0+'/'+path+' '+rparams['DATA_DIR']+'/'
        os.system(cmd)
        print(cmd)

    rparams['DATA_SOURCE'] = 'CADC'

    #objs_to_do = []
    #objs_original_name = []
    #batchnames = []

    print(rparams['DATA_DIR']+'/science/')
    #print(objects_science)

    path_to_summary = '/space/spirou/LBL-PCA/wraps/batch_summary/'+rparams['INSTRUMENT']
    path_to_summary = path_to_summary+'/'+now_str
    sleep(3)
    


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
    ## Advanced settings
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

    path1 = rparams['DATA_DIR']

    tmp_dir = path1+'/lblstats_tmp'

    if not os.path.exists(tmp_dir):
        os.system('mkdir '+tmp_dir)
    else:
        os.system('rm -r '+tmp_dir)
        os.system('mkdir '+tmp_dir)

    for iobj in range(len(objs_to_do)):
        obj = objs_to_do[iobj]
        batchname = batchnames[iobj]
        obj_original = objs_original_name[iobj]

        print(obj, batchname, obj_original)

        if batchname == 'vanilla':
            continue

        name1 =f'{path1}/lblrdb/lbl_{obj}_{obj}.rdb'
        path2 = f'{tmp_dir}/{obj_original}'
        if not os.path.exists(path2):
            os.system('mkdir '+path2)

        name2 =f'{path2}/lbl_{obj_original}_{obj_original}_{batchname}.rdb'


        cmd = f'cp {name1} {name2}'
        print(cmd)
        os.system(cmd)

    path3 = f'/cosmos99/{inst_short}/apero-data/{inst_folder}_online/lbl/lblstats/'
    cmd_rsync = f'rsync -av {tmp_dir}/* {rsync_machine}:{path3}'
    print(cmd_rsync)
    os.system(cmd_rsync)
    os.system('rm -r ' + tmp_dir)


if __name__ == "__main__":
    #if len(sys.argv) ==2:
    #yaml_file = sys.argv[1]
    wrapper_paul()
    #else:
    #    print('You must pass a yaml file name')


# =============================================================================
# End of code
# =============================================================================
