#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-02-02 at 11:16

@author: cook
"""
import os

# =============================================================================
# Define variables
# =============================================================================
modes = ['HE', 'HA']
TEST = False

aparams = dict()
aparams['HE'] = dict()
aparams['HA'] = dict()

aparams['HE']['inpath'] = '/nirps_raw/nirps/raw-data/helios_he'
aparams['HE']['outpath1'] = '/nirps_raw/nirps/apero-data/common/rawsym202211-HE'
aparams['HE']['outpath2'] = '/nirps_raw/nirps/apero-data/common/rawsym202211-HE-pernight'
aparams['HA']['inpath'] = '/nirps_raw/nirps/raw-data/helios_ha'
aparams['HA']['outpath1'] = '/nirps_raw/nirps/apero-data/common/rawsym202211-HA'
aparams['HA']['outpath2'] = '/nirps_raw/nirps/apero-data/common/rawsym202211-HA-pernight'


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # loop around modes
    for mode in modes:
        # get params
        params = aparams[mode]
        directories = []
        # find directories
        for root, dirs, files in os.walk(params['inpath']):
            for _dir in dirs:
                directory = os.path.join(root, _dir)
                if os.path.isdir(directory):
                    directories.append(directory)
        # make sym links
        for directory in directories:
            # set base directory name
            basedir = os.path.basename(directory)
            # make outpath
            outpath1 = os.path.join(params['outpath1'], basedir + '_helios')
            outpath2 = os.path.join(params['outpath2'], basedir + '_helios')
            # print and make outpath1
            print(f'{directory}-->{outpath1}')
            if not TEST:
                os.symlink(directory, outpath1, target_is_directory=True)
            # print and make outpath2
            print(f'{directory}-->{outpath2}')
            if not TEST:
                os.symlink(directory, outpath2, target_is_directory=True)

# =============================================================================
# End of code
# =============================================================================
