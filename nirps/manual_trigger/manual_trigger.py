#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-03-08 at 11:59

@author: cook
"""
from typing import Any, Dict
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import warnings
import argparse
import yaml
import os


# =============================================================================
# Define variables
# =============================================================================
PROFILE_FILE = 'profiles.yaml'
# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def get_args():
    parser = argparse.ArgumentParser(description='Rene\'s magic trigger')
    # add obs dir
    parser.add_argument('--obsdir', type=str, default='None',
                        help='Observation directory name(s) separated by '
                             'commas')
    # test mode - do not run apero recipes just test
    parser.add_argument('--test', type=bool, default=False,
                        help='Do not run apero recipes just test')
    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return args


def get_settings():
    """
    Get the settings from the command line and yaml file

    :return:
    """
    # storage of parameters
    settings = dict()
    # get arguments
    args = get_args()
    # ----------------------------------------------------------------------
    # must have an observation directory
    if args.obsdir in ['None', None, 'Null', '']:
        raise ValueError('Must provide an observation directory')
    else:
        obs_dirs = args.obsdir.split(',')
    # ----------------------------------------------------------------------
    # add this to settings
    settings['OBS_DIRS'] = obs_dirs
    # add the test mode
    settings['TEST'] = args.test
    # ----------------------------------------------------------------------
    # read the yaml file and push into settings
    settings = read_yaml(PROFILE_FILE, settings)
    # ----------------------------------------------------------------------
    # return the settings
    return settings


def read_yaml(yaml_filename: str, settings: Dict[str, Any]) -> Dict[str, Any]:
    """
    Read the yaml file and add to settings

    :param yaml_filename: str, yaml file name
    :param settings: dict, settings dictionary

    :return: dict, updated settings dictionary
    """
    # read the yaml file
    with open(yaml_filename, 'r') as f:
        yaml_data = yaml.load(f, Loader=yaml.FullLoader)
    # add a profiles sub-dictionary
    settings['PROFILES'] = dict()
    # loop around yaml data
    for key, value in yaml_data.items():
        # if key is in settings
        settings['PROFILES'][key] = value
    # return settings
    return settings


def print_process(msg: str):
    """
    print a headering message
    :param msg: str, print message

    :return: None, prints a message
    """
    print('=' * 50)
    print(msg)
    print('=' * 50)


def update_apero_profile(profile: dict):
    """
    Update the apero profile with the correct paths
    :param profile: dict, the profile to update
    :return:
    """
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = profile['apero profile']
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # make sure parameters is reloaded (and not cached)
    return constants.load(cache=False)


def make_sym_links(settings: Dict[str, Any]):
    """
    Make the symbolic links for the raw data

    :param settings: dict, settings dictionary

    :return: None, creates symbolic links
    """
    # loop around profiles
    for profile in settings['PROFILES']:
        # print progress
        print(f'\tRunning profile: {profile}')
        # update the apero profile
        params = update_apero_profile(settings['PROFILES'][profile])
        # get raw directory path from profile
        inpath = settings['PROFILES'][profile]['raw dir']
        # get the raw directory from params
        outpath = params['DRS_DATA_RAW']
        # loop around obs dirs
        for obs_dir in settings['OBS_DIRS']:
            # get the full path
            full_inpath = os.path.join(inpath, obs_dir)
            full_outpath = os.path.join(outpath, obs_dir)
            # print sym link creation
            print(f'\t\tCreating symlink {full_outpath}')
            # deal with test mode
            if settings['TEST']:
                continue
            # check if the symlink exists
            if os.path.exists(full_outpath):
                # if it does then remove it
                os.remove(full_outpath)
            # make the symlink
            os.symlink(full_inpath, full_outpath)


def run_precheck(settings: Dict[str, Any]):
    """
    Run the precheck on all profiles

    :param settings: dict, settings dictionary
    """
    _  = settings
    return


def run_processing(settings: Dict[str, Any]):
    """
    Run the processing on all profiles

    :param settings: dict, settings dictionary
    """
    # loop around profiles
    for profile in settings['PROFILES']:
        # print progress
        print(f'\tRunning profile: {profile}')
        # update the apero profile
        _ = update_apero_profile(settings['PROFILES'][profile])
        # get the run file
        runfile = profile['apero run file']
        # need to import apero_processing
        from apero.tools.recipes.bin import apero_processing
        # run apero processing
        apero_processing.main(runfile=runfile,
                              include_obs_dirs=','.join(settings['OBS_DIRS']),
                              test=settings['TEST'])


def run_apero_get(settings: Dict[str, Any]):
    """
    Run the apero get on all profiles

    :param settings: dict, settings dictionary
    """
    # loop around profiles
    for profile in settings['PROFILES']:
        # print progress
        print(f'\tRunning profile: {profile}')
        # update the apero profile
        _ = update_apero_profile(settings['PROFILES'][profile])
        # get the output types
        outtypes = profile['apero out types']
        # get the dpr types
        dprtypes = profile['apero dpr types']
        # get the output path
        outpath = profile['lbl in path']
        # need to import apero_processing
        from apero.tools.recipes.bin import apero_get
        # run apero processing
        apero_get.main(objnames='*',
                       dprtypes=dprtypes, outtypes=outtypes,
                       outpath=outpath, fibers='A', symlinks=True,
                       test=settings['TEST'])


def run_apero_reduction_interface(settings: Dict[str, Any]):
    """
    Run the apero reduction interface

    :param settings: dict, settings dictionary
    """
    _  = settings
    return


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get settings
    trigger_settings = get_settings()
    # ----------------------------------------------------------------------
    # make symbolic links
    print_process('Making symbolic links')
    make_sym_links(trigger_settings)
    # ----------------------------------------------------------------------
    # run apero precheck on all profiles
    print_process('Running apero precheck')
    run_precheck(trigger_settings)
    # ----------------------------------------------------------------------
    # run apero processing on all profiles
    print_process('Running apero processing')
    run_processing(trigger_settings)
    # ----------------------------------------------------------------------
    # run apero get on all profiles
    print_process('Running apero get')
    run_apero_get(trigger_settings)
    # ----------------------------------------------------------------------
    # run the apero reduction interface
    print_process('Running apero reduction interface')
    run_apero_reduction_interface(trigger_settings)




# =============================================================================
# End of code
# =============================================================================
