#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-03-08 at 11:59

@author: cook
"""
import argparse
import os
import shutil
import sys
from typing import Any, Dict

import astropy.units as uu
import numpy as np
import yaml
from astropy.time import Time

# =============================================================================
# Define variables
# =============================================================================
# start time
START_TIME = Time.now()
# -----------------------------------------------------------------------------


# =============================================================================
# Define functions
# =============================================================================
def get_args():
    parser = argparse.ArgumentParser(description='Rene\'s magic trigger')
    # add obs dir
    parser.add_argument('profile', type=str, default='None',
                        help='The profile yaml to use')
    parser.add_argument('--obsdir', type=str, default='*',
                        help='Observation directory name(s) separated by '
                             'commas')
    # test mode - do not run apero recipes just test
    parser.add_argument('--test', type=bool, default=False,
                        help='Do not run apero recipes just test')
    # link switch
    parser.add_argument('--links', type=bool, default=True)
    parser.add_argument('--only_links', type=bool, default=False)
    # precheck switch
    parser.add_argument('--precheck', type=bool, default=True)
    parser.add_argument('--only_precheck', type=bool, default=False)
    # process switch
    parser.add_argument('--process', type=bool, default=True)
    parser.add_argument('--only_process', type=bool, default=False)
    # apero get switch
    parser.add_argument('--get', type=bool, default=True)
    parser.add_argument('--only_get', type=bool, default=False)
    # apero reduction interface switch
    parser.add_argument('--ari', type=bool, default=True)
    parser.add_argument('--only_ari', type=bool, default=False)
    # apero get --since parameter
    parser.add_argument('--since', type=str, default='None',
                        help='APERO get - only copy files processed since this date')
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
    if args.obsdir in ['None', 'Null', '', None, '*']:
        obs_dirs = '*'
    else:
        obs_dirs = args.obsdir.split(',')
    # ----------------------------------------------------------------------
    # add this to settings
    settings['OBS_DIRS'] = obs_dirs
    # add the test mode
    settings['TEST'] = args.test
    # deal with since parameter
    if args.since in [None, 'None', 'Null']:
        settings['SINCE'] = None
    else:
        settings['SINCE'] = args.since
    # get switches
    settings['MAKELINKS'] = args.links
    settings['ONLYLINKS'] = args.only_links
    settings['PRECHECK'] = args.precheck
    settings['ONLYPRECHECK'] = args.only_precheck
    settings['PROCESSING'] = args.process
    settings['ONLYPROCESSING'] = args.only_process
    settings['APERO_GET'] = args.get
    settings['ONLYAPEROGET'] = args.only_get
    settings['REDUCTION_INTERFACE'] = args.ari
    settings['ONLYREDUCTIONINTERFACE'] = args.only_ari
    # deal with only switches (turn off all other switches)
    if settings['ONLYLINKS']:
        settings['MAKELINKS'] = True
        settings['PRECHECK'] = False
        settings['PROCESSING'] = False
        settings['APERO_GET'] = False
        settings['REDUCTION_INTERFACE'] = False
    if settings['ONLYPRECHECK']:
        settings['MAKELINKS'] = False
        settings['PRECHECK'] = True
        settings['PROCESSING'] = False
        settings['APERO_GET'] = False
        settings['REDUCTION_INTERFACE'] = False
    if settings['ONLYPROCESSING']:
        settings['MAKELINKS'] = False
        settings['PRECHECK'] = False
        settings['PROCESSING'] = True
        settings['APERO_GET'] = False
        settings['REDUCTION_INTERFACE'] = False
    if settings['ONLYAPEROGET']:
        settings['MAKELINKS'] = False
        settings['PRECHECK'] = False
        settings['PROCESSING'] = False
        settings['APERO_GET'] = True
        settings['REDUCTION_INTERFACE'] = False
    if settings['ONLYREDUCTIONINTERFACE']:
        settings['MAKELINKS'] = False
        settings['PRECHECK'] = False
        settings['PROCESSING'] = False
        settings['APERO_GET'] = False
        settings['REDUCTION_INTERFACE'] = True
    # ----------------------------------------------------------------------
    # read the yaml file and push into settings
    settings = read_yaml(args.profile, settings)
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
    os.environ['DRS_UCONFIG'] = profile['general']['apero profile']
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
        # get the profile dictionary
        pdict = settings['PROFILES'][profile]
        # print progress
        print(f'\tRunning profile: {profile}')
        # update the apero profile
        params = update_apero_profile(pdict)
        # get raw directory path from profile
        inpath = pdict['general']['raw dir']
        # get the raw directory from params
        outpath = params['DRS_DATA_RAW']
        # get obs dirs
        obs_dirs = settings['OBS_DIRS']
        # deal with getting all obs_dirs
        if obs_dirs == '*':
            obs_dirs = get_obs_dirs(pdict)
        # loop around obs dirs
        for obs_dir in obs_dirs:
            # get the full path
            full_inpath = os.path.join(inpath, obs_dir)
            full_outpath = os.path.join(outpath, obs_dir)

            # deal with test mode
            if settings['TEST']:
                print(f'\t\tSkipping symlink {full_outpath} [TEST]')
                continue
            # check if rawlink exists
            if not os.path.exists(full_inpath):
                print(f'\t\tSkipping {full_inpath} [does not exist]')
                continue

            # print symlink creation
            print(f'\t\tCreating symlink {full_outpath}')
            # check if the symlink exists
            if os.path.islink(full_outpath):
                continue
            # check if directory is a real link
            if os.path.exists(full_outpath):
                continue

            # make the symlink
            os.symlink(full_inpath, full_outpath)


def get_obs_dirs(profile):
    # get raw directory path from profile
    path = profile['general']['raw dir']
    # set start date
    start_date = str(profile['general']['start date'])
    # store obs dirs
    obs_dirs = []
    # convert start date to a astropy time
    start = Time(start_date)
    end = Time.now()
    # get time delta
    delta = (end - start).to(uu.day).round()
    # get a list of days
    days = start + np.arange(delta.value) * uu.day
    # loop around days and check for directory in path
    for day in days:
        # get the directory name
        obs_dir = day.iso.split(' ')[0]
        # check if the directory exists
        if os.path.exists(os.path.join(path, obs_dir)):
            obs_dirs.append(obs_dir)
    # return the obs dirs
    return obs_dirs


def run_precheck(settings: Dict[str, Any]):
    """
    Run the precheck on all profiles

    :param settings: dict, settings dictionary
    """
    _ = settings
    print('\tNot implemented.')
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
        # get the yaml dictionary for this profile
        pdict = settings['PROFILES'][profile]
        # update the apero profile
        _ = update_apero_profile(pdict)
        # get the run file
        runfile = pdict['processing']['run file']
        # get the obs dirs
        obs_dirs = settings['OBS_DIRS']
        # need to import apero_processing
        from apero.tools.recipes.bin import apero_processing
        # run apero processing
        if obs_dirs == '*':
            apero_processing.main(runfile=runfile,
                                  test=settings['TEST'])
        else:
            apero_processing.main(runfile=runfile,
                                  include_obs_dirs=','.join(obs_dirs),
                                  test=settings['TEST'])


def run_apero_get(settings: Dict[str, Any]):
    """
    Run the apero get on all profiles

    :param settings: dict, settings dictionary
    """
    from apero.core import constants
    # loop around profiles
    for profile in settings['PROFILES']:
        # print progress
        print(f'\tRunning profile: {profile}')
        # get the yaml dictionary for this profile
        pdict = settings['PROFILES'][profile]
        # update the apero profile
        pparams = update_apero_profile(pdict)
        pconst = constants.pload()
        # get the output types
        red_outtypes = ','.join(pdict['get']['science out types'])
        lbl_outtypes = ','.join(pdict['get-lbl']['science out types'])
        # get the dpr types
        red_dprtypes = ','.join(pdict['get']['science dpr types'])
        lbl_dprtypes = ','.join(pdict['get-lbl']['science dpr types'])
        simfp_dprtypes = ','.join(pdict['get-lbl']['simfp dprtypes'])
        # get fibers from settings
        scifibers, calfiber = pconst.FIBER_KINDS()
        calfibers = [calfiber]
        # template output types
        template_outtypes = pdict['get-lbl']['template out types']
        # get the object dir in the apero reduction path
        red_path = pparams['DRS_DATA_REDUC']
        obj_path = os.path.join(os.path.dirname(red_path), 'objects')
        # get the output path
        lbl_in_path = pdict['general']['lbl path']
        outpath_objects = os.path.join(lbl_in_path, 'science')
        outpath_templates = os.path.join(lbl_in_path, 'templates')
        outpath_calib = os.path.join(lbl_in_path, 'calib')
        outpath_fp = os.path.join(lbl_in_path, 'science/FP')
        # whether we want symlinks
        red_symlinks = pdict['get']['symlinks']
        lbl_symlinks = pdict['get-lbl']['symlinks']
        # ----------------------------------------------------------
        # reset reduced directories
        # ----------------------------------------------------------
        if red_symlinks:
            directories = [outpath_objects]
            dtypes = ['objects']
            reset = pdict['get']['reset']
            for it, directory in enumerate(directories):
                if dtypes[it] is None:
                    continue
                if reset is None:
                    continue
                if dtypes[it] in reset:
                    reset_directory(directory)
        # ----------------------------------------------------------
        # reset lbl directories
        # ----------------------------------------------------------
        if lbl_symlinks:
            directories = [obj_path, outpath_templates, outpath_calib]
            dtypes = ['science', 'templates', 'calib']
            reset = pdict['get-lbl']['reset']
            for it, directory in enumerate(directories):
                if dtypes[it] is None:
                    continue
                if reset is None:
                    continue
                if dtypes[it] in reset:
                    reset_directory(directory)

        # ---------------------------------------------------------------------
        # need to import apero_get (for this profile)
        from apero.tools.recipes.bin import apero_get
        # --------------------------------------------------------------
        # Copy to reduced 'objects' directory
        # --------------------------------------------------------------
        # run apero get to make the objects dir in apero dir
        apero_get.main(objnames='*', dprtypes=red_dprtypes,
                       outtypes=red_outtypes, outpath=obj_path,
                       fibers=scifibers, symlinks=red_symlinks,
                       test=settings['TEST'], since=settings['SINCE'])
        # --------------------------------------------------------------
        # Copy to LBL directory
        # --------------------------------------------------------------
        # run apero get for objects for lbl
        apero_get.main(objnames='*', dprtypes=lbl_dprtypes,
                       outtypes=lbl_outtypes,
                       outpath=outpath_objects, fibers=scifibers,
                       symlinks=lbl_symlinks,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for templates (no DPRTYPE as they could be different)
        apero_get.main(objnames='*', outtypes=template_outtypes,
                       outpath=outpath_templates, fibers=scifibers,
                       symlinks=False, nosubdir=True,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for simultaneous FP
        apero_get.main(objnames='*', dprtypes=simfp_dprtypes,
                       outtypes='EXT_E2DS_FF', nosubdir=True,
                       outpath=outpath_fp, fibers=calfibers,
                       symlinks=lbl_symlinks,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for extracted FP_FP
        apero_get.main(objnames='*', dprtypes='FP_FP',
                       outtypes='EXT_E2DS_FF',
                       outpath=outpath_fp, fibers=calfibers,
                       symlinks=lbl_symlinks, nosubdir=True,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for calibs (wave + blaze) science fiber
        apero_get.main(objnames='*',
                       outtypes='FF_BLAZE,WAVE_NIGHT',
                       outpath=outpath_calib, fibers=scifibers,
                       symlinks=lbl_symlinks, nosubdir=True,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for calibs (wave + blaze) science fiber
        apero_get.main(objnames='*',
                       outtypes='FF_BLAZE,WAVE_NIGHT',
                       outpath=outpath_calib, fibers=calfibers,
                       symlinks=lbl_symlinks, nosubdir=True,
                       test=settings['TEST'], since=settings['SINCE'])


def reset_directory(directory: str, reset: bool = False):
    # if directory does not exist create it
    if not os.path.exists(directory):
        os.makedirs(directory)
        return
    # if we are here and want to reset we should remove the directory tree
    if reset:
        shutil.rmtree(directory)
        if not os.path.exists(directory):
            os.makedirs(directory)


def run_apero_reduction_interface(settings: Dict[str, Any]):
    """
    Run the apero reduction interface

    :param settings: dict, settings dictionary
    """
    _ = settings
    print('\tNot implemented.')
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
    if trigger_settings['MAKELINKS']:
        print_process('Making symbolic links')
        make_sym_links(trigger_settings)
        # deal with only creating links
        if trigger_settings['ONLYLINKS']:
            print_process('Only making symbolic links')
            sys.exit(0)
    # ----------------------------------------------------------------------
    # run apero precheck on all profiles
    if trigger_settings['PRECHECK']:
        print_process('Running apero precheck')
        run_precheck(trigger_settings)
        # deal with only running precheck
        if trigger_settings['ONLYPRECHECK']:
            print_process('Only running apero precheck')
            sys.exit(0)
    # ----------------------------------------------------------------------
    # run apero processing on all profiles
    if trigger_settings['PROCESSING']:
        print_process('Running apero processing')
        run_processing(trigger_settings)
        # deal with only running processing
        if trigger_settings['ONLYPROCESSING']:
            print_process('Only running apero processing')
            sys.exit(0)
    # ----------------------------------------------------------------------
    # run apero get on all profiles
    if trigger_settings['APERO_GET']:
        print_process('Running apero get')
        run_apero_get(trigger_settings)
    # ----------------------------------------------------------------------
    # run the apero reduction interface
    if trigger_settings['REDUCTION_INTERFACE']:
        print_process('Running apero reduction interface')
        run_apero_reduction_interface(trigger_settings)


# =============================================================================
# End of code
# =============================================================================
