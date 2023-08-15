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
from typing import Any, Dict, List

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
# define messages
MANUAL_START = 'MANUAL_START'
MANUAL_END = 'MANUAL_END'
APERO_START = 'APERO_START'
APERO_ERR = 'APERO_ERR'
APERO_END = 'APERO_END'
ARI_START = 'ARI_START'
ARI_END = 'ARI_END'


MESSAGES = [MANUAL_START, MANUAL_END, APERO_START, APERO_ERR, APERO_END,
            ARI_START, ARI_END]

# =============================================================================
# Define classes
# =============================================================================
class TriggerLog:
    def __init__(self, filename, profile, obsdirs):
        self.filename = filename
        self.profile = profile
        self.obsdirs = obsdirs

    def write(self, logkind: str, message: str = 'None'):
        # log kind must be correct
        if logkind not in MESSAGES:
            raise ValueError(f'logkind must be one of {MESSAGES}')
        # set obs dir string
        obsdir_str = '|'.join(self.obsdirs)
        # get time now
        timenow = Time.now()
        # replace any double speech marks with single ones
        message = message.replace('"', "'")
        # construct the line
        elements = [timenow.fits, self.profile, logkind, f"{obsdir_str}",
                    f'"{message}"']
        # get line as single string
        line = ', '.join(elements)
        # open file and append line
        with open(self.filename, 'a') as logfile:
            logfile.write(line + '\n')



# =============================================================================
# Define functions
# =============================================================================
def get_args():
    """
    Define the command line arguments

    :return: argparse namespace object containing arguments
    """
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
    # apero process switch
    parser.add_argument('--apero_process', type=bool, default=True)
    parser.add_argument('--only_apero_process', type=bool, default=False)
    # apero get switch
    parser.add_argument('--get', type=bool, default=True)
    parser.add_argument('--only_aperoget', type=bool, default=False)
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
    settings['ONLY_LINKS'] = args.only_links

    settings['APERO_PROCESSING'] = args.apero_process
    settings['ONLY_APEROPROCESSING'] = args.only_apero_process

    settings['APERO_GET'] = args.get
    settings['ONLY_APEROGET'] = args.only_aperoget

    settings['REDUCTION_INTERFACE'] = args.ari
    settings['ONLY_REDUCTIONINTERFACE'] = args.only_ari
    # -------------------------------------------------------------------------
    # deal with only switches (turn off all other switches)
    # -------------------------------------------------------------------------
    # only links
    if settings['ONLY_LINKS']:
        settings['MAKELINKS'] = True
        settings['APERO_PROCESSING'] = False
        settings['APERO_GET'] = False
        settings['REDUCTION_INTERFACE'] = False
    # only APERO processing
    if settings['ONLY_APEROPROCESSING']:
        settings['MAKELINKS'] = False
        settings['APERO_PROCESSING'] = True
        settings['APERO_GET'] = False
        settings['REDUCTION_INTERFACE'] = False
    # only apero get
    if settings['ONLY_APEROGET']:
        settings['MAKELINKS'] = False
        settings['APERO_PROCESSING'] = False
        settings['APERO_GET'] = True
        settings['REDUCTION_INTERFACE'] = False
    # only ARI
    if settings['ONLY_REDUCTIONINTERFACE']:
        settings['MAKELINKS'] = False
        settings['APERO_PROCESSING'] = False
        settings['APERO_GET'] = False
        settings['REDUCTION_INTERFACE'] = True
    # ----------------------------------------------------------------------
    # read the yaml file and push into settings
    settings = read_yaml(args.profile, settings)
    # ----------------------------------------------------------------------
    # construct log path
    homedir = os.path.expanduser('~')
    # construct log path
    logpath = os.path.join(homedir, '.apero', 'manual_trigger')
    # make sure logpath exists
    if not os.path.exists(logpath):
        os.makedirs(logpath)
    # construct log filename
    logfile = os.path.join(logpath, args.profile.replace('.yaml', '.log'))
    # make a log for each profile
    settings['LOG'] = dict()
    # loop around profiles
    for profile in settings['PROFILES']:
        # get the obs dirs
        obs_dirs = settings['OBS_DIRS']
        # get log class
        logclass = TriggerLog(logfile, profile=profile, obsdirs=obs_dirs)
        # push into settings
        settings['LOG'][profile] = logclass
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


def get_obs_dirs(profile: Dict[str, Any]) -> List[str]:
    """
    Get the observation directories from the raw directory
    This is used when no obs_dirs are specified in the yaml file

    :param profile: dict, the profile dictionary from the yaml file

    :return: list of observation directories
    """
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


def run_processing(settings: Dict[str, Any]):
    """
    Run the processing on all profiles

    :param settings: dict, settings dictionary
    """
    # loop around profiles
    for profile in settings['PROFILES']:
        # get the obs dirs
        obs_dirs = settings['OBS_DIRS']
        # get the obsdirs as a string
        obsdir_str = ','.join(obs_dirs)
        # log that APERO started
        trigger_settings['LOG'][profile].write(APERO_START)
        # print progress
        print(f'\tRunning profile: {profile}')
        # get the yaml dictionary for this profile
        pdict = settings['PROFILES'][profile]
        # update the apero profile
        params = update_apero_profile(pdict)
        # deal with a reset
        if pdict['processing']['reset'] not in ['None', None]:
            # ask user because this is dangerous
            msg = 'Are you sure you want to reset {0}? [Y]es/[N]o: '
            uinput = input(msg.format(pdict['processing']['reset']))
            if uinput.upper() == 'Y' or uinput.upper() == 'YES':
                apero_reset(params, pdict)
        # get the run file
        runfile = pdict['processing']['run file']

        # need to import apero_processing
        from apero.tools.recipes.bin import apero_processing
        # run apero processing
        if obs_dirs == '*':
            ll = apero_processing.main(runfile=runfile,
                                       test=settings['TEST'])
        else:
            ll = apero_processing.main(runfile=runfile,
                                       include_obs_dirs=obsdir_str,
                                       test=settings['TEST'])
        # TODO: log apero errors using ll
        _ = ll
        # log that APERO ended
        trigger_settings['LOG'][profile].write(APERO_END)



def apero_reset(params: Any, pdict: Dict[str, Any]):
    """
    Reset the processing folders

    :param params: ParamDict, parmaeters dictionary of constants for this
                   profile
    :param pdict: Dict, the yaml dictionary for this profile

    :return: None, removes files and updates database
    """
    # import the drs_reset module
    from apero.tools.module.setup import drs_reset
    from apero.tools.recipes.bin import apero_reset
    # get the reset directories
    reset_dirs = pdict['processing']['reset']
    # define the function to use for each type (via a dictionary)
    reset_funcs = dict()
    reset_funcs['tmp'] = drs_reset.reset_tmp_folders
    reset_funcs['red'] = drs_reset.reset_reduced_folders
    reset_funcs['out'] = drs_reset.reset_out_folders
    # deal with None in reset dirs
    if 'None' in reset_dirs:
        return
    # print progress
    print('\t\tResetting processing')
    # deal with all in reset dirs
    if 'all' in reset_dirs:
        apero_reset.main(warn=False)
    # loop around directories to be reset
    for reset_dir in reset_dirs:
        if reset_dir in reset_funcs:
            reset_funcs[reset_dir](params)


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
        # scifiber and calfiber must be strings (comma separated)
        red_scifibers = ','.join(pdict['get']['science fibers'])
        # red_calfibers = ','.join(pdict['get']['calib fibers'])
        lbl_scifibers = ','.join(pdict['get-lbl']['science fibers'])
        lbl_calfibers = ','.join(pdict['get-lbl']['calib fibers'])
        # template output types
        red_template_outtypes = ','.join(pdict['get']['template out types'])
        lbl_template_outtypes = ','.join(pdict['get-lbl']['template out types'])
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
        # check directories exist - try to make them if they don't
        # ----------------------------------------------------------
        directories = [obj_path, outpath_templates, outpath_calib,
                       outpath_objects]
        for directory in directories:
            if not os.path.exists(directory):
                os.makedirs(directory)
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
                       fibers=red_scifibers, symlinks=red_symlinks,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for templates (no DPRTYPE as they could be different)
        apero_get.main(objnames='*', outtypes=red_template_outtypes,
                       outpath=obj_path, fibers=red_scifibers,
                       symlinks=red_symlinks,
                       test=settings['TEST'], since=settings['SINCE'])
        # --------------------------------------------------------------
        # Copy to LBL directory
        # --------------------------------------------------------------
        # run apero get for objects for lbl
        apero_get.main(objnames='*', dprtypes=lbl_dprtypes,
                       outtypes=lbl_outtypes,
                       outpath=outpath_objects, fibers=lbl_scifibers,
                       symlinks=lbl_symlinks,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for templates (no DPRTYPE as they could be different)
        apero_get.main(objnames='*', outtypes=lbl_template_outtypes,
                       outpath=outpath_templates, fibers=lbl_scifibers,
                       symlinks=False, nosubdir=True,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for simultaneous FP
        apero_get.main(objnames='None', dprtypes=simfp_dprtypes,
                       outtypes='EXT_E2DS_FF', nosubdir=True,
                       outpath=outpath_fp, fibers=lbl_calfibers,
                       symlinks=lbl_symlinks,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for extracted FP_FP
        apero_get.main(objnames='None', dprtypes='FP_FP',
                       outtypes='EXT_E2DS_FF',
                       outpath=outpath_fp, fibers=lbl_calfibers,
                       symlinks=lbl_symlinks, nosubdir=True,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for calibs (wave + blaze) science fiber
        apero_get.main(objnames='None', outtypes='FF_BLAZE,WAVE_NIGHT',
                       outpath=outpath_calib, fibers=lbl_scifibers,
                       symlinks=lbl_symlinks, nosubdir=True,
                       test=settings['TEST'], since=settings['SINCE'])
        # run apero get for calibs (wave + blaze) science fiber
        apero_get.main(objnames='None',
                       outtypes='FF_BLAZE,WAVE_NIGHT',
                       outpath=outpath_calib, fibers=lbl_calfibers,
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
    # get the current working directory
    cwd = os.getcwd()
    # loop around profiles
    for profile in settings['PROFILES']:
        # get ari path
        ari_path = settings[profile]['ari']['path']
        # get ari profile
        ari_profile = settings[profile]['ari']['profile']
        # deal with no ari code
        if 'simple_ari.py' not in os.listdir(ari_path):
            print('\t\tERROR: simple_ari.py not found in {0}'.format(ari_path))
            continue
        # log that APERO started
        settings['LOG'][profile].write(ARI_START)
        # change to ari path
        os.chdir(ari_path)
        # set up command
        ari_cmd = 'python simple_ari.py {0} --filter="{1}"'
        # run simple ari interface
        # TODO: This is terrible - do not use os.system
        os.system(ari_cmd.format(ari_profile, profile))
        # log that APERO started
        settings['LOG'][profile].write(ARI_END)
    # change back to original path
    os.chdir(cwd)


def run_lbl_processing(settings: Dict[str, Any]):
    """
    Run the LBL processing

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
    # log that we have started manual trigger
    for profile in trigger_settings['PROFILES']:
        trigger_settings['LOG'][profile].write(MANUAL_START)
    # ----------------------------------------------------------------------
    # make symbolic links
    if trigger_settings['MAKELINKS']:
        print_process('Making symbolic links')
        make_sym_links(trigger_settings)
        # deal with only creating links
        if trigger_settings['ONLY_LINKS']:
            print_process('Only making symbolic links')
            sys.exit(0)
    # ----------------------------------------------------------------------
    # run apero processing on all profiles
    if trigger_settings['APERO_PROCESSING']:
        print_process('Running apero processing')
        run_processing(trigger_settings)
        # deal with only running processing
        if trigger_settings['ONLY_APEROPROCESSING']:
            print_process('Only running apero processing')
            sys.exit(0)
    # ----------------------------------------------------------------------
    # run apero get on all profiles
    if trigger_settings['APERO_GET']:
        print_process('Running apero get')
        run_apero_get(trigger_settings)
        # deal with only running processing
        if trigger_settings['ONLY_APEROGET']:
            print_process('Only running apero get')
            sys.exit(0)
    # ----------------------------------------------------------------------
    # run the apero reduction interface
    if trigger_settings['REDUCTION_INTERFACE']:
        print_process('Running apero reduction interface')
        # TODO: flag only to do this night
        run_apero_reduction_interface(trigger_settings)
        # deal with only running processing
        if trigger_settings['ONLY_REDUCTIONINTERFACE']:
            print_process('Only running apero get')
            sys.exit(0)
    # ----------------------------------------------------------------------
    # log that we have finished
    # log that we have started manual trigger
    for profile in trigger_settings['PROFILES']:
        trigger_settings['LOG'][profile].write(MANUAL_END)

# =============================================================================
# End of code
# =============================================================================
