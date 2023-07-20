#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-13 at 13:07

@author: cook
"""
import argparse
import os
import shutil
import sys
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import yaml
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
__version__ = '0.0.1'
# prefix for TMP directory to copy to (before renaming)
TMP_DIR_PREFIX = 'APERO_COPY_TMP_'


# -----------------------------------------------------------------------------


# =============================================================================
# Define functions
# =============================================================================
class AperoCopyError(Exception):
    """
    Error class for apero_raw_tests
    """
    pass


def update_apero_profile(params: Dict[str, Any], profile_path: str) -> Any:
    """
    Update the apero profile with the correct paths

    :param params: dict, the parameters from the yaml file
    :param profile_path: dict, the profile to update
    :return:
    """
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = profile_path
    # allow getting apero
    if 'apero install' in params and params['apero install'] is not None:
        sys.path.append(params['apero install'])
    # load apero modules
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    from apero.core.utils import drs_startup
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    apero_params = constants.load(cache=False)
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # set apero pid
    apero_params['PID'], apero_params['DATE_NOW'] = drs_startup.assign_pid()
    # no inputs
    apero_params['INPUTS'] = dict()
    apero_params['OBS_DIR'] = None
    # make sure parameters is reloaded (and not cached)
    return apero_params


def load_params():
    """
    Load parameters from yaml file
    :return:
    """
    # start arg parser
    parser = argparse.ArgumentParser(description='Apero Copy')
    # add arguments
    parser.add_argument('yaml', help='yaml file to use', type=str,
                        default=None)
    parser.add_argument('--overwrite', help='overwrite existing files',
                        default=False, type=bool)
    parser.add_argument('--test', help='run in test mode (no copy or remove)',
                        default=False, type=bool)
    # get arguments
    args = vars(parser.parse_args())
    # ------------------------------------------------------------------
    # deal with getting yaml file
    # ------------------------------------------------------------------
    # deal with no yaml file set
    if args['yaml'] is None:
        emsg = 'yaml file must be set'
        raise AperoCopyError(emsg)
    # deal with bad path
    if not os.path.exists(args['yaml']):
        emsg = 'yaml file {0} does not exist'
        eargs = [args['yaml']]
        raise AperoCopyError(emsg.format(*eargs))
    # otherwise we read the yaml file
    params = read_yaml(args['yaml'])
    # add args to params
    for arg in args:
        if arg not in params:
            params[arg] = args[arg]
    # ------------------------------------------------------------------
    # return the parameters
    return params


def read_yaml(yaml_filename: Union[str, None]) -> Dict[str, Any]:
    """
    Read the yaml file and add to settings

    :param yaml_filename: str, yaml file name

    :return: dict, updated settings dictionary
    """
    # deal with yaml_filename being None
    if yaml_filename is None:
        emsg = 'yaml_filename must be set to a valid file'
        raise AperoCopyError(emsg)
    # deal with yaml_filename not existing
    if not os.path.exists(yaml_filename):
        emsg = 'yaml_filename {0} does not exist'
        eargs = [yaml_filename]
        raise AperoCopyError(emsg.format(*eargs))
    # read the yaml file
    with open(yaml_filename, 'r') as f:
        yaml_data = yaml.load(f, Loader=yaml.FullLoader)
    # add a profiles sub-dictionary
    settings = dict()
    # loop around yaml data
    for key, value in yaml_data.items():
        # if key is in settings
        settings[key] = value
    # return settings
    return settings


def get_files_profile1(params) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    # load parameters from profile 1
    apero_params = update_apero_profile(params, params['profile1'])
    # import database
    from apero.core.core import drs_database
    from apero.core.constants import path_definitions
    # get block definitions
    blocks = path_definitions.BLOCKS
    # get file index database
    findexdb = drs_database.FileIndexDatabase(apero_params)
    findexdb.load_db()
    calibdb = drs_database.CalibrationDatabase(apero_params)
    calibdb.load_db()
    telludb = drs_database.TelluricDatabase(apero_params)
    telludb.load_db()
    # storage for for files
    files1 = dict()
    paths1 = dict()
    # loop around block definitions
    for block in blocks:
        # get block instance
        block_inst = block(apero_params)
        # get block
        block_name = block_inst.name
        # get block path
        block_path = block_inst.path
        # do not copy raw files
        if block_name == 'raw':
            continue
        # deal with calib and tellu blocks
        if block_inst.name in ['calib', 'tellu']:
            files1[block_name] = []
            paths1[block_name] = block_inst.path
            continue
        # only keep those that have an observation directory
        if not block_inst.has_obs_dirs:
            continue
        # condition for database
        condition = 'BLOCK_KIND="{0}"'.format(block_name)
        # get files from database
        files1[block_name] = findexdb.get_entries(columns='ABSPATH',
                                                  condition=condition)
        # add to path dict
        paths1[block_name] = block_path

    # need to add the calibration and telluric database files differently
    for block_name in ['calib', 'tellu']:
        # deal with calibration database
        if block_name == 'calib':
            query = calibdb.database.get('FILENAME', return_pandas=True)
        # deal with telluric database
        else:
            query = telludb.database.get('FILENAME', return_pandas=True)
        # push query into a list of basenames (as a character array)
        basenames = np.char.array(list(query['FILENAME']))
        # add path to basenames
        filenames = paths1[block_name] + os.sep + basenames
        # get files from database
        files1[block_name] = filenames
    # add to path dict
    return files1, paths1


def get_files_profile2(params, files1: Dict[str, List[str]],
                       paths1: Dict[str, str]
                       ) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], Dict[str, str], Dict[str, str]]:
    # load parameters from profile 1
    apero_params = update_apero_profile(params, params['profile2'])
    # import database
    from apero.core.constants import path_definitions
    # get block definitions
    blocks = path_definitions.BLOCKS
    # storage for for files
    files2 = dict()
    files3 = dict()
    paths2 = dict()
    paths3 = dict()
    # loop around block definitions
    for block in blocks:
        # skip if not in files1
        if block.name not in files1:
            continue
        # get block instance
        block_inst = block(apero_params)
        # get block name
        block_name = block.name
        # get block path
        path2 = block_inst.path
        # get path3 the temporary path to copy things to first
        block_base = os.path.basename(path2)
        block_dir = os.path.dirname(path2)
        path3 = os.path.join(block_dir, TMP_DIR_PREFIX + block_base)
        # get path1 from paths1
        path1 = paths1[block_name]
        # add a list for each block kind to files
        files2[block_name] = []
        files3[block_name] = []
        # loop around files
        for infilename in files1[block_name]:
            # replace path1 with path2
            outfilename2 = infilename.replace(path1, path2)
            outfilename3 = infilename.replace(path1, path3)
            # add to files
            files2[block_name].append(outfilename2)
            files3[block_name].append(outfilename3)
        # add to path dict
        paths2[block_name] = path2
        paths3[block_name] = path3
    # return files
    return files2, files3, paths2, paths3


def copy_files(params, files1: Dict[str, List[str]],
               files2: Dict[str, List[str]],
               files3: Dict[str, List[str]]):
    # loop around each block kind and copy files
    for block_name in files1:

        # print progress
        print('\n\n' + '-' * 50)
        print('Copying files for block={0}'.format(block_name))
        print('-' * 50 + '\n')

        # get a list of iterations
        iterations = list(range(len(files1[block_name])))
        # loop around files
        for it in tqdm(iterations):
            # get infilename
            infilename = files1[block_name][it]
            # get outfilename
            outfilename2 = files2[block_name][it]
            outfilename3 = files3[block_name][it]
            # make sure in file exists (it really should)
            if not os.path.exists(infilename):
                emsg = 'Input file does not exist'
                emsg += '\n\t infilename = {0}'
                emsg += '\n\t tmpfilename = {1}'
                emsg += '\n\t outfilename = {1}'
                eargs = [infilename, outfilename3, outfilename2]
                raise AperoCopyError(emsg.format(*eargs))
            # check that infilename and outfilename are not the same
            if infilename == outfilename2 or infilename == outfilename3:
                emsg = 'Cannot copy file input and output filename are the same'
                emsg += '\n\t infilename = {0}'
                emsg += '\n\t tmpfilename = {1}'
                emsg += '\n\t outfilename = {1}'
                eargs = [infilename, outfilename3, outfilename2]
                raise AperoCopyError(emsg.format(*eargs))
            # copy file
            if params['test']:
                print('Copying {0} to {1}'.format(infilename, outfilename3))
            else:
                # deal with file existing
                if os.path.exists(outfilename2) and not params['overwrite']:
                    msg = ('Skipping {0} (already exists) use '
                           '--overwrite to copy')
                    margs = [outfilename2]
                    print(msg.format(*margs))
                    continue
                # deal with tmp file existing
                if os.path.exists(outfilename3):
                    msg = ('Skipping {0} (tmp file already exists) use '
                           '--overwrite to copy')
                    margs = [outfilename3]
                    print(msg.format(*margs))
                    continue
                # need to deal with no tmp directory
                if not os.path.exists(os.path.dirname(outfilename3)):
                    os.makedirs(os.path.dirname(outfilename3))
                # copy file
                shutil.copy(infilename, outfilename3)


def remove_files(params: Dict[str, Any], paths2: Dict[str, str]):
    if params['overwrite']:
        # loop around each block
        for block_name in paths2:
            msg = 'Removing {0}: ({1})'
            margs = [block_name, paths2[block_name]]
            print(msg.format(*margs))
            if not params['test']:
                shutil.rmtree(paths2[block_name])


def rename_directories(params: Dict[str, Any], paths2: Dict[str, str],
                       paths3: Dict[str, str]):
    # loop around each block
    for block_name in paths2:
        msg = 'Renaming {0}: ({1} -> {2})'
        margs = [block_name, paths3[block_name], paths2[block_name]]
        print(msg.format(*margs))
        if not params['test']:
            os.rename(paths3[block_name], paths2[block_name])


def update_databases_profile2(params):
    # print progress
    print('Updating APERO databases')
    # if test mode do nothing here
    if params['test']:
        return
    # load parameters from profile 1
    apero_params = update_apero_profile(params, params['profile2'])
    # apero imports
    from apero.tools.module.database import database_update
    # update the databases
    dbkind = 'all'
    database_update.update_database(apero_params, dbkind=dbkind)


def main():
    """
    Get files from profile 1, get blocks from profile 2
    :return:
    """
    # load parameters
    params = load_params()
    # check that both profiles exist (we assume that this means the directories
    # are found on disk)
    if not os.path.exists(params['profile1']):
        emsg = 'profile1={0} does not exist'
        eargs = [params['profile1']]
        raise AperoCopyError(emsg.format(*eargs))
    if not os.path.exists(params['profile2']):
        emsg = 'profile2={0} does not exist'
        eargs = [params['profile2']]
        raise AperoCopyError(emsg.format(*eargs))
    # get a list of files from profile 1 for each block kinds
    files1, paths1 = get_files_profile1(params)
    # get the output files for profile 2 for each block kind
    files2, files3, paths2, paths3 = get_files_profile2(params, files1, paths1)
    # copy files from profile 1 to profile 2 for each block kind
    # must copy files to a temporary path first (copying can be slow)
    copy_files(params, files1, files2, files3)
    # remove all old files from profile 2 blocks
    remove_files(params, paths2)
    # rename the directories in profile 2 (this is quicker than copying)
    rename_directories(params, paths2, paths3)
    # update databases for profile 2
    update_databases_profile2(params)
    # return to __main__
    return


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # retun main code
    main()

# =============================================================================
# End of code
# =============================================================================
