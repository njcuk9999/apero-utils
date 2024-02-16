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

import git
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


def update_apero_profile(params: Dict[str, Any], profile: int) -> Any:
    """
    Update the apero profile with the correct paths

    :param params: dict, the parameters from the yaml file
    :param profile: int, the profile (either 1 or 2)
    :return:
    """
    # deal with profile 1 or profile 2
    if profile == 1:
        profile_path = params['profile1']
        install_path = params.get('apero install 1', None)
    elif profile == 2:
        profile_path = params['profile2']
        install_path = params.get('apero install 2', None)
    else:
        emsg = 'profile must be 1 or 2'
        raise AperoCopyError(emsg)
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = profile_path
    # allow getting apero
    if install_path is not None:
        sys.path.append(install_path)
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
                        action='store_true', default=False)
    parser.add_argument('--test', help='run in test mode (no copy or remove)',
                        action='store_true', default=False)
    parser.add_argument('--symlinks', action='store_true', default=False,
                        help='Copy in symlink mode (not recommended)')
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
    apero_params = update_apero_profile(params, profile=1)
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
    apero_params = update_apero_profile(params, profile=2)
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


def fast_copy(params, files1: Dict[str, List[str]],
              paths1: Dict[str, str]):

    # deal with symlinks --> return
    if params['symlinks']:
        return False
    # flag for fast copy
    did_fast_copy = True
    # load parameters from profile 2
    apero_params = update_apero_profile(params, profile=2)
    # import database
    from apero.core.constants import path_definitions
    # get block definitions
    blocks = path_definitions.BLOCKS
    # loop around block definitions
    for block in blocks:
        # skip if not in files1
        if block.name not in files1:
            continue
        # get block instance
        block_inst = block(apero_params)
        # get block name
        block_name = block.name
        # get path1 from paths1
        path1 = paths1[block_name]
        # get block path
        path2 = block_inst.path
        # get path3 the temporary path to copy things to first
        block_base = os.path.basename(path2)
        block_dir = os.path.dirname(path2)

        # get the input and output paths
        inpath = str(path1)
        outpath = os.path.join(block_dir, TMP_DIR_PREFIX + block_base)
        # check if the input path exists
        if not os.path.exists(outpath):
            # create path
            os.makedirs(outpath)
        # now add os.sep to inpath and outpath
        if not inpath.endswith(os.sep):
            inpath += os.sep
        if not outpath.endswith(os.sep):
            outpath += os.sep
        # get the rsync command
        rsync_cmd = 'rsync -avu {0} {1}'.format(inpath, outpath)
        # print progress
        print('\n\n' + '-' * 50)
        print('Copying files for block={0}'.format(block_name))
        print('Using rsync command:')
        print('\t', rsync_cmd)
        print('-' * 50 + '\n')
        # run the command
        try:
            os.system(rsync_cmd)
        except Exception as e:
            print('Skipping rsync. Error: {0}'.format(e))
            did_fast_copy = False
    # return whether we did fast copy for all directories
    return did_fast_copy


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
                # need to deal with directory not existing
                if not os.path.exists(os.path.dirname(outfilename3)):
                    os.makedirs(os.path.dirname(outfilename3))
                # deal with tmp file existing
                if os.path.exists(outfilename3):
                    msg = ('Skipping {0} (tmp file already exists)')
                    margs = [outfilename3]
                    print(msg.format(*margs))
                    continue
                # copy file
                if params['symlinks']:
                    os.symlink(infilename, outfilename3)
                else:
                    shutil.copy(infilename, outfilename3)


def update_git_profile2(params: Dict[str, Any]):
    # get install path for profile 1
    install_path1 = params.get('apero install 1', None)
    # get install path for profile 2
    install_path2 = params.get('apero install 2', None)
    # deal with no install path (something is installed with pip)
    if install_path1 is None or install_path2 is None:
        return
    # deal with install paths being the same
    # in this case we don't need to do anything
    if install_path1 == install_path2:
        return
    # print that we are updating profile 2
    print('Using git to update profile 2 to match profile 1')
    # get the git repo for profile 1
    repo1 = git.Repo(install_path1)
    # get the git repo for profile 2
    repo2 = git.Repo(install_path2)
    # get the active branch for profile 1
    active1 = repo1.active_branch.name
    # move to that branch in repo2
    msg = repo2.git.checkout(active1)
    print(msg)
    # update repo2
    msg = repo2.git.pull()
    print(msg)


def reset_profile2(params: Dict[str, Any], paths3: Dict[str, str]):
    # if any block paths in path3 don't exist do not reset anything
    for block_name in paths3:
        if not os.path.exists(paths3[block_name]):
            print('Error {0} does not exist'.format(paths3[block_name]))
            return False
    # load parameters from profile 1
    apero_params = update_apero_profile(params, profile=2)
    # import modules
    from apero.tools.module.setup import drs_reset
    # get database timeout
    database_timeout = 0
    # assets folder
    drs_reset.reset_assets(apero_params, dtimeout=database_timeout)
    # tmp folder
    drs_reset.reset_tmp_folders(apero_params, dtimeout=database_timeout)
    # reduced folder
    drs_reset.reset_reduced_folders(apero_params, dtimeout=database_timeout)
    # calibration folder
    drs_reset.reset_calibdb(apero_params, dtimeout=database_timeout)
    # telluric folder
    drs_reset.reset_telludb(apero_params, dtimeout=database_timeout)
    # log folder
    drs_reset.reset_log(apero_params)
    # plot folder
    drs_reset.reset_plot(apero_params)
    # run folder
    drs_reset.reset_run(apero_params)
    # out folder
    drs_reset.reset_out_folders(apero_params, dtimeout=database_timeout)
    # reset correct
    return True


def rename_directories(params: Dict[str, Any], paths2: Dict[str, str],
                       paths3: Dict[str, str]):
    # loop around each block
    for block_name in paths2:
        if not os.path.exists(paths3[block_name]):
            print('Error {0} does not exist'.format(paths3[block_name]))
            return False
        msg = 'Renaming {0}: ({1} -> {2})'
        margs = [block_name, paths3[block_name], paths2[block_name]]
        print(msg.format(*margs))
        if not params['test']:
            if os.path.exists(paths2[block_name]):
                if os.path.islink(paths2[block_name]):
                    os.remove(paths2[block_name])
                else:
                    shutil.rmtree(paths2[block_name])
            os.rename(paths3[block_name], paths2[block_name])
    return True


def update_databases_profile2(params):
    # print progress
    print('Updating APERO databases')
    # if test mode do nothing here
    if params['test']:
        return
    # load parameters from profile 1
    apero_params1 = update_apero_profile(params, profile=1).copy()
    apero_params2 = update_apero_profile(params, profile=2).copy()
    # apero imports
    from apero.tools.module.database import database_update
    # update the databases
    copy_database(apero_params1, apero_params2)


def copy_database(params_old: Any, params_new: Any):
    """
    Copy the calib/tellu/log and index databases from one location to another

    :param params:

    :return:
    """
    print('stop')
    # get table names for each old table

    # get table names for each new table

    # get old block paths

    # get new block paths

    # delete new tables

    # duplicate all databases

    # loop around block kinds

    # Deal with findex database



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

    # try a fast copy using rsync
    did_fast_copy = fast_copy(params, files1, paths1)
    # copy files from profile 1 to profile 2 for each block kind
    # must copy files to a temporary path first (copying can be slow)
    if not did_fast_copy:
        copy_files(params, files1, files2, files3)
    # ask after copying whether we are ready to continue
    uinput = ''
    while uinput.lower() not in ['yes', 'no']:
        uinput = input('Are you ready to remove old data and replace with '
                       'the new data. Note this is a one way process and there '
                       'is no undo \n\tType "yes" or "no"?\t')
    if uinput.lower() ==  'no':
        print('\n\tStopping copy.')
        return
    # may need to update profile 2 (via git) to match profile 1
    update_git_profile2(params)
    # remove all old files from profile 2 blocks
    success = reset_profile2(params, paths3)
    # rename the directories in profile 2 (this is quicker than copying)
    if success:
        success = rename_directories(params, paths2, paths3)
    # update databases for profile 2
    if success:
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
