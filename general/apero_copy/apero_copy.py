#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-13 at 13:07

@author: cook
"""
import os
import argparse
from typing import Any, Dict, List, Tuple, Union
import shutil
import yaml
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
__version__ = '0.0.1'

TEST = True
# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
class AperoCopyError(Exception):
    """
    Error class for apero_raw_tests
    """
    pass


def update_apero_profile(profile_path: str) -> Any:
    """
    Update the apero profile with the correct paths
    :param profile: dict, the profile to update
    :return:
    """
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = profile_path
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # make sure parameters is reloaded (and not cached)
    return constants.load(cache=False)


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
    apero_params = update_apero_profile(params['profile1'])

    # import database
    from apero.core.core import drs_database
    from apero.core.constants import path_definitions
    # get block definitions
    blocks = path_definitions.BLOCKS
    # get file index database
    findexdb = drs_database.FileIndexDatabase(params['profile1'])
    # storage for for files
    files1 = dict()
    paths1 = dict()
    # loop around block definitions
    for block in blocks:
        # get block instance
        block_inst = block(apero_params)
        # only keep those that have an observation directory
        if not block_inst.has_obs_dirs:
            continue
        # get block
        block_name = block_inst.name
        # get block path
        block_path = block_inst.path
        # condition for database
        condition = 'BLOCK_KIND="{0}"'.format(block_name)
        # get files from database
        files1[block_name] = findexdb.get_entries(columns='ABSPATH',
                                                 condition=condition)
        # add to path dict
        paths1[block_name] = block_path

    return files1, paths1


def get_blocks_profile2(params, files1: Dict[str, List[str]],
                        paths1: Dict[str, str]) -> Dict[str, List[str]]:
    # load parameters from profile 1
    apero_params = update_apero_profile(params['profile2'])
    # import database
    from apero.core.constants import path_definitions
    # get block definitions
    blocks = path_definitions.BLOCKS
    # storage for for files
    files2 = dict()
    # loop around block definitions
    for block in blocks:
        # get block instance
        block_inst = block(apero_params)
        # only keep those that have an observation directory
        if not block_inst.has_obs_dirs:
            continue
        # get block name
        block_name = block.name
        # get block path
        path2 = block_inst.path
        # get path1 from paths1
        path1 = paths1[block_name]
        # add a list for each block kind to files
        files2[block_name] = []
        # loop around files
        for infilename in files1[block_name]:
            # replace path1 with path2
            outfilename = infilename.replace(path1, path2)
            # add to files
            files2[block_name].append(outfilename)
    # return files
    return files2


def copy_files(files1: Dict[str, List[str]], files2: Dict[str, List[str]]):

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
            outfilename = files2[block_name][it]
            # copy file
            if TEST:
                print('Copying {0} to {1}'.format(infilename, outfilename))
            else:
                shutil.copy(infilename, outfilename)


def update_databases_profile2(params):
    pass


def main():
    """
    Get files from profile 1, get blocks from profile 2
    :return:
    """
    # load parameters
    params = load_params()
    # get a list of files from profile 1 for each block kinds
    files1, paths1 = get_files_profile1(params)
    # get the output files for profile 2 for each block kind
    files2 = get_blocks_profile2(params, files1, paths1)
    # copy files from profile 1 to profile 2 for each block kind
    copy_files(files1, files2)
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
