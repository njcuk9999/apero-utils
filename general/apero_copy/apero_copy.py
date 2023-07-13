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
from typing import Any, Dict, Tuple, Union

import yaml

# =============================================================================
# Define variables
# =============================================================================
__version__ = '0.0.1'

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
class AperoCopyError(Exception):
    """
    Error class for apero_raw_tests
    """
    pass



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


def load_params():
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



def main():

    # load parameters
    params = load_params()

    # get a list of files from profile 1 for each block kinds

    # get the block kinds from profile 2

    # copy files from profile 1 to profile 2

    # update databases for profile 2

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
