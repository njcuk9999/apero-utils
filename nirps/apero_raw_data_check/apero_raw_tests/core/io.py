#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:51

@author: cook
"""
import os
import platform
from typing import Any, Dict, Union

import yaml

from apero_raw_tests.core import base

# =============================================================================
# Define variables
# =============================================================================
__VERSION__ = base.__VERSION__
__AUTHOR__ = base.__AUTHOR__
# define relative path to google token files
PARAM1 = ('241559402076-vbo2eu8sl64ehur7'
          'n1qhqb0q9pfb5hei.apps.googleusercontent.com')
PARAM2 = ('apero-data-manag-', '1602517149890')
PARAM3 = ''.join(base.GSPARAM)
PARAM4 = ('1//0dBWyhNqcGHgdCgYIARAAGA0SNwF-L9IrhXoPCjWJtD4f0EDxA',
          'gFX75Q-f5TOfO1VQNFgSFQ_89IW7trN3B4I0UYvvbVfrGRXZZg')
PATH1 = 'gspread_pandas/google_secret.json'
PATH2 = 'gspread_pandas/creds/default'
TEXT1 = ('{{"installed":{{"client_id":"{0}","project_id":"{1}","auth_uri":'
         '"https://accounts.google.com/o/oauth2/auth","token_uri":'
         '"https://oauth2.googleapis.com/token","auth_provider_x509_cert'
         '_url":"https://www.googleapis.com/oauth2/v1/certs","client_'
         'secret":"{2}","redirect_uris":["urn:ietf:wg:oauth:2.0:oob",'
         '"http://localhost"]}}}}')
TEXT2 = ('{{"refresh_token": "{0}", "token_uri": "https://oauth2.googleap'
         'is.com/token", "client_id": "{1}", "client_secret": "{2}", '
         '"scopes": ["openid", "https://www.googleapis.com/auth/drive", '
         '"https://www.googleapis.com/auth/userinfo.email", '
         '"https://www.googleapis.com/auth/spreadsheets"]}}')

# =============================================================================
# Define functions
# =============================================================================
def read_yaml(yaml_filename: Union[str, None]) -> Dict[str, Any]:
    """
    Read the yaml file and add to settings

    :param yaml_filename: str, yaml file name
    :param settings: dict, settings dictionary

    :return: dict, updated settings dictionary
    """
    # deal with yaml_filename being None
    if yaml_filename is None:
        emsg = 'yaml_filename must be set to a valid file'
        raise base.AperoRawTestsError(emsg)
    # deal with yaml_filename not existing
    if not os.path.exists(yaml_filename):
        emsg = 'yaml_filename {0} does not exist'
        eargs = [yaml_filename]
        raise base.AperoRawTestsError(emsg.format(*eargs))
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


def get_obs_dirs(params: Dict[str, Any]):
    """
    Get all raw observation directories from params

    :param params:
    :return:
    """
    # get the raw path
    raw_path = params['raw dir']
    # test that raw directory exists
    if not os.path.exists(raw_path):
        emsg = 'Raw directory {0} does not exist'
        eargs = [raw_path]
        raise base.AperoRawTestsError(emsg.format(*eargs))
    # get a list of all fits files in raw path
    valid_files = []
    for root, dirs, files in os.walk(raw_path):
        for filename in files:
            if filename.endswith('.fits'):
                valid_files.append(os.path.join(root, filename))
    # now we have a list of files we need to get the unique directories
    obs_dirs = []
    for filename in valid_files:
        obs_dir = os.path.dirname(filename)
        # we need to remove raw_path from the obs_dir
        obs_dir = obs_dir.replace(raw_path, '')
        # strip any leading separators
        while obs_dir.startswith(os.sep):
            obs_dir = obs_dir[len(os.sep):]
        # add to the list of observation directories
        if obs_dir not in obs_dirs:
            obs_dirs.append(obs_dir)
    # sort alphabetically
    obs_dirs.sort()
    # return the list of observation directories
    return obs_dirs


def gsp_setup():
    # make sure token is in correct directory
    if platform.system() == 'Windows':
        outpath = os.path.join(os.path.expanduser('~'), 'AppData', 'Roaming')
    else:
        outpath = os.path.join(os.path.expanduser('~'), '.config/')
    # make sure .config exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    # construct paths
    path1 = os.path.join(outpath, PATH1)
    path2 = os.path.join(outpath, PATH2)
    # make sure paths exist
    if not os.path.exists(os.path.dirname(path1)):
        os.makedirs(os.path.dirname(path1))
    if not os.path.exists(os.path.dirname(path2)):
        os.makedirs(os.path.dirname(path2))
    # make file
    with open(path1, 'w') as file1:
        file1.write(TEXT1.format(PARAM1, ''.join(PARAM2), PARAM3))
    with open(path2, 'w') as file2:
        file2.write(TEXT2.format(''.join(PARAM4), PARAM1, PARAM3))



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # print 'Hello World!'
    print("Hello World!")

# =============================================================================
# End of code
# =============================================================================
