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
import time
from typing import Any, Dict, Union

import yaml

from apero_request.core import base

# =============================================================================
# Define variables
# =============================================================================
# version, date, author
__VERSION__ = base.__VERSION__
__DATE__ = base.__DATE__
__AUTHOR__ = base.__AUTHOR__
# define relative path to google token files
# TODO: move param1, param2, param3, param4 to local file
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
def read_yaml(yaml_filename: Union[str, None],
              profile_mode: bool = False) -> Dict[str, Any]:
    """
    Read the yaml file and add to settings

    :param yaml_filename: str, yaml file name

    :return: dict, updated settings dictionary
    """
    # deal with yaml_filename being None
    if yaml_filename is None:
        emsg = 'yaml_filename must be set to a valid file'
        raise base.AperoRequestError(emsg)
    # deal with yaml_filename not existing
    if not os.path.exists(yaml_filename):
        emsg = 'yaml_filename {0} does not exist'
        eargs = [yaml_filename]
        raise base.AperoRequestError(emsg.format(*eargs))
    # read the yaml file
    with open(yaml_filename, 'r') as f:
        yaml_data = yaml.load(f, Loader=yaml.FullLoader)
    # add a profiles sub-dictionary
    settings = dict()
    # deal with profile mode
    if profile_mode:
        # add a profiles sub-dictionary
        settings['PROFILES'] = dict()
        # loop around yaml data
        for key, value in yaml_data.items():
            # if key is in settings
            settings['PROFILES'][key] = value
    else:
        # loop around yaml data
        for key, value in yaml_data.items():
            # if key is in settings
            settings[key] = value
    # return settings
    return settings

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


def pull_from_googlesheet(google_sheet, logger=None, **kwargs):
    # push dataframe back to server
    fail_count, error = 0, ''
    while fail_count < 10:
        try:
            return google_sheet.sheet_to_df(**kwargs)
        except Exception as e:
            msg = ('Attempt {0}: Could not pull from googlesheets. '
                   'Trying again in 45 s')
            margs = [fail_count + 1]
            # log the message
            if logger is None:
                print(msg.format(*margs))
            else:
                logger(msg.format(*margs))
            fail_count += 1
            error = e
            time.sleep(45)
    # if we still have not succeeded raise exception here
    emsg = ('Could not pull from google sheets. Tried 10 times. '
            'Error {0}: {1}')
    raise base.AperoRequestError(emsg.format(type(error), str(error)))


def push_to_googlesheet(google_sheet, dataframe, logger=None, **kwargs):
    """
    Push to googlesheets using the gspread class
    Has a while loop to try multiple times before raising an error

    :param google_sheet: gspread instance
    :param dataframe: pandas dataframe to add
    :param logger: logger instance (None means use print)
    :param kwargs: passed to google_sheets.df_to_sheet

    :return:
    """
    # push dataframe back to server
    fail_count, error = 0, ''
    while fail_count < 10:
        try:
            google_sheet.df_to_sheet(dataframe, **kwargs)
            return
        except Exception as e:
            msg = ('Attempt {0}: Could not upload to googlesheets. '
                   'Trying again in 45 s')
            margs = [fail_count + 1]
            # log the message
            if logger is None:
                print(msg.format(*margs))
            else:
                logger(msg.format(*margs))
            fail_count += 1
            error = str(e)
            time.sleep(45)
    # if we still have not succeeded raise exception here
    emsg = ('Could not upload to google sheets. Tried 10 times. '
            'Error {0}: {1}')
    raise base.AperoRequestError(emsg.format(type(error), str(error)))


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
