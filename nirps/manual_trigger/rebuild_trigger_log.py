#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-03-08 at 11:59

@author: cook
"""
import glob
import os
import platform
from typing import Any, Dict, Union

import gspread_pandas as gspd
import numpy as np
import yaml

import manual_trigger

# =============================================================================
# Define variables
# =============================================================================
# -----------------------------------------------------------------------------
# define messages
MANUAL_START = manual_trigger.MANUAL_START
MANUAL_END = manual_trigger.MANUAL_END
APERO_START = manual_trigger.APERO_START
APERO_ERR = manual_trigger.APERO_ERR
APERO_END = manual_trigger.APERO_END
ARI_START = manual_trigger.ARI_START
ARI_END = manual_trigger.ARI_END
# Requests sheet ID
SHEET_ID = '1zvU_XFA1ZOJE111qZKiav7v6ptYveWpDkMHjdhiN06M'
# TODO: move param1, param2, param3, param4 to local file
PARAM1 = ('241559402076-vbo2eu8sl64ehur7'
          'n1qhqb0q9pfb5hei.apps.googleusercontent.com')
PARAM2 = ('apero-data-manag-', '1602517149890')
PARAM3 = ''.join(manual_trigger.GSPARAM)
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


def read_yaml(yaml_filename: Union[str, None]) -> Dict[str, Any]:
    """
    Read the yaml file and add to settings

    :param yaml_filename: str, yaml file name

    :return: dict, updated settings dictionary
    """
    # deal with yaml_filename being None
    if yaml_filename is None:
        emsg = 'yaml_filename must be set to a valid file'
        raise ValueError(emsg)
    # deal with yaml_filename not existing
    if not os.path.exists(yaml_filename):
        emsg = 'yaml_filename {0} does not exist'
        eargs = [yaml_filename]
        raise ValueError(emsg.format(*eargs))
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


def load_sheet(sheet_name: str):
    # add gspread directory and auth files
    gsp_setup()
    # load google sheet instance
    google_sheet = gspd.spread.Spread(SHEET_ID)
    # convert google sheet to pandas dataframe
    return google_sheet.sheet_to_df(index=0, sheet=sheet_name)

def main():
    # construct log path
    homedir = os.path.expanduser('~')
    # construct log path
    logpath = os.path.join(homedir, '.apero', 'manual_trigger')
    # load yaml files
    yaml_files = glob.glob('*.yaml')
    # loop around yaml files
    for yaml_file in yaml_files:
        # load yaml dictionary
        yaml_dict = read_yaml(yaml_file)
        # get the profiles
        profiles = list(yaml_dict.keys())
        # storage of log lines
        lines = dict()
        # loop around profiles
        for profile in profiles:

            if 'check' not in yaml_dict[profile]:
                continue
            checks = yaml_dict[profile]['check']
            if 'run_check' not in checks:
                continue
            if not checks['run_check']:
                continue
            if 'red sheet name' not in checks:
                continue
            # construct log filename
            logfile = os.path.join(logpath, yaml_file.replace('.yaml', '.log'))
            # deal with a different log file
            if (logfile, profile) not in lines:
                lines[(logfile, profile)] = []

            red_sheet = checks['red sheet name']
            # load raw sheet
            red_df = load_sheet(red_sheet)
            # loop around observation directories
            for row in range(len(red_df)):
                # get obs dir for this row
                obsdir = red_df['obsdir'].iloc[row]
                # get date submitted for this row
                date_mod = f'{obsdir} 00:00:00'
                # print progress
                margs = [profile, obsdir, row+1, len(red_df), date_mod]
                print('Processing {0} {1} [{2}/{3}]: {4}'.format(*margs))
                # load trigger log
                logger = manual_trigger.TriggerLog(logfile, profile, [obsdir])
                # go through all keys and add them
                for logkind in manual_trigger.MESSAGES:
                    if logkind not in red_df.columns:
                        continue
                    if red_df[logkind].iloc[row] == 'TRUE':
                        line = logger.write(logkind, message='None',
                                            timestr=date_mod, return_line=True)
                        lines[(logfile, profile)].append(line)
        if len(lines) == 0:
            continue
        # loop around line keys
        for key in lines:
            logfile, profile = key
            # load trigger log
            logger = manual_trigger.TriggerLog(logfile, profile, None)
            # sort lines
            sorted_lines = list(np.sort(lines[key]))[::-1]
            # print progress
            print('Writing log file {0}'.format(logger.filename))
            # write these lines
            logger.write_all(sorted_lines)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    main()

# =============================================================================
# end of code
# =============================================================================