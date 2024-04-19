#!/usr/bin/env python
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-02-20 at 9:35

@author: cook
"""
import argparse
import os
from typing import Optional

# =============================================================================
# Define variables
# =============================================================================
# PATH ON SERVER
PATH = '/cosmos99/{INSTRUMENT}/apero-data/{PROFILE}/objects/{OBJECT}/'
# SERVER TO USE
SERVER = 'maestria'
# rsync command
RSYNC_CMD = 'rsync --copy-links -avu -e "ssh -oport=5822" {REMOTE} {LOCAL}'


# =============================================================================
# Define classes
# =============================================================================
class AperoProfile:
    def __init__(self, name: str, instrument: str,
                 user: str, host: str, mode: Optional[str] = None):
        self.name = name
        self.user = user
        self.host = host
        self.instrument = instrument
        self.mode = mode

    def get_data(self, objname: str, local_path: str, filestring: str,
                 test: bool = False):
        remote_path = PATH.format(INSTRUMENT=self.instrument,
                                  PROFILE=self.name,
                                  OBJECT=objname)

        remote_path += filestring

        remote_path = f'{self.user}@{self.host}:{remote_path}'

        rsync_cmd = RSYNC_CMD.format(REMOTE=remote_path, LOCAL=local_path)

        print(rsync_cmd)
        if not test:
            os.system(rsync_cmd)


# Define the apero profiles here (after class definition)
APROFS = dict()
APROFS['nirps_he_online'] = AperoProfile(name='nirps_he_online',
                                         instrument='nirps', mode='he',
                                         user='nirps-client', host=SERVER)

APROFS['nirps_ha_online'] = AperoProfile(name='nirps_ha_online',
                                         instrument='nirps', mode='ha',
                                         user='nirps-client', host=SERVER)

APROFS['spirou_online'] = AperoProfile(name='spirou_online',
                                       instrument='spirou',
                                       user='spirou-client', host=SERVER)

APROFS['spirou_offline'] = AperoProfile(name='spirou_offline',
                                        instrument='spirou',
                                        user='spirou-client', host=SERVER)


# =============================================================================
# Define functions
# =============================================================================
def get_args():
    """
    Define the command line arguments

    :return: argparse namespace object containing arguments
    """

    examples = """
    ==================================
    Examples of use
    ==================================
    
    # Get all PROXIMA s1d for NIRPS HE online for fiber A
    #    downloads into current working directory
    >> local_get.py -o PROXIMA -p nirps_he_online -f *_s1d_v_A.fits
    
    # Get all GL699 e2dsff fiber=AB for SPIROU offline
    >> local_get.py -o GL699 -p spirou_offline -f *e2dsff_AB.fits -d /home/cook/downloads/
    
    # Set the username and host used in rsync
    # Prompts for all arguments
    >> local_get.py -u spirou-client -s rali.astro.umontreal.ca  
    """

    parser = argparse.ArgumentParser(description='Get APERO data from remote '
                                                 'server',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=examples)
    # add obs dir
    parser.add_argument('--object', '--obj', '-o', '--objname',
                        type=str, default='None',
                        help='The profile yaml to use')
    parser.add_argument('--profile', '-p', '--name', '-n',
                        type=str, default='None',
                        help='The specific apero profile to use')
    parser.add_argument('--filestring', '--string', '--file', '-f',
                        type=str, default='None',
                        help='Filename search string (with wildcards)')
    parser.add_argument('--path', '--dir', '-d',
                        type=str, default='None',
                        help='Set a path to copy the data to '
                             '(Created if it does not exist)')
    parser.add_argument('--test', action='store_true',
                        help='Run in test mode (do not download)')
    parser.add_argument('--host', '--server', '-s', default='None',
                        help='Define the server name i.e. user@server')
    parser.add_argument('--user', '-u', default='None',
                        help='Define the user name i.e. user@server')
    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return args


def main():
    # get arguments from command line
    params = get_args()
    # deal with object not set (ask)
    if str(params.object).lower() in ['', 'null', 'none']:
        question = 'Please choose an APERO object name:\t'
        params.object = input(question)
    # deal with profile not set (ask)
    if params.profile.lower() in ['', 'null', 'none']:
        question = 'Please choose an APERO profile name:\t'
        params.profile = input(question)
    # check that apero profile is valid
    if params.profile.lower() not in APROFS.keys():
        # need a valid number between 0 and
        uinput = -1
        while uinput not in list(range(len(APROFS))):
            print('Invalid apero profile. Please select from the following')
            for p_it, profile in APROFS.keys():
                print(p_it + 1, profile)
            try:
                uinput = int(input(f'Enter 0 - {len(APROFS) - 1}:\t'))
            except Exception as _:
                continue

        params.profile = list(APROFS.keys())[int(uinput)]

    # deal with no file search string given
    if params.filestring.lower() in ['', 'null', 'none']:
        question = 'Please provide a search string:\t'
        params.filestring = input(question)

    # deal with no path given
    if params.path.lower() in ['', 'null', 'none']:
        params.path = os.path.join(os.getcwd(), params.object)

    if not params.path.endswith(os.sep):
        params.path += os.sep

    # deal with invalid path
    if not os.path.exists(params.path):
        try:
            os.makedirs(str(params.path))
        except Exception as _:
            print('Could not make path')

    # get the apero profile instance
    apero_profile = APROFS[params.profile.lower()]

    if params.host.lower() not in ['', 'null', 'none']:
        apero_profile.host = params.host

    if params.user.lower() not in ['', 'null', 'none']:
        apero_profile.user = params.user

    # get the data
    apero_profile.get_data(objname=params.object,
                           local_path=params.path,
                           filestring=params.filestring,
                           test=params.test)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    main()


# =============================================================================
