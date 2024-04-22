#!/usr/bin/env python
"""
Locally get APERO files from the /cosmos99 drive using rsync

for any spirou or nirps profile

Requires the manual_trigger to be run to symlink to objects directory

Not all files are copied

Created on 2024-04-19 18:19

@author: cook
"""
import argparse
import os

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'local_apero_get.py'
__VERSION__ = '0.0.1'
__DATE__ = '2024-04-19'
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
                 user: str, host: str):
        """
        Construct the AperoProfile class
        :param name: str, the profile name of this APERO profile
        :param instrument: str, the instrument (spirou or nirps)
        :param user: str, the default username used for ssh (i.e. ssh user@host)
        :param host: str, the default host used for ssh  (i.e. ssh user@host)
        """
        self.name = name
        self.user = user
        self.host = host
        self.instrument = instrument

    def get_data(self, objname: str, local_path: str, filestring: str,
                 multi: bool = False, test: bool = False):
        """
        Create the remote_path and rsync command and run it

        :param objname: str, the object name (as in APERO astrometric database)
        :param local_path: str, the local path to save files to
        :param filestring: str, the file search string to use
        :param test: bool, whether we are in test mode, if True does not run
                     the rsync command
        :return: None
        """
        # deal with multi and localpath
        if multi:
            lpath = os.path.join(local_path, objname)
            # deal with local path not existing
            if not os.path.exists(lpath):
                os.makedirs(lpath)
        else:
            lpath = str(local_path)

        # construct the remove path
        remote_path = PATH.format(INSTRUMENT=self.instrument,
                                  PROFILE=self.name,
                                  OBJECT=objname)
        # add the search string to the remote path
        remote_path += filestring
        # the remote path must have format "user@host:path"
        remote_path = f'{self.user}@{self.host}:{remote_path}'
        # construct the rsync command
        rsync_cmd = RSYNC_CMD.format(REMOTE=remote_path, LOCAL=lpath)
        # print the rsync command
        print(rsync_cmd)
        # if not in test mode run the rsync command
        if not test:
            os.system(rsync_cmd)


# Define the apero profiles here (after class definition)
APROFS = dict()
APROFS['nirps_he_online'] = AperoProfile(name='nirps_he_online',
                                         instrument='nirps',
                                         user='nirps-client', host=SERVER)

APROFS['nirps_ha_online'] = AperoProfile(name='nirps_ha_online',
                                         instrument='nirps',
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
    description = """
    Locally get APERO files from the /cosmos99 drive using rsync
    
    for any spirou or nirps profile
    
    Requires the manual_trigger to be run to symlink to objects directory
    
    Not all files are copied
    """

    examples = """
    ==================================
    Examples of use
    ==================================
    
    # Get all PROXIMA s1d for NIRPS HE online for fiber A
    #    downloads into current working directory
    >> local_apero_get.py -o PROXIMA -p nirps_he_online -f *_s1d_v_A.fits
    
    # Get all GL699 e2dsff fiber=AB for SPIROU offline
    >> local_apero_get.py -o GL699 -p spirou_offline -f *e2dsff_AB.fits -d /home/cook/downloads/
    
    # Set the username and host used in rsync
    # Prompts for all arguments
    >> local_apero_get.py -u spirou-client -s rali.astro.umontreal.ca  
    """
    # create the parser class
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=examples)
    # add obs dir
    parser.add_argument('--object', '--obj', '-o', '--objname',
                        type=str, default='None',
                        help='The object name to search for (as in APERO '
                             'database). Note you can NOT use an alias here.')
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
        question = '\nPlease choose an APERO object name:\t'
        params.object = input(question)

    if ',' not in str(params.object):
        objnames = [params.object]
        multi = False
    else:
        objnames = params.object.split(',')
        multi = True

    # deal with profile not set (ask)
    if params.profile.lower() in ['', 'null', 'none']:
        question = '\nPlease choose an APERO profile name:\t'
        params.profile = input(question)
    # check that apero profile is valid
    if params.profile.lower() not in APROFS.keys():
        # need a valid number between 0 and
        uinput = -1
        while uinput not in list(range(len(APROFS))):
            print('\nInvalid apero profile. Please select from the following')
            for p_it, profile in enumerate(APROFS.keys()):
                print(f'{p_it + 1}. {profile}')
            # noinspection PyBroadException
            try:
                uinput = int(input(f'Enter 1 - {len(APROFS)}:\t'))
            except Exception as _:
                continue
        # set profile based on user input
        params.profile = list(APROFS.keys())[int(uinput)-1]

    # deal with no file search string given
    if params.filestring.lower() in ['', 'null', 'none']:
        question = '\nPlease provide a search string:\t'
        params.filestring = input(question)

    # deal with no path given
    if params.path.lower() in ['', 'null', 'none']:
        params.path = os.path.join(os.getcwd(), params.object)

    if not params.path.endswith(os.sep):
        params.path += os.sep

    # deal with invalid path
    if not os.path.exists(params.path):
        # noinspection PyBroadException
        try:
            os.makedirs(str(params.path))
        except Exception as _:
            print('\nPath (--path --dir -d) invalid. Could not make path.')
            return

    # get the apero profile instance
    apero_profile = APROFS[params.profile.lower()]
    # deal with a host override
    if params.host.lower() not in ['', 'null', 'none']:
        apero_profile.host = params.host
    # deal with a username override
    if params.user.lower() not in ['', 'null', 'none']:
        apero_profile.user = params.user

    # loop around object names
    for objname in objnames:
        # get the data
        apero_profile.get_data(objname=objname,
                               local_path=params.path,
                               filestring=params.filestring,
                               multi=multi,
                               test=params.test)


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # run the main function to process args and create/run the rsync
    main()

# =============================================================================
