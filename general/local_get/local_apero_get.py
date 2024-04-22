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
import string

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'local_apero_get.py'
__VERSION__ = '0.0.2'
__DATE__ = '2024-04-22'
# PATH ON SERVER
PATH = '/cosmos99/{INSTRUMENT}/apero-data/{PROFILE}/objects/{OBJECT}/'
# SERVER TO USE
SERVER = 'maestria'
# rsync command
RSYNC_CMD = 'rsync --copy-links -avu -e "ssh -oport=5822" {REMOTE} {LOCAL}'
# define bad characters for objects (alpha numeric + "_")
BAD_OBJ_CHARS = [' '] + list(string.punctuation.replace('_', ''))
# Define colours for printing
COLOURS = dict()
COLOURS['BLACK'] = '\033[90;1m'
COLOURS['RED'] = '\033[1;91;1m'
COLOURS['GREEN'] = '\033[92;1m'
COLOURS['YELLOW'] = '\033[1;93;1m'
COLOURS['BLUE'] = '\033[94;1m'
COLOURS['MAGENTA'] = '\033[1;95;1m'
COLOURS['CYAN'] = '\033[1;96;1m'
COLOURS['WHITE'] = '\033[97;1m'
COLOURS['ENDC'] = '\033[0;0m'
# logo
LOGO = ["\t\033[1;91;1m  █████\033[1;37m╗\033[1;91;1m ██████\033[1;37m╗"
        "\033[1;91;1m ███████\033[1;37m╗\033[1;91;1m██████\033[1;37m╗"
        "\033[1;91;1m  ██████\033[1;37m╗\033[1;91;1m  ",

        "\t ██\033[1;37m╔\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m██\033[1;37m╗\033[1;91;1m██\033[1;37m╔"
        "\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m██"
        "\033[1;37m╗\033[1;91;1m██\033[1;37m╔\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m"
        "\033[1;37m═\033[1;91;1m\033[1;37m╝\033[1;91;1m██\033[1;37m╔"
        "\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m██"
        "\033[1;37m╗\033[1;91;1m██\033[1;37m╔\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m██"
        "\033[1;37m╗\033[1;91;1m ",

        "\t ███████\033[1;37m║\033[1;91;1m██████\033[1;37m╔\033[1;91;1m"
        "\033[1;37m╝\033[1;91;1m█████\033[1;37m╗\033[1;91;1m  ██████"
        "\033[1;37m╔\033[1;91;1m\033[1;37m╝\033[1;91;1m██\033[1;37m║"
        "\033[1;91;1m   ██\033[1;37m║\033[1;91;1m ",

        "\t ██\033[1;37m╔\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m██\033[1;37m║\033[1;91;1m██\033[1;37m╔\033[1;91;1m"
        "\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m\033[1;37m╝\033[1;91;1m ██\033[1;37m╔\033[1;91;1m"
        "\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m╝"
        "\033[1;91;1m  ██\033[1;37m╔\033[1;91;1m\033[1;37m═\033[1;91;1m"
        "\033[1;37m═\033[1;91;1m██\033[1;37m╗\033[1;91;1m██\033[1;37m║"
        "\033[1;91;1m   ██\033[1;37m║\033[1;91;1m ",

        "\t ██\033[1;37m║\033[1;91;1m  ██\033[1;37m║\033[1;91;1m██"
        "\033[1;37m║\033[1;91;1m     ███████\033[1;37m╗\033[1;91;1m██"
        "\033[1;37m║\033[1;91;1m  ██\033[1;37m║\033[1;91;1m\033[1;37m╚"
        "\033[1;91;1m██████\033[1;37m╔\033[1;91;1m\033[1;37m╝"
        "\033[1;91;1m ",

        "\t \033[1;37m╚\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m╝"
        "\033[1;91;1m  \033[1;37m╚\033[1;91;1m\033[1;37m═\033[1;91;1m"
        "\033[1;37m╝\033[1;91;1m\033[1;37m╚\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m\033[1;37m╝\033[1;91;1m     \033[1;37m╚"
        "\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m"
        "\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m╝\033[1;91;1m"
        "\033[1;37m╚\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m╝"
        "\033[1;91;1m  \033[1;37m╚\033[1;91;1m\033[1;37m═\033[1;91;1m"
        "\033[1;37m╝\033[1;91;1m \033[1;37m╚\033[1;91;1m\033[1;37m═"
        "\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m"
        "\033[1;37m═\033[1;91;1m\033[1;37m═\033[1;91;1m\033[1;37m╝"
        "\033[1;91;1m  \033[0;0m"]

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
        cprint()
        cprint(rsync_cmd, colour='blue')
        cprint()
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
def cprint(text: str = '', colour: str = 'green', kind='print'):
    """
    Print text in a given color
    :param text: str, the text to print
    :param color: str, the color to print the text in
    :return: None
    """
    if colour.upper() in COLOURS:
        print_colour = COLOURS[colour.upper()]
        end_colour = COLOURS['ENDC']
    else:
        print_colour = ''
        end_colour = ''

    if kind == 'input':
        return input(f'{print_colour}{text}{end_colour}')
    else:
        print(f'{print_colour}{text}{end_colour}')


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
    # print arguments used
    cprint('Arguments used:', colour='blue')
    # flag for having arguments
    has_args = False
    # loop around arguments
    for arg in vars(args):
        # get value
        value = getattr(args, arg)
        # skip test if not in test mode
        if arg == 'test' and not value:
            continue
        if str(value).lower() not in ['', 'null', 'none']:
            cprint('\t- {0}: {1}'.format(arg, getattr(args, arg)),
                   colour='blue')
            has_args = True
    # deal with no arguments given
    if not has_args:
        cprint('\t- No arguments given', colour='blue')
    header_line()
    # return arguments
    return args


def header_line(symbol='*', colour='green'):
    """
    Print the header
    :return: None
    """
    cprint(symbol * 60, colour=colour)


def clean_object(rawobjname: str) -> str:
    """
    Clean a 'rawobjname' to allow it to be consistent

    :param rawobjname: str, the raw object name to clean

    :return: str, the cleaned object name
    """
    # if raw object name contains null text - return Null string
    if str(rawobjname).lower() in ['null', 'none', '']:
        return 'Null'
    # strip spaces off raw object
    objectname = rawobjname.strip()
    # replace + and - with "p" and "m"
    objectname = objectname.replace('+', 'p')
    objectname = objectname.replace('-', 'm')
    # now remove bad characters
    for bad_char in BAD_OBJ_CHARS:
        objectname = objectname.replace(bad_char, '_')
    objectname = objectname.upper()
    # deal with multiple underscores in a row
    while '__' in objectname:
        objectname = objectname.replace('__', '_')
    # strip leading / trailing '_'
    objectname = objectname.strip('_')
    # return cleaned object name
    return objectname


def splash():
    """
    Print the splash screen
    :return: None
    """
    # print end message
    header_line()
    for line in LOGO:
        cprint(line, colour='red')
    header_line()
    cprint(f'{__NAME__}\t\tVersion: {__VERSION__}', colour='red')
    header_line()


def main():
    """
    Main code

    :return:
    """
    # print splash screen
    splash()
    # get arguments from command line
    params = get_args()

    # deal with object not set (ask)
    if str(params.object).lower() in ['', 'null', 'none']:
        question = '\nPlease choose an APERO object name:\t'
        params.object = cprint(question, colour='magenta', kind='input')

    if ',' not in str(params.object):
        objnames = [params.object]
        multi = False
    else:
        objnames = params.object.split(',')
        multi = True

    # deal with profile not set (ask)
    if params.profile.lower() in ['', 'null', 'none']:
        question = '\nPlease choose an APERO profile name:\t'
        params.profile = cprint(question, colour='magenta', kind='input')
    # check that apero profile is valid
    if params.profile.lower() not in APROFS.keys():
        # need a valid number between 0 and
        uinput = -1
        while uinput not in list(range(len(APROFS))):
            cprint('\nInvalid apero profile. Please select from the following',
                   colour='magenta')
            for p_it, profile in enumerate(APROFS.keys()):
                cprint(f'{p_it + 1}. {profile}', colour='magenta')
            # noinspection PyBroadException
            try:
                question = f'Enter 1 - {len(APROFS)}:\t'
                uinput = int(cprint(question, colour='magenta', kind='input'))
            except Exception as _:
                continue
        # set profile based on user input
        params.profile = list(APROFS.keys())[int(uinput) - 1]

    # deal with no file search string given
    if params.filestring.lower() in ['', 'null', 'none']:
        question = '\nPlease provide a search string:\t'
        params.filestring = cprint(question, colour='magenta', kind='input')

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
            cprint('\nPath (--path --dir -d) invalid. Could not make path.',
                   colour='red')
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
        # clean object name
        objname = clean_object(objname)
        # print progress
        cprint('\n', colour='green')
        header_line(symbol='=')
        cprint(f'Getting data for object: {objname}', colour='green')
        header_line(symbol='=')
        # get the data
        apero_profile.get_data(objname=objname,
                               local_path=params.path,
                               filestring=params.filestring,
                               multi=multi,
                               test=params.test)

    # print end message
    header_line()
    cprint('Code ended successfully', colour='green')
    header_line()


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # run the main function to process args and create/run the rsync
    main()

# =============================================================================
