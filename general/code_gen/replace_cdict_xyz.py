#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-01-21 at 12:14

@author: cook
"""
import os
import re
from typing import Dict, List, Union
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
PACKAGE_PATH = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/'

INSTRUMENT_PATH = f'{PACKAGE_PATH}/apero/instruments/'

INSTRUMENTS = ['spirou', 'nirps_he', 'nirps_ha']

EXCLUDED_CODES = ['constants.py', 'config.py']

HEADER = '*' * 75

# list of constants with out entries (outside CDict)
MISSING = []
# -----------------------------------------------------------------------------
COLOURS = dict()
COLOURS['BLACK1'] = '\033[90;1m'
COLOURS['RED1'] = '\033[1;91;1m'
COLOURS['GREEN1'] = '\033[92;1m'
COLOURS['YELLOW1'] = '\033[1;93;1m'
COLOURS['BLUE1'] = '\033[94;1m'
COLOURS['MAGENTA1'] = '\033[1;95;1m'
COLOURS['CYAN1'] = '\033[1;96;1m'
COLOURS['WHITE1'] = '\033[97;1m'
COLOURS['BLACK2'] = '\033[1;30m'
COLOURS['RED2'] = '\033[1;31m'
COLOURS['GREEN2'] = '\033[1;32m'
COLOURS['YELLOW2'] = '\033[1;33m'
COLOURS['BLUE2'] = '\033[1;34m'
COLOURS['MAGENTA2'] = '\033[1;35m'
COLOURS['CYAN2'] = '\033[1;36m'
COLOURS['WHITE2'] = '\033[1;37m'
COLOURS['ENDC'] = '\033[0;0m'
COLOURS['BOLD'] = '\033[1m'
COLOURS['UNDERLINE'] = '\033[4m'

# =============================================================================
# Define functions
# =============================================================================
class Colors:
    def __init__(self, theme: Union[str, None] = None):
        """
        Constructor of the colour class (colours based on theme)
        :param theme: str, if set sets the theme ('DARK' or 'LIGHT') defaults
                      to 'DARK'
        """
        # Basic definition of colours to use in log to screen
        self.BLACK1 = COLOURS['BLACK1']
        self.RED1 = COLOURS['RED1']
        self.GREEN1 = COLOURS['GREEN1']
        self.YELLOW1 = COLOURS['YELLOW1']
        self.BLUE1 = COLOURS['BLUE1']
        self.MAGENTA1 = COLOURS['MAGENTA1']
        self.CYAN1 = COLOURS['CYAN1']
        self.WHITE1 = COLOURS['WHITE1']
        self.BLACK2 = COLOURS['BLACK2']
        self.RED2 = COLOURS['RED2']
        self.GREEN2 = COLOURS['GREEN2']
        self.YELLOW2 = COLOURS['YELLOW2']
        self.BLUE2 = COLOURS['BLUE2']
        self.MAGENTA2 = COLOURS['MAGENTA2']
        self.CYAN2 = COLOURS['CYAN2']
        self.WHITE2 = COLOURS['WHITE2']
        self.ENDC = COLOURS['ENDC']
        self.BOLD = COLOURS['BOLD']
        self.UNDERLINE = COLOURS['UNDERLINE']
        # if we have no theme set - set the default
        if theme is None:
            self.theme = 'DARK'
        # if anything else set the theme to them
        else:
            self.theme = theme
        # get inital definitions of themed objects
        self.header = self.MAGENTA1
        self.okblue = self.BLUE1
        self.okgreen = self.GREEN1
        self.ok = self.MAGENTA2
        self.warning = self.YELLOW1
        self.fail = self.RED1
        self.debug = self.BLACK1
        # define the end of string code (block to reset)
        self.endc = self.ENDC
        # define the bold string code
        self.bold = self.BOLD
        # define the underline string code
        self.underline = self.UNDERLINE
        # update all others via theme
        self.update_theme()

    def __getstate__(self) -> dict:
        """
        For when we have to pickle the class
        :return:
        """
        # set state to __dict__
        state = dict(self.__dict__)
        # return dictionary state (for pickle)
        return state

    def __setstate__(self, state):
        """
        For when we have to unpickle the class

        :param state: dictionary from pickle
        :return:
        """
        # update dict with state
        self.__dict__.update(state)

    def __str__(self) -> str:
        """
        Return string represenation of Const class
        :return:
        """
        return 'Colors[{0}]'.format(self.theme)

    def update_theme(self, theme: Union[str, None] = None):
        """
        Update themed object names
            header/okblue/okgreen/ok/warning/fail/debug
        based on theme

        :param theme: str, if set sets the theme ('DARK' or 'LIGHT') defaults
                      to 'DARK'
        :return:
        """
        # if we have no theme set - set the default
        if theme is not None:
            self.theme = theme
        # set the dark colours
        if self.theme == 'DARK':
            self.header = self.MAGENTA1
            self.okblue = self.BLUE1
            self.okgreen = self.GREEN1
            self.ok = self.MAGENTA2
            self.warning = self.YELLOW1
            self.fail = self.RED1
            self.debug = self.BLACK1
        # set the light colours
        else:
            self.header = self.MAGENTA2
            self.okblue = self.MAGENTA2
            self.okgreen = self.BLACK2
            self.ok = self.MAGENTA2
            self.warning = self.BLUE2
            self.fail = self.RED2
            self.debug = self.GREEN2

    def print(self, message: str, colour: str,
              highlight_words: List[str] = None,
              highlight_colour: str = 'green') -> str:
        """
        A basic coloured print mesage
        If colour is incorrect does nothing

        :param message: str, the message to print
        :param colour: str, the colour to print, colour must be one of the
                       following: b, r, h, y, m, k

        :return: a coloured string ready to be printed to stdout
        """
        if colour in ['b', 'blue']:
            start = self.BLUE1
        elif colour in ['r', 'red']:
            start = self.RED1
        elif colour in ['g', 'green']:
            start = self.GREEN1
        elif colour in ['y', 'yellow']:
            start = self.YELLOW1
        elif colour in ['m', 'magenta']:
            start = self.MAGENTA1
        elif colour in ['k', 'black', 'grey']:
            start = self.BLACK1
        else:
            start = self.endc


        if highlight_words is not None:

            if highlight_colour in ['b', 'blue']:
                hstart = self.BLUE1
            elif highlight_colour in ['r', 'red']:
                hstart = self.RED1
            elif highlight_colour in ['g', 'green']:
                hstart = self.GREEN1
            elif highlight_colour in ['y', 'yellow']:
                hstart = self.YELLOW1
            elif highlight_colour in ['m', 'magenta']:
                hstart = self.MAGENTA1
            elif highlight_colour in ['k', 'black', 'grey']:
                hstart = self.BLACK1
            else:
                hstart = self.endc
            # for all highlighed words replace word with:
            #       color + word + end + original colour
            for highlight_word in highlight_words:
                replacement = hstart + highlight_word + self.endc + start
                message = message.replace(highlight_word, replacement)

        # return colour message
        return start + message + self.endc

    def cprint(self, message, colour='green',
              highlight_words: List[str] = None,
              highlight_colour: str = 'green'):
        print(self.print(message, colour, highlight_words, highlight_colour))


CC = Colors()


def find_default_constants(file_path: str, method='add') -> Dict[str, str]:
    # Regular expression to capture the first argument of CDict.add
    pattern = (r'CDict\.add\(\s*[\'"]([^\'"]+)[\'"]'
               r'(?:[^)]*?group=\s*([^\s,]+))?')
    # Open and read the file
    with open(file_path, 'r') as file:
        content = file.read()
    # Find all matches using the regex
    matches = re.findall(pattern, content)

    from apero.instruments.default.consants import CDict


    constants_list = dict()

    for match in matches:
        constants_list[match[0]] = match[1] if match[1] else None
    # return constant list
    return constants_list


def find_constants_without_group(file_path: str, constants_list: Dict[str, str],
                                 method='set') -> Dict[str, str]:
    # Regular expression to capture the first argument of CDict.add
    # and check for the presence of the "group" argument
    pattern = r'CDict\.' + method + r'\(\s*[\'"]([^\'"]+)[\'"].*?(group=)'
    # Open and read the file
    with open(file_path, 'r') as file:
        content = file.read()
    # Find all matches using the regex
    matches = re.findall(pattern, content)
    # storage to return
    matches_with_group = []
    # loop around matches and only keep those without a group
    for match in matches:
        if match[0] in constants_list:
            matches_with_group.append(match[0])
    # keep a list of those in constants_list but not in matches_with_group
    matches_without_group = dict()
    for constant_name in constants_list.keys():
        if constant_name not in matches_with_group:
            matches_without_group[constant_name] = constants_list[constant_name]
    # return matches without a group
    return matches_without_group


def get_all_python_files(file_path: str):
    # list to return
    python_files = []
    # loop file_path
    for root, dir, files in os.walk(file_path):
        for filename in files:
            # skip filename if its in excluded codes
            if filename in EXCLUDED_CODES:
                continue
            # only consider python files
            if filename.endswith('.py'):
                python_files.append(os.path.join(root, filename))
    # return all python files
    return python_files


def read_all_python_files(python_files: List[str]) -> Dict[str, List[str]]:
    # python dictionary to return
    python_dict = dict()
    # loop around python files
    for python_file in python_files:
        # open file and read
        with open(python_file) as pfile:
            python_dict[str(python_file)] = pfile.readlines()

    return python_dict


def update_constant(constant_name, all_python_lines, group=None):

    # print the variable name
    CC.cprint(HEADER, colour='magenta')
    CC.cprint(f'Processing {constant_name}', colour='magenta')
    CC.cprint(HEADER, colour='magenta')


    CC.cprint('Finding instances')
    # storage python files
    python_files = dict()
    # find all instances of string in all python files
    for python_file in tqdm(all_python_lines.keys()):

        for l_it, line in enumerate(all_python_lines[python_file]):

            if f'\'{constant_name}\'' in line:
                python_files[python_file] = l_it

    # if none found we should keep a list
    if len(python_files) == 0:
        global MISSING
        MISSING.append(constant_name)
        CC.cprint('\tConstant not found outside definition. Skipping',
                  colour='yellow')
        return
    # -------------------------------------------------------------------------
    # print out these entries (removing the package path from python file
    for python_file in python_files:

        # truncated file name
        string_filename = python_file.replace(PACKAGE_PATH, '')

        CC.cprint(string_filename, 'blue')

        # get lines
        lines = all_python_lines[python_file]

        line_number = python_files[python_file]


        line_range = [max([0, line_number-3]),
                      min([line_number+3, len(lines)])]

        for line_it in range(*line_range):

            fmt_line = lines[line_it].replace('\n', '')
            CC.cprint(f'{str(line_it):5s}| {fmt_line}',
                      colour='green', highlight_colour='red',
                      highlight_words=[constant_name])
        CC.cprint('\n')
    # -------------------------------------------------------------------------






# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # step 1: identify constants in default constants file
    constants_list = find_default_constants(INSTRUMENT_PATH +
                                            'default/constants.py')

    # -------------------------------------------------------------------------
    # step 2: find matching constants without a group in the instrument files
    # -------------------------------------------------------------------------
    valid_constants_list = dict()
    # loop around instruments
    for instrument in INSTRUMENTS:
        _path = INSTRUMENT_PATH + f'{instrument}/constants.py'
        # find groups without constants
        constants_without_group = find_constants_without_group(_path,
                                                               constants_list)
        valid_constants_list.update(constants_without_group)
    # -------------------------------------------------------------------------
    # step 3: read all python files
    # -------------------------------------------------------------------------
    # get a list of all python files in the package
    all_python_files = get_all_python_files(PACKAGE_PATH)
    # read all python files and stora the lines in memory in a dictionary
    all_python_lines = read_all_python_files(all_python_files)
    # -------------------------------------------------------------------------
    # step 4: ask user for confirmation and edit files
    # -------------------------------------------------------------------------
    # loop around constants, display the variable, ask for the new name, and
    # then confirm changes, then write changes to files
    for constant_name in valid_constants_list.keys():

        update_constant(constant_name, all_python_lines,
                        group=valid_constants_list[constant_name])





# =============================================================================
# End of code
# =============================================================================
