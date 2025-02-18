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
import time
from typing import Dict, List, Tuple, Union
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
PACKAGE_PATH = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/'

INSTRUMENT_PATH = f'{PACKAGE_PATH}/apero/instruments/'

INSTRUMENTS = ['spirou', 'nirps_he', 'nirps_ha']

EXCLUDED_CODES = ['constants.py', 'config.py']

HEADER = '*' * 120

# list of constants with out entries (outside CDict)
MISSING = []
NO_GROUP = []
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
        print(self.print(message, colour, highlight_words, highlight_colour),
              flush=True)


CC = Colors()


class StopLoop(Exception):
    pass


def find_default_constants(file_path: str, method='add') -> Dict[str, str]:

    from apero.instruments.default.constants import CDict

    constants_list = dict()

    for constant_name in CDict.storage:
        # get group
        group = CDict.storage[constant_name].group
        # deal with no group
        if group is None:
            constants_list[constant_name] = None
            continue
        # get name without group
        if group + '.' in constant_name:
            name = constant_name.split(group + '.')[-1]
        else:
            name = constant_name
        # push into list
        constants_list[name] = group
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
    excluded_files = []
    # loop file_path
    for root, dir, files in os.walk(file_path):
        for filename in files:
            # skip filename if its in excluded codes
            if filename in EXCLUDED_CODES:
                excluded_files.append(os.path.join(root, filename))
                continue
            # only consider python files
            if filename.endswith('.py'):
                python_files.append(os.path.join(root, filename))
    # return all python files
    return python_files, excluded_files


def read_all_python_files(python_files: List[str]) -> Dict[str, List[str]]:
    # python dictionary to return
    python_dict = dict()
    # loop around python files
    for python_file in python_files:
        # open file and read
        with open(python_file) as pfile:
            lines = pfile.readlines()
            vlines = [line.replace('\n', '') for line in lines]
            python_dict[str(python_file)] = vlines

    return python_dict


def add_group_to_constant(old_name: str, new_name: str, content: List[str],
                          file_path: str) -> Tuple[List[str], int, int]:

    # Regular expression to match the specific constant
    pattern = rf"(CDict\.\w+\(\s*['\"]{old_name}['\"],.*?)(\))"
    # Join content into a single string for regex processing
    content_string = '\n'.join(content)
    # truncated file name
    string_filename = file_path.replace(PACKAGE_PATH, '')
    # do a basic match to see if we need to process
    matches = re.findall(pattern, content_string, flags=re.DOTALL)
    # don't both if not found
    if len(matches) == 0:
        CC.cprint(f'\n\nRegex failed for {old_name} in {string_filename}', colour='magenta')
        return content, 0, 0
    else:
        CC.cprint(f'\n\nProcessing {old_name} in {string_filename}', colour='magenta')

    state = dict(match_found=False, start_pos=0, end_pos=0)

    # Function to handle replacements
    def replace_callback(match):
        # Set the flag when a match is found
        state['match_found'] = True
        # The full matched call
        full_match = match.group(1)
        # start and end position
        state['start_pos'] = match.start()
        state['end_pos'] = match.end()
        # Check if "group=" is already in the arguments
        if "group=" in full_match:
            # Only rename the constant if needed
            return full_match.replace(old_name, new_name) + match.group(2)

        state['end_pos'] += len(", group=cgroup" + match.group(2))

        # Add the `group=cgroup` and rename the constant
        return (full_match.replace(old_name, new_name) +
                ", group=cgroup" + match.group(2))

    # Replace the specific constant and conditionally add `group=cgroup`
    updated_content_string = re.sub(pattern, replace_callback, content_string, flags=re.DOTALL)

    if not state['match_found']:
        return content, 0, 0
    # split content back into lines
    updated_content = updated_content_string.splitlines()
    # get the line positions of start and end
    start_line, end_line = string_pos_to_line_pos(updated_content,
                                                  start=int(state['start_pos']),
                                                  end=int(state['end_pos']))
    # return these values
    return updated_content, start_line, end_line



def update_constant(constant_name, all_python_lines,
                    const_python_lines, updated_lines, group=None):

    if group is None:
        CC.cprint('\tNo group found for constant. Skipping', colour='yellow')
        global NO_GROUP
        NO_GROUP.append(constant_name)
        return all_python_lines, const_python_lines, updated_lines
    # print progress
    CC.cprint('Finding instances...', colour='magenta')
    # storage python files
    python_files = dict()
    group_python_files = dict()
    # find all instances of string in all python files
    for python_file in all_python_lines.keys():
        for l_it, line in enumerate(all_python_lines[python_file]):

            if f'\'{constant_name}\'' in line:
                python_files[python_file] = l_it
            elif f'\'{group}.{constant_name}\'' in line:
                group_python_files[python_file] = l_it

    # if none found we should keep a list
    if len(python_files) == 0:
        CC.cprint('\tConstant not found outside definition. Skipping',
                  colour='yellow')

        if len(group_python_files) > 0:
            CC.cprint(f'\t{group}.{constant_name} found in: ', colour='magenta')
            for python_file in group_python_files:
                lines = all_python_lines[python_file]
                line_number = group_python_files[python_file]
                print_entry(f'{group}.{constant_name}', python_file,
                            lines, line_number, colour='green')
            return all_python_lines, const_python_lines, updated_lines
        else:
            accept = input('Skip target? [Y/N]>>\t')
            if accept.strip().upper() in ['Y', 'YES']:

                return all_python_lines, const_python_lines, updated_lines
            else:
                global MISSING
                MISSING.append(constant_name)

    CC.cprint('')
    # -------------------------------------------------------------------------
    CC.cprint('Found instances:', colour='magenta')
    # print out these entries (removing the package path from python file
    for python_file in python_files:
        # get lines
        lines = all_python_lines[python_file]

        line_number = python_files[python_file]

        print_entry(constant_name, python_file, lines, line_number,
                    colour='green')
    CC.cprint('')
    # -------------------------------------------------------------------------
    # propose changes
    # -------------------------------------------------------------------------
    # 1: ask user to change variable name
    # -------------------------------------------------------------------------
    while True:
        qmsg = (f'Enter new constant name for "{constant_name}" '
                f'(leave blank to skip)>>\t')
        new_constant_name = input(qmsg)
        if len(new_constant_name) == 0:
            new_constant_name = str(constant_name)
            # add the group name onto the new_constant_name
            new_constant_name1 = f'{group}.{new_constant_name}'
            break

        if str(new_constant_name).upper() in ['Y', 'N', 'YES', 'NO']:
            CC.cprint('\tInvalid constant name', colour='red')
            continue
        # add the group name onto the new_constant_name
        new_constant_name1 = f'{group}.{new_constant_name}'

        CC.cprint(f'\n\nChanging {constant_name} to {new_constant_name1}',
                  colour='magenta')
        accept = input('Accept changes? (y/n)>>\t')
        if accept.strip().upper() in ['Y', 'YES']:
            break

    # -------------------------------------------------------------------------
    # 2: rename variables in all python codes
    # -------------------------------------------------------------------------
    # print out these entries (removing the package path from python file
    for python_file in python_files:
        CC.cprint('Before:', colour='magenta')

        # get lines
        old_lines = all_python_lines[python_file]

        old_line_number = python_files[python_file]

        print_entry(constant_name, python_file, old_lines, old_line_number,
                    colour='green', indent=4)

        CC.cprint('After:', colour='magenta')
        # get lines
        lines = all_python_lines[python_file]
        # get line number
        line_number = python_files[python_file]
        # get line
        line = lines[line_number]
        # new line
        new_line = line.replace(f'\'{constant_name}\'', f'\'{new_constant_name1}\'')
        # update line
        lines[line_number] = new_line
        # print entry
        print_entry(new_constant_name1, python_file, lines, line_number,
                    colour='blue', indent=4)
        # ask to accept changes
        qmsg = 'Accept changes? (y/n)>>\t'
        accept = input(qmsg)

        if accept.strip().upper() in ['Y', 'YES']:
            all_python_lines[python_file] = lines
            updated_lines[python_file] = lines
        else:
            raise StopLoop()

    # -------------------------------------------------------------------------
    # 3: rename variables in the const files
    # -------------------------------------------------------------------------
    found_in_const = False


    for python_file in const_python_lines:
        # get context
        context = const_python_lines[python_file]
        # skip files not found in context
        if constant_name not in ''.join(context):
            continue
        else:
            found_in_const = True

        # get the updated line
        uout = add_group_to_constant(constant_name, new_constant_name,
                                     context, python_file)
        updated_line, start_line, end_line = uout

        # deal with not having found the constant
        if start_line == 0 and end_line == 0:
            CC.cprint(f'\tNo entry found for {constant_name} in {python_file}',
                      colour='yellow')
            time.sleep(0.1)
            continue
        # get old entry
        old_entry = ('\n'.join(context).splitlines())

        if ''.join(old_entry) == ''.join(updated_line):
            CC.cprint(f'\tNo changes needed for {constant_name} in {python_file}',
                      colour='green')
            time.sleep(0.1)
            continue

        # print new entry
        print_entry(constant_name, python_file, old_entry,
                    line_start=start_line, line_end=end_line + 1, colour='green')
        # print new entry
        print_entry(new_constant_name, python_file, updated_line,
                    line_start=start_line, line_end=end_line + 1, colour='blue')

        # ask to accept changes
        qmsg = 'Accept changes? (y/n)>>\t'
        accept = input(qmsg)

        if accept.strip().upper() in ['Y', 'YES']:
            const_python_lines[python_file] = updated_line
            updated_lines[python_file] = updated_line
        else:
            raise StopLoop()

    # deal with never finding
    if not found_in_const:
        CC.cprint(f'\tNo entry found for {constant_name} in any const file',
                  colour='red')
        input('Press enter to skip')
        raise StopLoop()

    return all_python_lines, const_python_lines, updated_lines


def print_entry(entry, python_file, lines, line_number=0,
                line_start=None, line_end=None, colour='green',
                indent=0):
    # deal with indent
    if indent > 0:
        prefix = ' ' * indent
    else:
        prefix = ''

    # deal with not line start and end
    if line_start is None:
        line_start = max([0, line_number - 3])
    if line_end is None:
        line_end = min([line_number + 4, len(lines)])

    # truncated file name
    string_filename = python_file.replace(PACKAGE_PATH, '')

    CC.cprint(prefix + '>> ' + string_filename, colour)

    for line_it in range(line_start, line_end):
        # add an indicator for the line we are changing
        if line_it == line_number:
            sep = '+'
        else:
            sep = '|'
        # format the line
        fmt_line = lines[line_it].replace('\n', '')
        # print the message
        CC.cprint(prefix + f'{str(line_it):5s}{sep} {fmt_line}',
                  colour=colour, highlight_colour='red',
                  highlight_words=[entry])


def string_pos_to_line_pos(lines, start, end):
    # Initialize variables to keep track of the current character position
    current_pos = 0

    # Find the lines that contain the start and end positions
    start_line = end_line = None
    for i, line in enumerate(lines):
        line_start = current_pos
        line_end = current_pos + len(line)

        if start_line is None and line_start <= start <= line_end:
            start_line = i
        if end_line is None and line_start <= end <= line_end:
            end_line = i

        current_pos = line_end + 1  # Account for the '\n'

    return start_line, end_line


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # step 1: identify constants in default constants file
    constants_list = find_default_constants(INSTRUMENT_PATH +
                                            'default/constants.py')
    # total number of constants
    total_num_constants = len(constants_list)
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
    all_python_files, const_python_files = get_all_python_files(PACKAGE_PATH)
    # read all python files and stora the lines in memory in a dictionary
    all_python_lines = read_all_python_files(all_python_files)
    # read all const python files
    const_python_lines = read_all_python_files(const_python_files)
    # -------------------------------------------------------------------------
    # step 4: ask user for confirmation and edit files
    # -------------------------------------------------------------------------
    updated_lines = dict()
    updated_constants = []
    # get the number of constants
    num_constants = len(valid_constants_list.keys())
    # loop around constants, display the variable, ask for the new name, and
    # then confirm changes, then write changes to files
    for c_it, constant_name in enumerate(constants_list.keys()):

        # get percentage done
        perc = ((c_it + 1) / num_constants) * 100
        # reset next and stop
        next, stop = False, False
        # loop around so we can redo constant if needed
        while not next:
            # print the variable name
            CC.cprint('\n\n')
            CC.cprint(HEADER, colour='magenta')
            CC.cprint(f'Processing {constant_name} ({c_it+1} of '
                      f'{total_num_constants} ({perc:.2f} %)', colour='magenta')

            if constant_name not in valid_constants_list:
                CC.cprint(HEADER, colour='magenta')

                CC.cprint('No group/constant already fixed. Skipping',
                          colour='magenta')
                # reset next and stop
                next, stop = True, False
                continue
            # get the group name
            group_name = valid_constants_list[constant_name]
            CC.cprint(f'\tGroup: {group_name}', colour='magenta')
            CC.cprint(HEADER, colour='magenta')
            try:
                uout = update_constant(constant_name, all_python_lines,
                                       const_python_lines, updated_lines,
                                       group=group_name)
                # we need to update the input dictionaries so we can change the
                # next constant
                all_python_lines, const_python_lines, updated_lines = uout
                # this constant has been updated
                updated_constants.append(constant_name)
                # go to next entry
                next = True
            except (KeyboardInterrupt, StopLoop, Exception):
                # we need to ask what to do next
                qtime = True
                # loop until valid response
                while qtime:
                    # ask user if they want to continue
                    qmsg = ('\n\nRedo loop [R], Continue [C] or Stop [S]'
                            '\n>>\t')
                    accept = input(qmsg)
                    if accept.strip().upper() == 'R':
                        next, stop, qtime = False, False, False
                    elif accept.strip().upper() == 'C':
                        next, stop, qtime = True, False, False
                    elif accept.strip().upper() == 'S':
                        next, stop, qtime = True, True, False
                    else:
                        continue
        # deal with stopping
        if stop:
            break

    # -------------------------------------------------------------------------
    # step 5: save all updated line files
    # -------------------------------------------------------------------------
    for python_file in updated_lines:
        # print progress
        CC.cprint(f'\n\nWriting changes to {python_file}', colour='magenta')
        # get the lines
        lines = updated_lines[python_file]
        # write the lines to the file
        with open(python_file, 'w') as f:
            f.write('\n'.join(lines))



# =============================================================================
# End of code
# =============================================================================
