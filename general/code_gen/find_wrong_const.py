#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-02-21 at 09:57

@author: cook
"""
from astropy.table import Table
import os
import re


# =============================================================================
# Define variables
# =============================================================================
HISTORY_FILE = '/home/cook/good_matches.fits'
# -----------------------------------------------------------------------------
PATH = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX'

EXCLUDE_FILES = ['constants.py', 'config.py', 'default_help.py',
                 'default_text.py', 'keywords.py',
                 'nirps_ha_help.py', 'nirps_ha_text.py',
                 'nirps_he_help.py', 'nirps_he_text.py',
                 'spirou_help.py', 'spirou_text.py']

# =============================================================================
# Define functions
# =============================================================================
COLOURS = dict()
COLOURS['BLACK'] = '\033[90;1m'
COLOURS['RED'] = '\033[1;91;1m'
COLOURS['GREEN'] = '\033[92;1m'
COLOURS['YELLOW'] = '\033[1;93;1m'
COLOURS['BLUE'] = '\033[94;1m'
COLOURS['MAGENTA'] = '\033[1;95;1m'
COLOURS['CYAN'] = '\033[1;96;1m'
COLOURS['WHITE'] = '\033[97;1m'
COLOURS['END'] = '\033[0;0m'


def cprint(message, kwargs=None):
    if kwargs is None:
        kwargs = dict()
    for colour in COLOURS:
        kwargs[colour.lower()] = COLOURS[colour]
    print(message.format(**kwargs))


def find_and_replace_in_file(file_path, search_string, replacement_string):
    """Finds and replaces occurrences of a search string in a file with user confirmation."""

    # Regex pattern to match both 'search_string' and "search_string"
    pattern = re.compile(rf'(["\']){re.escape(search_string)}\1')

    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    modified = False  # Track if the file is modified
    new_lines = lines[:]  # Copy original lines

    for i, line in enumerate(lines):
        matches = list(pattern.finditer(line))  # Find all matches in the line

        if matches:
            # Show context (previous, current, next lines)
            before = lines[i - 1] if i > 0 else ""
            after = lines[i + 1] if i < len(lines) - 1 else ""

            sline = line.strip()
            sline = sline.replace('{', '{{').replace('}', '}}')
            sline = sline.replace(search_string,
                                  '{red}' + search_string + '{green}')

            kwargs = dict(file_path=file_path,
                          before=before.strip(), after=after.strip(),
                          search_string=search_string,
                          replacement_string=replacement_string)

            message = ("{green}\n" + "-" * 40)
            message += ("\nFile: {yellow}{file_path}{green}")
            message += ("\n\nContext:")
            message += ("\n  {before}")
            message += ("\n> " + sline)
            message += ("\n  {after}")
            message += ("\n\nReplace {red}'{search_string}'{green} with ")
            message += ("{blue}'{replacement_string}'{green}?{end}")

            cprint(message, kwargs)

            # Ask for confirmation
            user_input = input("(y/n): ").strip().lower()
            if user_input.lower() in ['y', 'yes']:
                # Replace all occurrences in this line
                new_lines[i] = pattern.sub(rf'\1{replacement_string}\1', line)
                modified = True

    if modified:
        # Write the modified content back to the file
        with open(file_path, 'w', encoding='utf-8') as file:
            file.writelines(new_lines)


def get_all_module_python_files(path):
    out_files = []
    for root, dirs, files in os.walk(path):
        # get the root minus the path1 + os.sep
        rel_root = root.replace(path, '')
        if rel_root.endswith(os.sep):
            rel_root = rel_root[:-len(os.sep)]
        if rel_root.startswith(os.sep):
            rel_root = rel_root[len(os.sep):]

        for filename in files:
            if filename.endswith('.py'):
                out_files.append(os.path.join(rel_root, filename))
    return out_files


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load the history file
    history_table = Table.read(HISTORY_FILE)

    # get the old names
    old_names = history_table['OLD']
    # get the new names
    new_names = history_table['NEW']
    # get the comments
    comments = history_table['COMMENT']


    # get all python files
    pyfiles = get_all_module_python_files(PATH)

    # loop around all old names
    for it in range(len(old_names)):
        # get old name
        old_name = old_names[it]
        # get new name
        new_name = new_names[it]
        # print progress
        pargs = [old_name, it + 1, len(old_names),
                 (it + 1) / len(old_names) * 100]
        print('\n\nProcessing old name: {0} [{1}/{2}] {3:3f}%'.format(*pargs))

        for jt, pyfile in enumerate(pyfiles):
            # skip these files
            if os.path.basename(pyfile) in EXCLUDE_FILES:
                continue
            print('\tProcessing file {0}/{1}'.format(jt + 1, len(pyfiles)))
            # get fullpath
            fullpath = os.path.join(PATH, pyfile)
            # find and replace in file
            find_and_replace_in_file(fullpath, old_name, new_name)


# =============================================================================
# End of code
# =============================================================================
