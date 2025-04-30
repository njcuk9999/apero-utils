#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-02-19 at 15:10

@author: cook
"""
import os
import re
import string
from typing import List

import numpy as np
from astropy.table import Table

# =============================================================================
# Define variables
# =============================================================================
PATH1 = '/scratch2/apero/apero_9c7fb3a17768ed4c49f63c8088c8fac2f69f653c'
PATH2 = '/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX'
# default constant files
CFILES = ['apero-drs/apero/instruments/default/constants.py',
          'apero-drs/apero/instruments/default/config.py']

EXCLUDE_FILES = ['constants.py', 'config.py', 'default_help.py',
                 'default_text.py', 'keywords.py',
                 'nirps_ha_help.py', 'nirps_ha_text.py',
                 'nirps_he_help.py', 'nirps_he_text.py',
                 'spirou_help.py', 'spirou_text.py']

# -----------------------------------------------------------------------------
# we'll get this recursively later
FILEPATH = 'apero-drs/apero/core/drs_file.py'
# string match pattern
STRING_PATTERN = re.compile(r"(\".*?\"|'.*?')")
# allowed punctuation
allowed_punc = list(string.punctuation.replace('_', '').replace('.', ''))

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


def extract_strings_from_code(code: List[str], codename: str):
    # Regex pattern to match single ('...') and double ("...") quoted strings
    extracted_data = []

    for line_number, line in enumerate(code, start=1):
        strings = STRING_PATTERN.findall(line)# Extract all strings
        clean_line = STRING_PATTERN.sub('', line)  # Remove strings from line
        # only add if we have a string
        if len(strings) > 0:
            for string_ in strings:
                string_ = string_.strip('\'').strip('\"')
                if len(string_) == 0:
                    continue
                if len(string_) > 0:
                    my_rtn = []
                    my_rtn += [line_number]
                    my_rtn += [string_]
                    my_rtn += [clean_line.strip()]
                    my_rtn += [codename]
                    extracted_data.append(tuple(my_rtn))

    return extracted_data


def get_line(line_number, code):
    if line_number < 0:
        line_number = 0
    if line_number >= len(code):
        line_number = len(code) - 1
    return code[line_number]


def bad_constant(mystring):
    # check if string is a bad constant
    if not str(mystring).isupper():
        return True
    if np.sum(np.in1d(list(mystring), allowed_punc)) > 0:
        return True
    if mystring.startswith('KW_'):
        return True
    if len(mystring) < 2:
        return True
    if mystring.upper() in ['NONE', 'TRUE', 'FALSE']:
        return True
    if mystring[0].isdigit():
        return True
    return False


def get_code_matches(filepath, constants):
    # ---------------------------------------------------------------------
    codefile1 = os.path.join(PATH1, filepath)
    codefile2 = os.path.join(PATH2, filepath)

    # load both files
    with open(codefile1, 'r') as f:
        code1 = f.read()
    with open(codefile2, 'r') as f:
        code2 = f.read()

    codelines1 = code1.splitlines()
    codelines2 = code2.splitlines()

    results1 = extract_strings_from_code(codelines1, codefile1)
    results2 = extract_strings_from_code(codelines2, codefile2)

    matches1 = []
    matches2 = []

    for result1 in results1:
        # get the string(s) and clean_line1
        line_number1, string1, clean_line1, _ = result1
        # get the line before and line after
        line_before1 = get_line(line_number1, codelines1)
        line_after1 = get_line(line_number1, codelines1)
        # need to remove strings from the line before and after
        line_before1 = STRING_PATTERN.sub('', line_before1)
        line_after1 = STRING_PATTERN.sub('', line_after1)
        # conditions on just result 1
        if bad_constant(string1):
            continue
        # string1 must be in constants
        if string1 not in constants:
            continue
        # find any clean_lines that match a line in result1
        for result2 in results2:
            line_number2, string2, clean_line2, _ = result2
            # get the line before and line after
            line_before2 = get_line(line_number2, codelines2)
            line_after2 = get_line(line_number2, codelines2)
            # need to remove strings from the line before and after
            line_before2 = STRING_PATTERN.sub('', line_before2)
            line_after2 = STRING_PATTERN.sub('', line_after2)
            # Only match code blocks that are the same but with strings
            # replaced
            if bad_constant(string2):
                continue
            if string1 == string2:
                continue
            if clean_line1.strip() != clean_line2.strip():
                continue
            if line_before1.strip() != line_before2.strip():
                continue
            if line_after1.strip() != line_after2.strip():
                continue

            matches1.append(result1)
            matches2.append(result2)

    return matches1, matches2


def get_all_module_python_files(path) -> List[str]:
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
    # compile a list of known constants
    constants = []
    for cfile in CFILES:
        for path in [PATH1, PATH2]:
            with open(os.path.join(PATH1, cfile), 'r') as f:
                code = f.read()
            codelines = code.splitlines()
            results = extract_strings_from_code(codelines, cfile)
            for result in results:
                _, string_, _, _ = result
                # only keep upper case strings
                if not (str(string_).isupper()):
                    continue
                # only keep strings that are not punctuation
                if np.sum(np.in1d(list(string_), allowed_punc)) > 0:
                    continue
                # add to constants
                constants.append(string_.strip())
    # remove duplicates
    constants = list(set(constants))
    # -------------------------------------------------------------------------
    # find all python files in the directories that match between path1 and
    # path2
    files1 = get_all_module_python_files(PATH1)
    files2 = get_all_module_python_files(PATH2)
    # find those files1 that match files2
    file_matches = []
    for file1 in files1:
        # ignore files in CFILES
        if file1 in CFILES:
            continue
        # ignore excluded files
        if os.path.basename(file1) in EXCLUDE_FILES:
            continue
        # otherwise we have a match to test
        if file1 in files2:
            file_matches.append(file1)
    # -------------------------------------------------------------------------
    # get matches
    matches1 = []
    matches2 = []
    # loop around all files and get code matches
    for f_it, file_match in enumerate(file_matches):
        pargs = [file_match, f_it+1, len(file_matches)]
        print('Analysing file: {0} [{1}/{2}]'.format(*pargs))
        matches1_, matches2_ = get_code_matches(file_match, constants)
        print('\t\tFound {0} matches'.format(len(matches1_)))
        matches1 += matches1_
        matches2 += matches2_
    # -------------------------------------------------------------------------
    # filter out duplicates (where string1 and string2 match)
    matches1_ = []
    matches2_ = []
    used = []
    # loop around all matches
    for it in range(len(matches1)):
        string_match = (matches1[it][1], matches2[it][1])
        if string_match in used:
            continue
        used.append(string_match)
        matches1_.append(matches1[it])
        matches2_.append(matches2[it])
    # -------------------------------------------------------------------------
    used1 = []
    used2 = []
    pos_match = []
    comment = []
    bad_match = []
    loop = True
    # go through matches
    counter = 0
    last_operation = None

    while counter < len(matches1_):
        match1 = matches1_[counter]
        match2 = matches2_[counter]
        line_number1, string1, clean_line1, code_file1 = match1
        line_number2, string2, clean_line2, code_file2 = match2

        # do not ask again about good matches
        if (string1, string2) in pos_match:
            counter += 1
            continue
        # do not ask again about bad matches
        if (string1, string2) in bad_match:
            counter += 1
            continue
        # if we have already assigned the before or after then we don't do it
        # again
        if string1 in used1:
            counter += 1
            continue
        if string2 in used2:
            counter += 1
            continue

        kwargs = dict(string1=string1, string2=string2,
                      it=counter+1, total=len(matches1_),
                      percentage=100*(counter+1)/len(matches1_))

        cprint('\n\n{green}Match: {red}{string1}{green}-->{blue}{string2}'
               ' {yellow}[{it}/{total}] {percentage:.3f}%{end}',
               kwargs=kwargs)

        while True:
            # ask the user if this match is good
            answer = input('Is this a good match? (y/n/e/+/a/s/b): ')
            answer = answer.strip().upper()
            # deal with back
            if answer in ['B', 'BACK']:
                counter = max([0, counter-2])
                if last_operation == 'pos':
                    pos_match = pos_match[:-1]
                    comment = comment[:-1]
                    used1 = used1[:-1]
                    used2 = used2[:-1]
                elif last_operation == 'bad':
                    bad_match = bad_match[:-1]
                break
            # deal with skip
            if answer in ['S', 'SKIP']:
                answer_skip = input('\tSkip this match? (y/n): ')
                if answer_skip.upper().strip() in ['Y', 'YES']:
                    break
                else:
                    continue
            # deal with a (add)
            if answer in ['A', 'ADD']:
                answer_add1 = input('\tAdd BERFORE: ')
                answer_add2 = input('\tAdd AFTER: ')
                answer_comment = input('\tAdd COMMENT: ')
                answer_confirm = input('Confirm {0}-->{1} (y/n): ')

                if answer_confirm.upper().strip() in ['Y', 'YES']:
                    if len(answer_add1) > 0 and len(answer_add2) > 0:
                        if answer_add1 not in used1 and answer_add2 not in used2:
                            used1.append(answer_add1)
                            used2.append(answer_add2)
                            pos_match.append((answer_add1, answer_add2))
                            comment.append(answer_comment)
                            # Now ask again about previous match
                            continue
            # deal with yes
            elif answer in ['Y', 'YES']:
                pos_match.append((string1, string2))
                comment.append('')
                used1.append(string1)
                used2.append(string2)
                last_operation = 'pos'
                break
            # deal with exit
            elif answer in ['E', 'EXIT']:
                answer_exit = input('\tExit? (y/n): ')
                if answer_exit.upper().strip() in ['Y', 'YES']:
                    counter = np.inf
                    break
                else:
                    continue
            # deal with more info requested
            elif answer in ['+']:
                # open code1
                with open(os.path.join(PATH1, code_file1), 'r') as f:
                    code1 = f.read()
                codelines1 = code1.splitlines()
                codestr = '\n'.join(codelines1[line_number1-10:line_number1+10])
                kwargs = dict(code=codestr)
                cprint('{green}{code}{end}', kwargs=kwargs)
                continue
            # deal with no
            elif answer in ['N', 'NO']:
                bad_match.append((string1, string2))
                last_operation = 'bad'
                break
            # otherwise try again
            else:
                print('Invalid input')
        # update the counter
        counter += 1
    # -------------------------------------------------------------------------
    table_dict = dict(OLD=[], NEW=[], COMMENT=[])
    # write good matches to a fits table
    for it in range(len(pos_match)):
        string1, string2 = pos_match[it]
        table_dict['OLD'].append(string1)
        table_dict['NEW'].append(string2)
        table_dict['COMMENT'].append(comment[it])
    # convert to astropy table
    table = Table(table_dict)
    # write table to file
    table_path = os.path.abspath(os.path.realpath('good_matches.fits'))
    print('Writing good matches to: {0}'.format(table_path))
    table.write(table_path, overwrite=True)


# =============================================================================
# End of code
# =============================================================================
