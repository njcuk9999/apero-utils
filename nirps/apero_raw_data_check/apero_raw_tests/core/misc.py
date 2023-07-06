#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 14:51

@author: cook
"""
import argparse
import copy
from typing import Any, Dict, Optional

from apero_raw_tests.core import base
from apero_raw_tests.core import io
from apero_raw_tests.core import parameters

# =============================================================================
# Define variables
# =============================================================================
__VERSION__ = base.__VERSION__
__AUTHOR__ = base.__AUTHOR__


# =============================================================================
# Define functions
# =============================================================================
def get_args() -> Dict[str, Any]:
    """
    Define the command line arguments

    :return: argparse namespace object containing arguments
    """
    parser = argparse.ArgumentParser(description='Run apero raw tests')
    # add obs dir
    parser.add_argument('yaml', type=str, default='None',
                        help='The profile yaml to use')
    parser.add_argument('--obsdir', type=str, default='None',
                        help='Observation directory name(s) separated by '
                             'commas')
    parser.add_argument('--test_run', type=str, default='None',
                        help='If running a single test define the name of that'
                             ' test here (the column name)')
    parser.add_argument('--today', action='store_true', default=False,
                        help='If we do not give an --obsdir we can give '
                             'this argument to always use today\'s date as '
                             'the observation directory (only works when '
                             'observation directories are named by date) in '
                             'the form YYYY-MM-DD')
    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return vars(args)


def add_cmd_arg(args, key, value, null=None):
    # define the null values
    null_values = {null, str(null), None}
    # if value is equal to the null value then we check "args"
    if value in null_values:
        # if the "args" is not in the null values then we use the "args" value
        if args[key] not in null_values:
            value = args[key]
    # force value to null if in any of the null values (i.e. a string version
    #   of the null)
    if value in null_values:
        return null
    # otherwise return our value
    else:
        return value


def load_params(yaml_file: Optional[str] = None,
                obsdir: Optional[str] = None,
                test_name: Optional[str] = None,
                today: bool = False) -> Dict[str, Any]:
    # set up the return dictionary
    params = dict()
    # load default parameters
    for key in parameters.parameters:
        params[key] = copy.deepcopy(parameters.parameters[key].value)
    # -------------------------------------------------------------------------
    # read from command line
    args = get_args()
    # get yaml file from cmd args
    yaml_file = add_cmd_arg(args, 'yaml', yaml_file)
    # get obs dir from cmd args
    params['obsdir'] = add_cmd_arg(args, 'obsdir', obsdir)
    # deal with no observation directory (and having today defined
    if params['obsdir'] is None:
        # deal with getting today (bool) from input or cmd args
        today = add_cmd_arg(args, 'today', today, null=False)
        # if today is True then we set obsdir to today
        if today:
            # we set it to now YYYY-MM-DD
            # TODO: this will not be correct for anyone who uses
            #       observation directories that are not YYYY-MM-DD
            params['obsdir'] = base.AstropyTime.now().iso.split(' ')[0]
    # get test name from cmd args
    params['test_name'] = add_cmd_arg(args, 'test_run', test_name)
    # -------------------------------------------------------------------------
    # load from yaml file
    yaml_params = io.read_yaml(yaml_file)
    # push into params
    for key in yaml_params:
        # deal with invalid key
        if key not in parameters.parameters:
            continue
        # check value from yaml file
        value = parameters.parameters[key].check(yaml_params[key])
        # otherwise add key to params
        params[key] = copy.deepcopy(value)
    # return params
    return params


def log_msg(message, level: str = '', color: Optional[str] = None):
    # set up print string
    print_string = f'{base.AstropyTime.now()} | {message}'
    # deal with level
    if color is None and level is not None:
        if level == 'error':
            color = 'red'
        elif level == 'warning':
            color = 'yellow'
        elif level == 'info':
            color = 'cyan'
        else:
            color = 'green'
    # add color (if required)
    if color in base.COLOURS:
        start = base.COLOURS[color]
        end = base.COLOURS['ENDC']
    else:
        start, end = '', ''
    # print message
    if level == 'error':
        raise base.AperoRawTestsError(print_string)
    else:
        print(start + message + end)


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
