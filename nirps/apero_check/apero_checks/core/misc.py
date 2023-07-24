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
from typing import Any, Dict, Optional, Tuple

from astropy.time import TimeDelta
from astropy import units as uu

from apero_checks.core import base
from apero_checks.core import io
from apero_checks.core import parameters

# =============================================================================
# Define variables
# =============================================================================
# version, date, author
__VERSION__ = base.__VERSION__
__DATE__ = base.__DATE__
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
    parser.add_argument('--test_run', '--test', '--test_name', '--testrun',
                        '--testname', type=str, default='None',
                        help='If running a single test define the name of that'
                             ' test here (the column name)')
    parser.add_argument('--today', action='store_true', default=False,
                        help='If we do not give an --obsdir we can give '
                             'this argument to always use today\'s date as '
                             'the observation directory (only works when '
                             'observation directories are named by date) in '
                             'the form YYYY-MM-DD')
    parser.add_argument('--yesterday', action='store_true', default=False,
                        help='If we do not give an --obsdir we can give '
                             'this argument to always use yesterdays\'s date '
                             'as the observation directory (only works when '
                             'observation directories are named by date) in '
                             'the form YYYY-MM-DD. Note --today takes '
                             'priority over --yesterday')
    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return vars(args)


def add_cmd_arg(args, key, value, null=None) -> Tuple[Any, str]:
    # set source
    source = 'FUNC'
    # define the null values
    null_values = {null, str(null), None}
    # if value is equal to the null value then we check "args"
    if value in null_values:
        # if the "args" is not in the null values then we use the "args" value
        if args[key] not in null_values:
            value = args[key]
            source = 'CMD'
    # force value to null if in any of the null values (i.e. a string version
    #   of the null)
    if value in null_values:
        return null, 'NULL'
    # otherwise return our value
    else:
        return value, source


def add_source(name: str, value: Any, source: str):
    if value is None:
        return ''
    if source in [None, 'NULL']:
        return ''
    else:
        return f'{name}={value}   [{source}]'

def load_params(yaml_file: Optional[str] = None,
                obsdir: Optional[str] = None,
                test_name: Optional[str] = None,
                today: bool = False,
                yest: bool = False) -> Dict[str, Any]:
    # set up the return dictionary
    params = dict()
    # string for printing variables used
    sources = dict()
    # load default parameters
    for key in parameters.parameters:
        params[key] = copy.deepcopy(parameters.parameters[key].value)
    # -------------------------------------------------------------------------
    # read from command line
    args = get_args()
    # get yaml file from cmd args
    yaml_file , yaml_source = add_cmd_arg(args, 'yaml', yaml_file)
    sources['yaml'] = add_source('yaml', yaml_file, yaml_source)
    # get obs dir from cmd args
    params['obsdir'] , obsdir_source = add_cmd_arg(args, 'obsdir', obsdir)
    sources['obsdir'] = add_source('obsdir', params['obsdir'], obsdir_source)
    # deal with no observation directory (and having today defined
    if params['obsdir'] is None:
        # deal with getting today (bool) from input or cmd args
        today , today_source = add_cmd_arg(args, 'today', today, null=False)
        sources['today'] = add_source('today', today, today_source)
        # if today is True then we set obsdir to today
        if today:
            # we set it to now YYYY-MM-DD
            # TODO: this will not be correct for anyone who uses
            #       observation directories that are not YYYY-MM-DD
            params['obsdir'] = base.AstropyTime.now().iso.split(' ')[0]

    if params['obsdir'] is None:
        # deal with getting today (bool) from input or cmd args
        yest, yest_source = add_cmd_arg(args, 'yesterday', yest, null=False)
        sources['yest'] = add_source('yest', yest, yest_source)
        # if today is True then we set obsdir to today
        if yest:
            # we set it to now YYYY-MM-DD
            # TODO: this will not be correct for anyone who uses
            #       observation directories that are not YYYY-MM-DD
            timenow = base.AstropyTime.now()
            timeyest = timenow - TimeDelta(1 * uu.day)
            params['obsdir'] = timeyest.iso.split(' ')[0]

    # get test name from cmd args
    params['test_name'], tname_source = add_cmd_arg(args, 'test_run', test_name)
    sources['test_name'] = add_source('today', today, tname_source)
    # -------------------------------------------------------------------------
    # load from yaml file
    yaml_params = io.read_yaml(yaml_file)
    # push into params
    for key in yaml_params:
        # deal with invalid key
        if key not in parameters.parameters:
            continue
        # check value from yaml file
        value, value_source = parameters.parameters[key].check(yaml_params[key])
        # add sources to keys we want to report
        if parameters.parameters[key].report:
            sources[key] = add_source(key, value, 'YAML')
        # otherwise add key to params
        params[key] = copy.deepcopy(value)
    # -------------------------------------------------------------------------
    # print message on which observation directory we are processing
    msg = f'Parameters used:'
    log_msg(msg, level='info')
    # print which yaml we are using
    for source in sources:
        msg = '\t - ' + sources[source]
        if len(sources[source]) > 0:
            log_msg(msg, level='info')
    # -------------------------------------------------------------------------
    # return params
    return params


def splash(codename):
    msg1 = '*' * 50
    log_msg(msg1, color='cyan')
    msg2 = '*  ' + codename
    log_msg(msg2, color='magenta')
    msg3 = '*' * 50
    msg3 += f'\nVersion={__VERSION__}   Date={__DATE__}'
    msg3 += '\n' + '*' * 50
    log_msg(msg3, color='cyan')


def end_msg():
    # print message on which observation directory we are processing
    msg = '*' * 50
    msg += f'\nCode finished successfully'
    msg += '\n' + '*' * 50
    log_msg(msg, level='info')


def log_msg(message, level: str = '', color: Optional[str] = None):

    # we want to print all lines with prefix
    if '\n' in message:
        messages = message.split('\n')
    else:
        messages = [message]

    # get time now
    timenow = base.AstropyTime.now()
    # loop around all messages
    print_string = ''
    for it, msg in enumerate(messages):
        if it > 0:
            prefix = '\n'
        else:
            prefix = ''
        # set up print string
        print_string += f'{prefix}{timenow} | {msg}'

    # deal with level
    if color is None and level is not None:
        if level == 'error':
            color = 'red'
        elif level == 'warning':
            color = 'yellow'
        elif level == 'test':
            color = 'magenta'
        elif level == 'info':
            color = 'green'
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
        raise base.AperoChecksError(print_string)
    else:
        print(start + print_string + end)


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
