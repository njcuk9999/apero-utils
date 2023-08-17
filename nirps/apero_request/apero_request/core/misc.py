#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
import argparse
from typing import Any, Dict, Optional

from apero_request.core import base
from apero_request.core import parameters
from apero_request.core import io

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
    parser = argparse.ArgumentParser(description='Run apero requests')
    # add obs dir
    parser.add_argument('yaml', type=str, default='None',
                        help='The profiles yaml to use')
    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return vars(args)



def load_params():
    """
    Load the parameters from defaults, command line and yaml file
    :return:
    """
    # start with defaults
    params = parameters.parameters.copy()
    # read from command line
    args = get_args()
    # get yaml file
    yaml = args['yaml']
    # load yaml file
    yaml_params = io.read_yaml(yaml)
    # push into params
    for key in yaml_params:
        params[key] = yaml_params[key].copy()
    # return parameters
    return params


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
        raise base.AperoRequestError(print_string)
    else:
        print(start + print_string + end)


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



# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
