#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on {DATE}

@author: cook
"""
import os
import argparse
import copy
from hashlib import blake2b
from typing import Any, Dict, List, Optional, Union
import smtplib
from email.mime.text import MIMEText
from socket import gethostname

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
    parser.add_argument('profiles', type=str, default='None',
                        help='Only use these profiles (ignore any others)')
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
        params[key] = copy.deepcopy(yaml_params[key])
    # push the profiles argument into params
    if args['profiles'] != 'None':
        profiles = args['profiles'].split(',')
        # clean up profile names
        params['profiles'] = []
        # loop around profiles
        for profile in profiles:
            params['filter profiles'].append(profile.strip())
        # if for some reason we have no profiles set it back to None
        if len(params['profiles']) == 0:
            params['filter profiles'] = None

    # return parameters
    return params


def log_msg(params, message, level: str = '', color: Optional[str] = None,
            log_only: bool = False):

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

    if not log_only:
        # print message
        if level == 'error':
            raise base.AperoRequestError(print_string)
        else:
            print(start + print_string + end)

    # construct log file path
    log_file = os.path.join(params['local path'], 'apero_request.log')
    # save to log file
    with open(log_file, 'a') as f:
        f.write(print_string + '\n')


def splash(params, codename):
    msg1 = '*' * 50
    log_msg(params, msg1, color='cyan')
    msg2 = '*  ' + codename
    log_msg(params, msg2, color='magenta')
    msg3 = '*' * 50
    msg3 += f'\nVersion={__VERSION__}   Date={__DATE__}'
    msg3 += '\n' + '*' * 50
    log_msg(params, msg3, color='cyan')


def end_msg(params):
    # print message on which observation directory we are processing
    msg = '*' * 50
    msg += f'\nCode finished successfully'
    msg += '\n' + '*' * 50
    log_msg(params, msg, level='info')


def generate_arg_checksum(source: Union[List[str], str],
                          ndigits: int = 10) -> str:
    """
    Take a list of strings or a string and generate a unique hash from
    them

    :param source: list of strings or string - the string to generate the hash
                   from
    :param ndigits: int, the size of the hash (in characters) default is 10

    :return: str, the hash
    """
    # set function name
    # _ = display_func('generate_arg_checksum', __NAME__)
    # flatten list into string
    if isinstance(source, list):
        source = ' '.join(source)
    # need to encode string
    encoded = source.encode('utf')
    # we want a hash of 10 characters
    digest = blake2b(encoded, digest_size=ndigits)
    # create hash
    _hash = digest.hexdigest()
    # return hash
    return str(_hash)


def send_email(message: str,
               people: Optional[Union[List[str], str]],
               sender_address: str,
               subject: Optional[str] = None,):
    """Simple wrapper around `MIMEText()` and `smtplib` to send an email.

    Note: For this to work, your account must be setup to send emails at
    the command line, so it will likely crash on a laptop without prior config.

    :param message: Text for the body of the email.
    :param people: People to email. No email sent if None.
    :param subject: Subject of the email
    """
    # deal with no people

    if people is None or len(people) == 0:
        return
    if len(people) == 1 and people[0] is None:
        return

    if subject is None:
        subject = "APERO REQUEST"

    msg = MIMEText(message)
    msg["Subject"] = subject
    msg["From"] = f"APERO on {gethostname()} <{sender_address}>"

    if isinstance(people, str):
        msg["To"] = people
    else:
        msg["To"] = ", ".join(people)

    with smtplib.SMTP("localhost") as s:
        s.send_message(msg)


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
