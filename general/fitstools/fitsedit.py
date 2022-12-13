#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-12-12 at 13:03

@author: cook
"""
import argparse
import os
import sys

from astropy.io import fits


# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def get_args():
    parser = argparse.ArgumentParser(description='Fits header edit')
    # file name is a positional argument
    parser.add_argument(type=str, dest='filename',
                        help='Filename to operate on')

    # mode
    parser.add_argument('--mode', type=str,
                        choices=['edit', 'add', 'remove'],
                        default='None',
                        help='The way we are manipulating the file')
    # add key
    parser.add_argument('--key', type=str,
                        default='None', help='The header key to edit')
    # add value
    parser.add_argument('--value', type=str, help='The value to set --key to',
                        default='None')
    # set extension number
    parser.add_argument('--ext', type=int, default=0,
                        help='The extension number of the fits file')

    # load arguments with parser
    args = parser.parse_args()
    # return arguments
    return args


def edit_key(params):
    """
    Edit a key and save the file

    :param params: arguments passed from command line

    :return:
    """
    with fits.open(params.filename) as hdu:

        hdr = hdu[params.ext].header

        if params.key in hdr:
            print(f'HDR[{params.key}]     {hdr[params.key]}-->{params.value}')
            hdr[params.key] = params.value

            hdu[params.ext].header = hdr
            print(f'File updated: {params.filename}')
            hdu.writeto(params.filename, overwrite=True)
        else:
            print(f'Key: {params.key} not found in header. Use --mode=add')


def add_key(params):
    """
    Add a key and save the file

    :param params: arguments passed from command line

    :return:
    """
    with fits.open(params.filename) as hdu:

        hdr = hdu[params.ext].header

        if params.key in hdr:
            print(f'Key: {params.key} found in header. Use --mode=edit')
        else:
            print(f'HDR[{params.key}]     Added: {params.value}')
            hdr[params.key] = params.value
            hdu[params.ext].header = hdr
            print(f'File updated: {params.filename}')
            hdu.writeto(params.filename, overwrite=True)


def remove_key(params):
    """
    Remove a key and save the file

    :param params: arguments passed from command line

    :return:
    """
    with fits.open(params.filename) as hdu:

        hdr = hdu[params.ext].header

        if params.key in hdr:
            print(f'HDR[{params.key}]     Removed: {params.value}')
            del hdr[params.key]
            hdu[params.ext].header = hdr
            print(f'File updated: {params.filename}')
            hdu.writeto(params.filename, overwrite=True)
        else:
            print(f'Key: {params.key} not found in header.')



# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    # get arguments
    params = get_args()

    # ----------------------------------------------------------------------
    # Deal with file not existing
    if not os.path.exists(params.filename):
        print(f'File not found: {params.filename}')
        sys.exit()

    # ----------------------------------------------------------------------
    # edit file name
    if params.key != 'None' and params.value != 'None':

        # evaluate value
        try:
            params.value = eval(params.value)
        except Exception as _:
            params.value = str(params.value).strip()

        if params.mode == 'edit':
            edit_key(params)
        if params.mode == 'add':
            add_key(params)

    if params.key != 'None' and params.mode == 'remove':
        remove_key(params)


# =============================================================================
# End of code
# =============================================================================
