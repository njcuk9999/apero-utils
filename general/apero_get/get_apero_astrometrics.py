#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Get apero astrometrics

Created on 2023-06-14
Last updated 2023-06-14

@author: cook
"""
import numpy as np
import requests
from astropy.table import Table, vstack

# =============================================================================
# define variables
# =============================================================================
# url to download main table
main_url = ('https://docs.google.com/spreadsheets/d/'
            '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/gviz/'
            'tq?tqx=out:csv&gid=0')
# url to download pending table
pend_url = ('https://docs.google.com/spreadsheets/d/'
            '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/gviz/'
            'tq?tqx=out:csv&gid=623506317')
# object name column
gl_objcol = 'OBJNAME'


# =============================================================================
# define functions
# =============================================================================
def get_apero_astrometrics():
    """
    Get the APERO astrometrics database as a table

    Deals with main and pending table
    """
    # get the main table
    try:
        rawdata = requests.get(main_url)
        main_table = Table.read(rawdata.text, format='ascii')
    except Exception as _:
        raise ValueError('Cannot get table from {0}'.format(main_url))
    # get the pending table
    try:
        rawdata = requests.get(pend_url)
        pend_table = Table.read(rawdata.text, format='ascii')
    except Exception as _:
        raise ValueError('Cannot get table from {0}'.format(pend_url))
    # deal with overlap in pending table
    for _table in [pend_table]:
        # only do this if this table has some entries
        if len(_table) != 0:
            # make sure we have the object name column
            if gl_objcol in _table.colnames:
                pmask = ~np.in1d(_table[gl_objcol], main_table[gl_objcol])
                # add new columns to main table
                main_table = vstack([main_table, _table[pmask]])
    # return the table
    return main_table


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    # get the apero astrometrics table
    apero_table = get_apero_astrometrics()
