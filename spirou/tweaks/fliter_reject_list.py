#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-04-25

@author: cook
"""
from astropy.table import Table
import numpy as np
import requests
from tqdm import tqdm

from apero.core import constants

# =============================================================================
# Define variables
# =============================================================================
# define standard google base url
GOOGLE_BASE_URL = ('https://docs.google.com/spreadsheets/d/{}/gviz/'
                   'tq?tqx=out:csv&gid={}')

OUTFILENAME = 'filtered_reject_list.csv'

# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # get params
    params = constants.load()
    # get parameters from params
    sheet_id = params['REJECT_LIST_GOOGLE_SHEET_URL']
    worksheet = params['REJECT_LIST_GSHEET_MAIN_LIST_ID']
    # construct url for worksheet
    url = GOOGLE_BASE_URL.format(sheet_id, worksheet)
    # get data using a request
    try:
        rawdata = requests.get(url)
    except Exception as e:
        raise ValueError(f'URL Error: {url}\n\t{type(e)}: {str(e)}')
    # open table
    table = Table.read(rawdata.text, format='ascii')
    # get unique entries
    uids = np.unique(table['IDENTIFIER'])
    # storage
    keep_rows = np.zeros(len(table), dtype=bool)
    # loop around uids
    for uid in tqdm(uids):
        # find first instance of uid
        first_pos = np.where(table['IDENTIFIER'] == uid)[0][0]
        # keep first position
        keep_rows[first_pos] = True

    # cut down table
    cuttable = table[keep_rows]
    # save as csv (for opening)
    cuttable.write(OUTFILENAME, format='csv')

# =============================================================================
# End of code
# =============================================================================
