#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Code to quickly get the astrometric database

Note: Do not spam the request too heavily you could get your IP blocked

@author: cook
"""
import requests
import pandas as pd
from astropy.table import Table
from io import StringIO

def fetch_astrom():
    # =============================================================================
    # Define variables
    # =============================================================================
    MAIN_URL = ('https://docs.google.com/spreadsheets/d/'
                '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
                'export?format=csv&gid=0')
    PENDING_URL = ('https://docs.google.com/spreadsheets/d/'
                   '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
                   'export?format=csv&gid=623506317')

    # =============================================================================
    # Start of code
    # =============================================================================

    # get main table
    main_request = requests.get(MAIN_URL)
    # open main table
    main_dataframe = pd.read_csv(StringIO(main_request.text))
    # get pending table
    pending_request = requests.get(PENDING_URL)
    # open pending table
    pending_dataframe = pd.read_csv(StringIO(pending_request.text))
    # merge these keeping all rows from main table and adding non-repeating
    #  rows from pending table
    astrom_dataframe = pd.concat([main_dataframe,
                                  pending_dataframe]).drop_duplicates()

    # if you prefer an astropy table
    astrom_table = Table.from_pandas(astrom_dataframe)

    # we add a column to indicate if the file is a telluric
    astrom_table['TELLURIC'] = False
    for i in range(len(astrom_table)):
        if 'TELLURIC' in astrom_table['KEYWORDS'][i]:
            astrom_table['TELLURIC'][i] = True
    # =============================================================================
    # End of code
    # =============================================================================

    return astrom_table