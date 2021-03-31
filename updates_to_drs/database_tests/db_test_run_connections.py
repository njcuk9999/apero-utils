#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-03-31

@author: cook
"""
from apero.base import drs_db
import numpy as np
import time


# =============================================================================
# Define variables
# =============================================================================
# add you login here
case = 2
if case == 1:
    host = 'localhost'
    username = 'cook'
    password = 'neilcook'
    tablename = 'ea_test'
    database = 'spirou'
else:
    host = 'rali'
    username = 'spirou'
    password = 'Covid19!'
    tablename = 'ea_test'
    database = 'spirou'

# number of add rows to try
N_ADDS = 9000


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # wait to start so we can run on multiple windows at once
    print('Waiting to start')
    time.sleep(2)
    # make mysql database
    database = drs_db.MySQLDatabase(host, username, password, database, tablename)
    # make sure database exists
    database.add_database()
    # add our test table
    database.add_table(tablename, ['X', 'Y', 'Z'], [str, str, str])
    # loop around N_ADDS times
    for it in range(N_ADDS):
        print('\nProcessing row {0}'.format(it + 1))
        # set up some values
        values = np.array([it, it * 2, it * 3])
        values = values.astype(str)
        # add rows
        database.add_row(values, tablename)

# =============================================================================
# End of code
# =============================================================================
