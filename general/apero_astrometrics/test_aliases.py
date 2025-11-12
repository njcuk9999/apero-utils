#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-09-10 at 12:00

@author: cook
"""
from apero.core import constants
from apero.core.core import drs_database

# =============================================================================
# Define variables
# =============================================================================
TARGET_NAME = 'DOTAU'


# =============================================================================
# Define functions
# =============================================================================


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # load params
    params = constants.load()
    # load astrometric database
    astrodb = drs_database.AstrometricDatabase(params)
    astrodb.load_db()
    # set up condition
    condition = f'OBJNAME="{TARGET_NAME}"'
    # get aliases for target
    aliases = astrodb.get_entries('ALIASES', condition=condition)[0].split('|')

    # print output
    print(f'Target name = {TARGET_NAME}')
    for a_it, alias in enumerate(aliases):
        print(f'Alias [{a_it + 1}/{len(aliases)}] = {alias}')




# =============================================================================
# End of code
# =============================================================================
