#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2021-11-10 at 16:04

@author: cook
"""
from astroquery.utils.tap.core import TapPlus
from astroquery.simbad import Simbad
import warnings

# =============================================================================
# Define variables
# =============================================================================
# simbad additional columns
SIMBAD_COLUMNS = ['ids']

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def query_object(rawobjname):
    # set up the TapPlus class
    # simbad = TapPlus(url=URL)

    # execute the job
    # sjob = simbad.launch_job(QUERY.format(rawobjname))

    # get the results
    # table = sjob.get_results()

    # get results
    with warnings.catch_warnings(record=True) as _:
        # add ids column
        for simbad_column in SIMBAD_COLUMNS:
            Simbad.add_votable_fields(simbad_column)
        # query simbad
        table = Simbad.query_object(rawobjname)

    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # step 1: Is object in database



    # ----------------------------------------------------------------------
    # step 2: if not in database see if we can find it in simbad
    _ = query_object('Gl699')




# =============================================================================
# End of code
# =============================================================================
