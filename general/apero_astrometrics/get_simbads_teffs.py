#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2024-01-30 at 14:39

@author: cook
"""
from astroquery.utils.tap.core import TapPlus

# =============================================================================
# Define variables
# =============================================================================
SIMBAD_QUERY = 'http://simbad.cds.unistra.fr/simbad/sim-tap'


TEFF_QUERY = """
SELECT 
       main_id AS "Main identifier",
       mesFe_H.teff as teff,
       mesFe_H.bibcode as bibcode,
       mesFe_H.log_g as log_g,
       mesFe_H.fe_h as fe_h
       
FROM basic JOIN ident ON ident.oidref = oid
INNER JOIN mesFe_H ON mesFe_H.oidref = oid
WHERE id = '{OBJNAME}';
"""

VSINI_QUERY = """
SELECT 
       main_id AS "Main identifier",
       mesRot.vsini as vsini,
       mesRot.vsini_err as vsini_err,
       mesRot.bibcode as bibcode
       
FROM basic JOIN ident ON ident.oidref = oid
INNER JOIN mesRot ON mesRot.oidref = oid
WHERE id = '{OBJNAME}';
"""

# -----------------------------------------------------------------------------
OBJNAME = 'Gl699'


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # create a TapPlus instance
    tap = TapPlus(url=SIMBAD_QUERY)
    # execute the teff query
    teff_job = tap.launch_job(TEFF_QUERY.format(OBJNAME=OBJNAME))
    # fetch the results
    teff_results = teff_job.get_results()
    # execute the vsini query
    vsini_job = tap.launch_job(VSINI_QUERY.format(OBJNAME=OBJNAME))
    # fetch the results
    vsini_results = vsini_job.get_results()


# =============================================================================
# End of code
# =============================================================================
