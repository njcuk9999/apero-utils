#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-10-28 at 13:32

@author: cook
"""
import sys

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import warnings


# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def function1():
    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":

    # properties dict (to be put in ParamDict)
    params = dict()
    # -------------------------------------------------------------------------
    # python version
    import sys

    major = sys.version_info.major
    minor = sys.version_info.minor
    micro = sys.version_info.micro

    pyversion = f'{major}.{minor}.{micro}'
    # add to params
    params['PYVERSION'] = pyversion
    # -------------------------------------------------------------------------
    # get all modules
    try:
        from pip._internal.operations import freeze
    except ImportError:  # pip < 10.0
        from pip.operations import freeze

    pkgs = freeze.freeze()

    for pkg in list(pkgs):
        if '==' in pkg:
            key, version = pkg.split('==', 1)
            params[f'PYTHONMOD_{key}'] = version.strip()
        elif '@' in pkg:
            key, version = pkg.split('@', 1)
            params[f'PYTHONOTHER_{key}'] = version.strip()
    # -------------------------------------------------------------------------

    # get git info
    # TODO: requires git

    import os

    # from apero.core import constants
    # params = constants.load()
    # drs_root = params['DRS_ROOT']

    drs_root = '/scratch2/spirou/drs-bin/apero-drs-spirou-07XXX/apero'

    try:
        from git import Repo
        repo = Repo(os.path.dirname(drs_root))
        branch_name = str(repo.active_branch.name)
        branch_hash = str(repo.head.object.hexsha)
        del repo
    except Exception as _:
        branch_name = 'Unknown'
        branch_hash = 'Unknown'

    params['GIT_BRANCH'] = branch_name
    params['GIT_HASH'] = branch_hash
    # -------------------------------------------------------------------------

# =============================================================================
# End of code
# =============================================================================
