#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-04-29 at 15:00

@author: cook
"""
from typing import List
import importlib.util
from pathlib import Path

from apero.core import constants
from apero.core.core import drs_database
from apero.tools.recipes.bin import apero_get
from apero.tools.recipes.bin import apero_remove
from apero.tools.recipes.bin import apero_processing

# =============================================================================
# Define variables
# =============================================================================
# backup directory
BACKUP_PATH = '/cosmos99/spirou/misc/persi-backup'
# List object names
OBJECT_NAMES = ['TOI2120', 'TRAPPIST1']
# path to persistence code
PERSISTENCE_CODE = '/cosmos99/spirou/misc/persi-code/persitools/persicorr.py'
# path to the persistence map
PERSIFILE = '/cosmos99/spirou/misc/persi-code/'
# path to the calibrations
CALIBDIR = '/cosmos99/spirou/apero-data/spirou_offline/calib/'
# -----------------------------------------------------------------------------
# test mode
TEST = True

# =============================================================================
# Define functions
# =============================================================================
def run_persicorr(targets: List[str], test: bool = False):
    """
    Run the persistence code for all targets
    :param targets:
    :param test:
    :return:
    """
    # get persi code
    persi_code = get_persi_code()
    # get parameters
    params = constants.load()
    # get the file index database
    findexdb = drs_database.FileIndexDatabase(params)
    findexdb.load_db()
    # set up the condition to get files
    condition = 'KW_OUTPUT="EXT_E2DS_FF" '
    # add object names
    target_conditions = []
    for target in targets:
        target_conditions.append(f'KW_OBJNAME="{target}"')
    # push target conditons
    condition += '(' + ' OR '.join(target_conditions) + ')'
    # get database table
    files = findexdb.get_entries('ABSPATH', condition=condition)

    if test:
        print('Running persicorr:')
        for filename in files:
            print(f'\t - {filename}')
        return

    # run persi code wrapper
    persi_code.main(files, path_to_persifile=PERSIFILE,
                    mode='e2dsff', do_plot=False, replace=True,
                    calibdir=CALIBDIR)


def get_persi_code():
    """
    Get the Persistence code module
    :return:
    """
    # get the module name from the persistence code variable
    file_path = Path(PERSISTENCE_CODE)
    module_name = file_path.stem
    # load the module
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # return the module
    return module


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # Step 1: Back up files
    # ----------------------------------------------------------------------
    apero_get.main(outpath=BACKUP_PATH, objnames=','.join(OBJECT_NAMES),
                   outtypes='EXT_E2DS_FF', fibers='AB', test=TEST)

    # ----------------------------------------------------------------------
    # Step 2: Remove all traces of target
    # ----------------------------------------------------------------------
    apero_remove.main(objnames=','.join(OBJECT_NAMES), test=TEST)

    # ----------------------------------------------------------------------
    # Step 3: Run apero processing up to and including extraction
    # ----------------------------------------------------------------------
    apero_processing.main(runfile='persi_part1.ini',
                          science_targets=','.join(OBJECT_NAMES), test=TEST)

    # ----------------------------------------------------------------------
    # Step 4: Run persistence correction
    # ----------------------------------------------------------------------
    run_persicorr(targets=OBJECT_NAMES, test=TEST)

    # ----------------------------------------------------------------------
    # Step 5: Run apero processing telluric correction onwards
    # ----------------------------------------------------------------------
    apero_processing.main(runfile='persi_part2.ini',
                          science_targets=','.join(OBJECT_NAMES), test=TEST)

# =============================================================================
# End of code
# =============================================================================
