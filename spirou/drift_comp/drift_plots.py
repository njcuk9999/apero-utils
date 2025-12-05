#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CODE DESCRIPTION HERE

Created on 2020-08-2020-08-14 12:48

@author: cook
"""
import os
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, join, vstack, MaskedColumn

# =============================================================================
# Define variables
# =============================================================================

PATH = '/data/spip/misc/download/andres_drift/'

DRIFT_FILES = 'cal_drift_FP_FP_{fiber}.fits'

SCI_FIBERS = ['AB', 'A', 'B']
REF_FIBERS = ['C']

DATE_KEY = 'MJDMID'


# =============================================================================
# Define functions
# =============================================================================
def read_drift_files() -> Dict[str, Table]:
    """
    Read drift files from disk
    :return:
    """
    # storage for output
    drift_files = dict()
    # loop around fibers
    for fiber in SCI_FIBERS + REF_FIBERS:
        # get drift file for fiber
        basename = DRIFT_FILES.format(fiber=fiber)
        # push into drift table
        drift_files[fiber] = Table.read(os.path.join(PATH, basename))

    return drift_files


def equalize_tables(tables: Dict[str, Table], date_column='date'):
    # 1. Get the full union of all dates
    all_dates = set()
    for tkey in tables:
        all_dates.update(tables[tkey][date_column])
    all_dates = sorted(all_dates)

    # 2. Make a master table with all dates
    master = Table({date_column: all_dates})

    # 3. Join each table to the master on date, filling missing values with np.nan
    result_tables = dict()
    for tkey in tables:
        # Use a left join to preserve all dates from the master
        merged = join(master, tables[tkey], keys=date_column, join_type='left')

        # Replace masked values with np.nan for numeric columns
        for col in merged.colnames:
            if col == date_column:
                continue
            if isinstance(merged[col], MaskedColumn):
                if np.issubdtype(merged[col].dtype, np.number):
                    merged[col] = merged[col].filled(np.nan)
                else:
                    merged[col] = merged[col].filled(None)
        # push into results dict
        result_tables[tkey] = merged

    return result_tables


def plot_drift(drift: Dict[str, Table]):

    plt.plot(drift['AB'][DATE_KEY],
             drift['AB']['RV'] - drift['C']['RV'])

    plt.show(block=True)


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # read drift files
    drift_files = read_drift_files()
    # equalize drift files (based on filename)
    drift_files = equalize_tables(drift_files, DATE_KEY)

    # plot drift
    plot_drift(drift_files)

# =============================================================================
# End of code
# =============================================================================
