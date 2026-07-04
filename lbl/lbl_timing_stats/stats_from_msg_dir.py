#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2026-07-03 at 11:28

@author: cook
"""
import glob
import os
import re
from typing import Tuple, Union

from astropy.table import Table
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
PATH = '/project/6102120/apero/spirou_data/internal/spirou_offline_07/msg/processing/APEROG-PID-00017824133614703310-TAMQ_apero_processing_group/other'
# -----------------------------------------------------------------------------
# generate table from all these files
FILE_PATTERN = '*_apero_lbl_*_spirou.log'
# Regex string line to look for in files
# Regex breakdown:
# (apero_lbl_.*?_spirou) -> Group 1: Captures the full recipe name
# \t\(                  -> Matches the tab and opening parenthesis
# ([0-9.]+)             -> Group 2: Captures the float number (time)
STRING_FILE_PATTERN = r"(apero_lbl_.*?_spirou) has been successfully completed\t\(([0-9.]+)"

# =============================================================================
# Define functions
# =============================================================================
def get_stats(log_file_path: str) -> Tuple[Union[str, None], Union[float, None]]:


    pattern = re.compile(STRING_FILE_PATTERN)

    # Table Header
    print(f"{'Recipe Name':<25} | {'Time (seconds)':<15}")
    print("-" * 43)

    try:
        with open(log_file_path, "r") as file:
            for line in file:
                match = pattern.search(line)
                if match:
                    full_recipe_name = match.group(1)
                    time_float = float(match.group(2))

                    # Split the recipe name at the *first* underscore only
                    # e.g., "apero_lbl_compute_spirou" -> "apero"
                    short_recipe_name = str(full_recipe_name.split("_", 1)[0])

        return short_recipe_name, time_float

    except FileNotFoundError:
        return None, None


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # get all valid files in directory
    files = glob.glob(os.path.join(PATH, FILE_PATTERN))

    # storage for output
    storage_dict: dict = dict(RECIPE=[], TIME_TAKEN=[])

    for filename in tqdm(files):

        _rname, _rtime = get_stats(filename)

        if _rname is not None:
            storage_dict['RECIPE'].append(_rname)
            storage_dict['TIME_TAKEN'].append(_rtime)

    # convert to astropy table
    stats_table = Table(storage_dict)



# =============================================================================
# End of code
# =============================================================================
