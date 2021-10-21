"""
Update information of remote objects.

@author: vandalt
"""
from argparse import ArgumentParser

import pandas as pd
import gspread_pandas as gspd

import utils as ut
import maintenance as mtn

# Set pattern to fits files
parser = ArgumentParser(
        description='Update remote info in the Google Sheet.'
        )

parser.add_argument(
        '-p',
        '--pattern',
        default=mtn.RAW_PATTERN,
        help='Pattern of fits files.',
        )
parser.add_argument(
        '-f',
        '--file',
        default=None,
        help='Path to local file with observation info.',
        )
parser.add_argument(
        '-s',
        '--save',
        action='store_true',
        help='Save (write or update) local file if true.',
        )
parser.add_argument(
        '-i',
        '--sheet-id',
        dest='sheet_id',
        default=mtn.SHEET_ID,
        help='ID of the Google sheet.',
        )
clargs = parser.parse_args()

# Fetch object info from headers (or load local csv file with this info)
if clargs.pattern is not None:
    df_loc = mtn.get_object_info(clargs.pattern, local_file=clargs.file)
    if clargs.save:
        df_loc.to_csv(mtn.LOCAL_FILE, index=False)
elif clargs.file is not None:
    df_loc = pd.read_csv(clargs.file)

# Get remote sheet and load main Info tab to dataframe
sh = gspd.spread.Spread(clargs.sheet_id)
df_sheet = ut.get_full_sheet(sh, ws='Info', index=0)

# Load Need Sorting sheet
df_sort = ut.get_full_sheet(sh, ws='Need Sorting')

# Split checked
df_new = df_sort[df_sort.CHECKED]
df_sort = df_sort[~df_sort.CHECKED]

# Update sorting sheet
ut.update_full_sheet(sh, df_sort, ws='Need Sorting')

# Add checked objects to sheet
df_sheet = ut.sheet_append(df_new, df_sheet)

# Get IDs and aliases
df_sheet = mtn.get_2mass_id(df_sheet)
df_sheet = mtn.find_aliases_simbad(df_sheet)
ut.update_full_sheet(sh, df_sheet, ws='Info')

# Update sheet with TEFF and RV information
mtn.teff_and_rv(sh, df_sheet, df_loc)
