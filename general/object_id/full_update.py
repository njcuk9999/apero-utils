"""
Perform full update of the Google sheet with local observations

@author: vandalt
"""
from argparse import ArgumentParser

import pandas as pd
import gspread_pandas as gspd

import maintenance as mtn
import utils as ut

# Set pattern to fits files
parser = ArgumentParser(
        description='Update the Google Sheet with new objects.'
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

# Get unique list of local objects
df_loc_unique = df_loc.drop_duplicates(subset=['OBJECT'], keep='last')
loc_names = df_loc_unique.OBJECT.copy()

# Handle rejected and bad names
loc_names = mtn.bad_and_rejected(loc_names, sh=sh)

# Filter out names or aliases already in sheet
loc_names = loc_names[~loc_names.isin(df_sheet.OBJECT)]
all_alias = df_sheet.ALIASES.str.cat(sep='|').split('|')
loc_names = loc_names[~loc_names.isin(all_alias)]

# Search names in Simbad and update based on Gaia ID
test_names, df_sheet = mtn.check_gaia(loc_names, df_sheet)

# For next steps, put test_names in a dataframe with same columns as sheet
df_test = pd.DataFrame(test_names, columns=df_sheet.columns)
df_test.CHECKED = False
df_test.FOUND = False
df_test.MANUAL_EDIT = False
df_loc_unique = df_loc_unique[df_loc_unique.OBJECT.isin(test_names)]

# Position cross-match for other gaia IDs
df_test, df_sheet = mtn.gaia_match_position(df_loc_unique, df_test, df_sheet)

# Update full Info sheet (reload to make sure have latest)
df_new = df_sheet.copy()
df_sheet = ut.get_full_sheet(sh, ws='Info', index=0)
df_sheet = ut.sheet_append(df_new, df_sheet)
df_sheet = mtn.get_2mass_id(df_sheet)
df_sheet = mtn.find_aliases_simbad(df_sheet)
ut.update_full_sheet(sh, df_sheet, ws='Info')

# Update sheet with TEFF and RV information
mtn.teff_and_rv(sh, df_sheet, df_loc)

# Do string comparisons to find possible matches
df_test = mtn.str_match(df_test, df_sheet)

# Upload remaining test objects to Need Sorting
try:
    df_sort = ut.get_full_sheet(sh, ws='Need Sorting')
    df_sort = ut.sheet_append(df_test, df_sort)
except AttributeError:
    df_sort = df_test
ut.update_full_sheet(sh, df_sort, ws='Need Sorting')
