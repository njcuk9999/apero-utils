"""
Perform full update of the Google sheet with local observations

@author: vandalt
"""
import pandas as pd
import gspread_pandas as gspd

import maintenance as mtn
import utils as ut

# Set pattern to fits files
# fits_pattern = mtn.RAW_PATTERN
fits_pattern = None
local_file = '/home/vandal/Documents/spirou/obj_db/object_info.csv'
save_loc = True
sheet_id = mtn.SHEET_ID

# Fetch object info from headers (or load local csv file with this info)
if fits_pattern is not None:
    df_loc = mtn.get_object_info(fits_pattern, local_file=local_file)
    if save_loc:
        df_loc.to_csv('object_info.csv', index=False)
elif local_file is not None:
    df_loc = pd.read_csv(local_file)

# Get remote sheet and load main Info tab to dataframe
sh = gspd.spread.Spread(sheet_id)
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
