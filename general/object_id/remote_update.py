"""
Update information of remote objects.

@author: vandalt
"""
import pandas as pd
import gspread_pandas as gspd

import utils as ut
import maintenance as mtn

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
