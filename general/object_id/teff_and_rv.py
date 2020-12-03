"""
For each object in list, load TEFF and RV, then upload to sheet
"""
import numpy as np
import pandas as pd
import gspread_pandas as gspd

import utils as ut

path_to_data = '/home/vandal/Documents/spirou/obj_db/object_info_full.csv'
sheet_id = '1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8'


def vals_and_counts(df, kwd, threshold=None):

    vals = df[kwd]
    counts = vals.value_counts()
    vals = vals.dropna().unique()

    # reject unphysical values
    if threshold is not None:
        vals = vals[np.abs(vals) < threshold]

    vals = vals[np.nonzero(vals)]
    counts = counts[vals].values

    return vals, counts


def pipelist(arr):
    arr = np.array(arr)
    return '|'.join(list(arr.astype(str)))


def tryaddlast(alist, arr):

    try:
        alist.append(arr[-1])
    except IndexError:
        alist.append(np.nan)


# Load google sheet
sh = gspd.spread.Spread(sheet_id)
sh.open_sheet('Info')
df_id = sh.sheet_to_df(index=0, unformatted_columns=ut.COLNAMES)
old_df_id = df_id.copy()


# Load local file with info of each data file
df = pd.read_csv(path_to_data)

rv_list = []
teff_list = []

rv_all = []
teff_all = []
rv_all_count = []
teff_all_count = []


# Loop through objects and get values
# There may be a better way with pandas, but fast enough with 121 object
for remote_object in df_id['OBJECT']:
    objrows = df[df['OBJECT'] == remote_object]
    rv, rv_count = vals_and_counts(objrows, 'OBJRV', threshold=(0, 2000))
    teff, teff_count = vals_and_counts(objrows,
                                       'OBJTEMP',
                                       # threshold=(1000, np.inf),
                                       )

    # Append latest value to "main" list (pandas unique used, so no sorting)
    tryaddlast(rv_list, rv)
    tryaddlast(teff_list, teff)

    # Keep pipe-separated string of all values and counts
    rv_all.append(pipelist(rv))
    rv_all_count.append(pipelist(rv_count))
    teff_all.append(pipelist(teff))
    teff_all_count.append(pipelist(teff_count))


# Update RV and TEFF in Info worksheet
update_mask = ~df_id['MANUAL_EDIT']
rv_teff_update = np.array([rv_list, teff_list]).T[update_mask]
df_id.loc[update_mask, ['RV', 'TEFF']] = rv_teff_update
rv_mask, teff_mask = ~np.isnan(rv_teff_update).T
df_id.loc[update_mask & rv_mask, 'RV_REF'] = 'HEADER[CFHT]'
df_id.loc[update_mask & teff_mask, 'TEFF_REF'] = 'HEADER[CFHT]'
sh.df_to_sheet(df_id, index=False, replace=True)

# Update the maintenance sheet (no mask for this one)
sh.open_sheet('Maintenance')
maint = sh.sheet_to_df(index=0)
cols = ['RV_ALL', 'RV_COUNTS', 'TEFF_ALL', 'TEFF_COUNTS']
maint[cols] = np.array([rv_all, rv_all_count, teff_all, teff_all_count]).T
sh.df_to_sheet(maint, index=False, replace=True)
