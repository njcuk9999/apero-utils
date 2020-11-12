import numpy as np
import gspread_pandas as gspd

import utils as ut

sheet_id = '1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8'

sh = gspd.spread.Spread(sheet_id, sheet='Info')

df = sh.sheet_to_df(index=0, unformatted_columns=ut.COLNAMES)
old_df = df.copy()

df['ALIASES'] = df['ALIASES'].str.split('|').apply(set).str.join('|')

# Double bracket to have Dataframe and not Series
sh.df_to_sheet(df[['ALIASES']], headers=True, start='H1', index=False)
