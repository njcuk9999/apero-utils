"""
Copy odometers from manually added list to automated list used by the DRS.

Important: Users need gspread_pandas with the credentials of an account that
has write access to both google sheets in order to run this code.
"""
import pandas as pd
import gspread_pandas as gspd

MAN_SHEET_ID = '13idu6tx-Zvzp-G9cjOSz-y5FH3E4a0T4Sdlw0235CNs'
AUTO_SHEET_ID = '1gvMp1nHmEcKCUpxsTxkx-5m115mLuQIGHhxJCyVoZCM'


skip_tabs = ['README']

sh_man = gspd.spread.Spread(MAN_SHEET_ID)

ws_names = [ws.title for ws in sh_man.sheets if ws.title not in skip_tabs]

colnames = ['ODOMETER', 'PP', 'RV', 'COMMENTS']
df_all = pd.DataFrame([], columns=colnames)

for ws_name in ws_names:

    # Load worksheet
    df_man = sh_man.sheet_to_df(index=0, sheet=ws_name)

    # If only one odometer, make sure not breaking code
    no_to = df_man['to'] == ''
    df_man.loc[no_to, 'to'] = df_man.loc[no_to, 'odometers']

    # Some sheets have total numbers, discard
    df_man = df_man[df_man.odometers != 'TOTAL']

    # Reformat into auto sheet columns
    zipped_odo = zip(
                     df_man['main reason'],
                     df_man['odometers'].astype(int),
                     df_man['to'].astype(int))
    df_auto = pd.DataFrame(
            [(o, False, True, c)
                for c, s, e in zipped_odo for o in range(s, e+1)],
            columns=colnames,
            )
    df_all = df_all.append(df_auto)

df_all = df_all.sort_values('ODOMETER')

# Load auto sheet
sh_auto = gspd.spread.Spread(AUTO_SHEET_ID)
df_auto = sh_auto.sheet_to_df(index=0, unformatted_columns=colnames)

# Add only new odometers
new = df_all.loc[~df_all.ODOMETER.isin(df_auto.ODOMETER)]
df_auto = df_auto.append(new)
df_auto = df_auto.sort_values('ODOMETER')

# Upload new list
# Putting replace=True such that whole sheet is overwritten
# Need to make sure that no other code edits the sheet anyway
# If it does, we must be careful to coordinate everything to avoid overwriting
sh_auto.df_to_sheet(df_auto, index=False, replace=True)
