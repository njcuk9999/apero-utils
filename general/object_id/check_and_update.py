"""
Check object identification Google sheet and update the remote document.
"""
import utils as ut

sheet_id = '1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8'

# Load dataframe corresponding to google sheet
df_full = ut.get_full_sheet(sheet_id)

# Only check new inputs
df = df_full[~df_full['CHECKED']]

# Run ID checks
df = ut.check_id(df, replace=True)

# Get 2mass id for each gaia id
df = ut.twomass_from_gaia(df)

# Some ids are not in the catalog. Try to get them from simbad aliases
df = ut.twomass_from_simbad(df)

# Add newly checked entries to full database
df_full[~df_full['CHECKED']] = df


# Update sheet
old_df = ut.update_full_sheet(sheet_id, df_full)

old_df.to_csv('backup_sheet.csv', index=False)
