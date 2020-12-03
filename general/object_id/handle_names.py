"""
Compare header names and CFHT
NOTE: Probably not optimized, but does the trick for a few 100 (unique) names.
"""
import re
import warnings
import numpy as np
import pandas as pd
import utils as ut
from astroquery.simbad import Simbad

# Known bad names
BAD_NAMES = [
        'sky.*_.*',
        '.*_sky.*',
        'sky',
        'SKY',
        'FakeTarget',
        'Test.*',
        '.*Test',
        'Engineering'
        ]
BAD_NAMES = '|'.join(BAD_NAMES)  # pipe symbol as OR in regex

# Known rejected names
REJ_NAMES = [
        'Moon.*',
        'Neptune.*',
        'Saturn.*',
        'Venus.*',
        ]
REJ_NAMES = '|'.join(REJ_NAMES)  # pipe symbol as OR in regex

path_to_data = '/home/vandal/Documents/spirou/obj_db/object_info_full.csv'
sheet_id = '1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8'

# Load all data
df_sheet = ut.get_full_sheet(sheet_id)
df_local = pd.read_csv(path_to_data)

# Get unique names
unames = df_local.OBJECT.drop_duplicates()

# Find bad names with regex and filter out
bad_names = unames[unames.str.fullmatch(BAD_NAMES)]
good_names = unames[~unames.isin(bad_names)]

# Match names that are already in list or aliases (can clear them)
new_names = good_names[~good_names.isin(df_sheet.OBJECT)]
all_alias = df_sheet.ALIASES.str.cat(sep='|').split('|')
new_names = new_names[~new_names.isin(all_alias)]

# When possible, find gaia dr2 id in simbad and cross check with sheet
df_new = pd.DataFrame([], columns=df_sheet.columns)
df_new['OBJECT'] = new_names.values
type_mask = df_new.OBJECT.str.match('\\([OBAFGKM]\\)')
df_new_no_types = df_new.copy()  # remove spectra type when there
df_no_type_names = df_new_no_types.OBJECT[type_mask].str.slice(start=4)
df_new_no_types.loc[type_mask, 'OBJECT'] = df_no_type_names
df_new_no_types = ut.check_simbad(df_new_no_types)
mod_cols = ['GAIADR2ID', 'FOUND', 'COMMENTS']
df_new[mod_cols] = df_new_no_types[mod_cols]  # put search results in real df
match_mask = df_new.GAIADR2ID.isin(df_sheet.GAIADR2ID)
df_alias = df_new[match_mask]  # Store matched IDs for alias later
df_new = df_new[~match_mask]
df_test = df_new[df_new.GAIADR2ID.isnull()]  # Store no ID for testing
df_new = df_new[df_new.GAIADR2ID.notnull()]  # Store unmatched ids as new

# Remove too short names as they will "over match"
all_new_names = df_test.OBJECT
# new_names = all_new_names[all_new_names.str.len() > 2]
new_names = all_new_names.copy()
new_modif = new_names.copy()

new_modif = new_modif.str.strip()

# If spectral type before, remove for comparison
type_mask = new_modif.str.match('\\([OBAFGKM]\\)')
new_modif.loc[type_mask] = new_modif[type_mask].str.slice(start=4)

# some formatting
new_modif = new_modif.str.upper().str.replace(' ', '').replace('_', '')
new_modif = new_modif.str.replace('*', '')

# Check potential new matches with string manipulation
cum_match_mask = np.zeros_like(new_modif)
print('This shows potential matches.')
print('The actual matches must be added to alias column manually')
print('False positives must be added to OBJECT column (if true object)')
for ind, row in df_sheet.iterrows():

    # Some basic formatting to increase matches
    aliases = row.ALIASES.upper().replace(' ', '').replace('_', '')
    aliases = aliases.replace('*', '')

    # Check if newname in alias
    match_mask = new_modif.isin(aliases.split('|'))
    cum_match_mask = match_mask | cum_match_mask
    found = new_names[match_mask]
    nomatch = new_names[~match_mask]

    if found.size > 0:
        print(
                'OBJECT {} has potential new aliases, formatted match: {}'
                .format(row.OBJECT, found.values)
              )

    # Need the re.escape to add \\ before [] and (), and others
    esc_aliases = '|'.join(re.escape(s) for s in aliases.split('|'))
    match_mask = new_modif.str.contains(esc_aliases)
    cum_match_mask = cum_match_mask | match_mask
    found_all = np.array(found)
    found = new_names[match_mask]
    found = found[~found.isin(found_all)]
    nomatch = new_names[~match_mask]

    if found.size > 0:
        print(
                'OBJECT {} has potential new aliases, new contains(alias): {}'
                .format(row.OBJECT, found.values)
              )

    # Somewhat arbitrary len > 2 to avoid over-matching short names
    match_mask = np.array([name in aliases
                           if len(name) > 2 else False
                           for name in new_modif])
    cum_match_mask = cum_match_mask | match_mask
    found_all = np.append(found_all, found)
    found = new_names[match_mask]
    found = found[~found.isin(found_all)]
    nomatch = new_names[~match_mask]

    if found.size > 0:
        print(
                'OBJECT {} has potential new aliases, new in any alias: {}'
                .format(row.OBJECT, found.values)
              )

    found_all = np.append(found_all, found)

    if found_all.size > 0:
        print()

# Objects without str match
no_match = new_names[~cum_match_mask].sort_values(ignore_index=True,
                                                  key=lambda v: v.str.upper(),
                                                  )
df_test = df_test[~cum_match_mask].sort_values('OBJECT',
                                               ignore_index=True,
                                               key=lambda v: v.str.upper(),
                                               )

# Append matched GAIA IDs as aliases automatically
al_arr = df_alias.GAIADR2ID.values
sh_arr = df_sheet.GAIADR2ID.values
inds = np.where(al_arr[:, None] == sh_arr)[-1]
new_al = df_sheet.loc[inds].ALIASES.str.cat(df_alias.OBJECT.values, sep='|')
df_sheet.loc[inds, 'ALIASES'] = new_al

# Among new gaia ids, make sure no duplicates
df_new['ALIASES'] = df_new['OBJECT']
# funs = dict(zip(df_new.columns, ['first'] * len(df_new.columns)))
# funs['ALIASES'] = '|'.join
# df_new = df_new.groupby('GAIADR2ID').agg(funs).reset_index(drop=True)
# df_new = df_new.sort_values('OBJECT',
#                             ignore_index=True,
#                             key=lambda v: v.str.upper(),
#                             )


# Among no_match, see if any are found in simbad
type_mask = df_test.OBJECT.str.match('\\([OBAFGKM]\\)')
df_test_no_types = df_test.copy()  # remove spectra type when there
df_no_type_names = df_test_no_types.OBJECT[type_mask].str.slice(start=4)
df_test_no_types.loc[type_mask, 'OBJECT'] = df_no_type_names
with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # mute simbad warnings
    for i, row in df_test_no_types.iterrows():
        tbl = Simbad.query_objectids(row.OBJECT)
        if tbl is None:
            continue
        aliases = '|'.join(tbl.to_pandas()['ID'].str.decode('utf-8').tolist())
        df_test.loc[i, 'ALIASES'] = aliases
        df_test.loc[i, 'COMMENTS'] = 'Found aliases in SIMBAD'

# RUN THIS IF WANT TO ADD FIELDS WITH EMTPY id
# df_new = df_new.append(df_test[df_test.ALIASES.notnull()])
# df_new = df_new.sort_values('OBJECT',
#                             ignore_index=True,
#                             key=lambda v: v.str.upper(),
#                             )
# df_test = df_test[df_test.ALIASES.isnull()]

# Update some columns for new objects
df_new['CHECKED'] = df_new.GAIADR2ID.notnull()
df_new['FOUND'] = df_new.GAIADR2ID.notnull()
df_new['MANUAL_EDIT'] = False

# Update sheet with new objects for which 100% sure
if df_new.shape[0] > 0:
    print('Updating info sheet with new objects')
    df_sheet = df_sheet.append(df_new, ignore_index=True).sort_values(
                                                'OBJECT',
                                                ignore_index=True,
                                                key=lambda v: v.str.upper(),
                                                )
    backup_df = ut.update_full_sheet(sheet_id, df_sheet)

# Rejected objects
rej_names = df_test.loc[df_test.OBJECT.str.fullmatch(REJ_NAMES), 'OBJECT']
df_test = df_test[~df_test.OBJECT.str.fullmatch(REJ_NAMES)]
rej_names = pd.DataFrame(rej_names)  # Make DF
df_rej = ut.get_full_sheet(sheet_id, ws='Rejected')  # Load known rejected
rej_mask = rej_names.OBJECT.isin(df_rej.OBJECT)
df_rej = df_rej.append(rej_names[~rej_mask], ignore_index=True).sort_values(
                                                'OBJECT',
                                                ignore_index=True,
                                                key=lambda v: v.str.upper(),
                                                )
if not rej_mask.all():
    print('Updating rejected sheet with new objects')
    df_back = ut.update_full_sheet(sheet_id, df_rej, ws='Rejected')

# Get gaia id from position match where possible
dup_mask = df_local.OBJECT.duplicated(keep='last')
df_local_one = df_local[~dup_mask]
test_mask = df_local.OBJECT.isin(df_test.OBJECT)
df_with_id = ut.gaia_crossmatch(df_local_one[test_mask])
ord_ids = df_with_id.set_index('OBJECT').loc[df_test.OBJECT].GAIADR2ID
df_test['GAIADR2ID'] = ord_ids.values  # sorted above so safe to use values
id_mask = df_test.GAIADR2ID.notnull()
# df_test['ALIASES'] = df_test.OBJECT
df_checked = df_test.copy()
df_checked[id_mask] = ut.add_aliases(df_checked[id_mask])
al_mask = df_checked.ALIASES.notnull()
df_checked['CHECKED'] = False
df_checked['FOUND'] = False
df_checked['MANUAL_EDIT'] = False
df_checked[al_mask] = ut.check_id(df_checked[al_mask])

# Split found and not found objects
df_found = df_checked[df_checked.FOUND]
df_not_found = df_checked[~df_checked.FOUND]

# Addd found to new
df_new = df_new.append(df_found)
df_new = df_new.sort_values('OBJECT',
                            ignore_index=True,
                            key=lambda v: v.str.upper(),
                            )

# Add new to sheet
if df_new.shape[0] > 0:
    print('Updating info sheet with new objects')
    df_sheet = df_sheet.append(df_new, ignore_index=True).sort_values(
                                                'OBJECT',
                                                ignore_index=True,
                                                key=lambda v: v.str.upper(),
                                                )
    id_mask = df_sheet.GAIADR2ID.notnull()
    sheet_no_id = df_sheet[~id_mask]  # Should be empty if sheet OK
    df_sheet = df_sheet[id_mask]
    funs = dict(zip(df_new.columns, ['first'] * len(df_new.columns)))
    funs['ALIASES'] = '|'.join
    df_sheet = df_sheet.groupby('GAIADR2ID').agg(funs).reset_index(drop=True)
    df_sheet = df_sheet.append(sheet_no_id)
    df_sheet = df_sheet.sort_values('OBJECT',
                                    ignore_index=True,
                                    key=lambda v: v.str.upper(),
                                    )
    unique_al = df_sheet['ALIASES'].str.split('|').apply(set).str.join('|')
    df_sheet['ALIASES'] = unique_al
    backup_df = ut.update_full_sheet(sheet_id, df_sheet)
    backup_df.to_csv('backup_sheet.csv', index=False)

# Update objects that need sorting to a separate sheet
need_sort = ut.get_full_sheet(sheet_id, ws='Need Sorting')
need_sort = need_sort.append(df_not_found)
backup_need = ut.update_full_sheet(sheet_id, need_sort, ws='Need Sorting')
