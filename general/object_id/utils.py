"""
Utility functions to manage object ID sheet.

@author: vandalt
"""
import warnings
import requests
import numpy as np
import pandas as pd
from astropy.table import Table
from astroquery.mast import Catalogs
from astroquery.simbad import Simbad
from astroquery.utils.tap.core import TapPlus

import astro_object as ao

COLNAMES = ['OBJECT',
            'GAIADR2ID',
            'TWOMASSID',
            'RV',
            'RV_REF',
            'TEFF',
            'TEFF_REF',
            'ALIASES',
            'CHECKED',
            'FOUND',
            'MANUAL_EDIT',
            'COMMENTS']
ID_COLS = ['GAIADR2ID', 'TWOMASSID']
VAL_COLS = ['RV', 'TEFF']
NA_COLS = ID_COLS + VAL_COLS
EXOFOP_URL = 'https://exofop.ipac.caltech.edu/'
GAIA_2MASS_URL = 'https://gea.esac.esa.int/tap-server/tap'
TOI_URL = EXOFOP_URL+'tess/download_toi.php?sort=toi&output=csv'


def hdr_to_dict(adict, kwd, ahdr):
    """
    Add value from fits header to list in dict with same keyword.

    Args:
        adict   (dict): dictionnary of list where values are added
        kwd      (str): key corresponding to value
        ahdr  (Header): fits header to fetch value from
    Returns:
        Nothing returned, adict is edited directly

    Notes:
        - If key is not in header, nan is added.
        - ahdr could in fact be any dict-like object.
    """
    try:
        adict[kwd].append(ahdr[kwd])
    except KeyError:
        adict[kwd].append(np.nan)


def get_full_sheet(sh, ws=0, index=0):
    """
    Load a full Google Sheet with Gspread pandas
    Args:
        sh (gspread_pandas.Spread): Spreadsheet object
        ws (int or str): Name or index of the worksheet (tab)
        index (int): Column number of index in the sheet
                     Default: 0, no index in the sheet
    Returns:
        df (pandas.Dataframe): dataframe with sheet content
    """

    # Open the worksheet (tab)
    sh.open_sheet(ws)

    if sh.sheet.title in ['Info', 'Need Sorting']:
        # Make sure true/false handled properly in Info sheet
        try:
            df = sh.sheet_to_df(index=index, unformatted_columns=COLNAMES)

            # Make sure nan for empty ids
            df[NA_COLS] = df[NA_COLS].replace('', np.nan)
        except IndexError:
            df = sh.sheet_to_df(index=index)

    else:
        try:
            df = sh.sheet_to_df(index=index, unformatted_columns=['OBJECT'])
        except (KeyError, AttributeError, IndexError):
            df = sh.sheet_to_df(index=index)

    return df


def update_full_sheet(sh, df, ws=0, index=0):
    """
    Load a google sheet and overwrite all of its content with a new dataframe
    using gspread_pandas.
    Args:
        sh (gspread_pandas.Spread): Spreadsheet object
        df (pandas.Dataframe): dataframe to push to sheet
        ws (int or str): Name or index of the worksheet (tab)
        index (int): Column number of index in the sheet
                     Default: 0, no index in the sheet

    Returns:
        old_df (pandas.Dataframe): dataframe of old sheet content.
    """

    # Open the worksheet
    sh.open_sheet(ws)

    old_df = sh.sheet_to_df(index=index)
    sh.df_to_sheet(df, index=index, replace=True)

    return old_df


def strip_spectral_type(series, return_mask=False):
    """
    Strip spectral type from series of string

    Args:
        series (pd.Series): series of object names (strings)
        return_mask (bool): returns boolean mask True where there is a type
    Returns:
        no_type (pd.Series): series without spectral types
        type_mask (pd.Series): boolean mask where type is given
    """
    type_mask = series.str.match('\\([OBAFGKM]\\)')
    no_type = series.copy()
    no_type[type_mask] = series[type_mask].str.slice(start=4)

    return (no_type, type_mask) if return_mask else no_type


def get_gaia_simbad(names, verbose=False):
    """
    Get Gaia DR2 ID from Simbad for a list of names
    Args:
        names (pd.Series): series of names to search in simbad
    Returns:
        gaia_id (pd.Series): series of gaia IDs with nan where no match
    """
    # Silent simbad warnings (thrown everytime no match)
    id_list = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for name in names:
            res = Simbad.query_objectids(name)
            count = 0
            try:
                for a in res['ID']:
                    if 'Gaia DR2' in a:
                        gaia_id = a.split(' ')[-1]
                        count += 1
            except TypeError:
                # If no match, keep count at 0
                pass
            if count == 0:
                if verbose:
                    print(f'No GAIA ID found for {name}')
                id_list.append(np.nan)
            elif count == 1:
                if verbose:
                    print(f'Found GAIA ID for {name}')
                id_list.append(gaia_id)
            elif count > 1:
                if verbose:
                    print(f'Multiple Gaia ID for {name}, returning no match')
                id_list.append(np.nan)
    gaia_id = pd.Series(id_list, name='GAIADR2ID', index=names.index)

    return gaia_id


def sheet_append(df_new, df_sheet):
    """
    Append new dataframe to sheet dataframe
    """
    df_sheet = df_sheet.copy()
    df_new = df_new.copy()

    df_sheet = df_sheet.append(df_new, ignore_index=True)
    df_sheet = unique_sheet(df_sheet)
    df_sheet = sort_sheet(df_sheet)

    return df_sheet


def unique_sheet(df_sheet):
    """
    Ensure no duplicates in objects of sheet
    """
    # OBJECT duplicates
    funs = dict(zip(df_sheet.columns, ['first'] * len(df_sheet.columns)))
    funs['ALIASES'] = '|'.join
    df_sheet = df_sheet.groupby('OBJECT').agg(funs).reset_index(drop=True)

    # Gaia ID duplicates (join aliases)
    id_mask = df_sheet.GAIADR2ID.notnull()
    sheet_no_id = df_sheet[~id_mask]
    df_sheet = df_sheet[id_mask]
    df_sheet = df_sheet.groupby('GAIADR2ID').agg(funs).reset_index(drop=True)
    df_sheet = df_sheet.append(sheet_no_id, ignore_index=True)

    # 2MASS ID duplicates
    id_mask = df_sheet['TWOMASSID'].notnull()
    sheet_no_id = df_sheet[~id_mask]
    df_sheet = df_sheet[id_mask]
    df_sheet = df_sheet.groupby('TWOMASSID').agg(funs).reset_index(drop=True)
    df_sheet = df_sheet.append(sheet_no_id, ignore_index=True)

    # Update aliases
    df_sheet = update_aliases(df_sheet)

    return df_sheet


def sort_sheet(df_sheet):
    """
    Sort sheet by object name.
    """
    df_sheet = df_sheet.copy()

    # Sort new dataframe
    # old_obj = df_sheet.OBJECT.copy()
    # df_sheet.OBJECT = old_obj.str.upper()
    # df_sheet = df_sheet.sort_values('OBJECT')
    # df_sheet.OBJECT = old_obj  # Will properly map index
    # df_sheet = df_sheet.reset_index(drop=True)

    # Pandas > 1.1 simpler way
    df_sheet = df_sheet.sort_values(
            'OBJECT',
            ignore_index=True,
            key=lambda v: v.str.upper(),
            )

    return df_sheet


def gaia_crossmatch(df):
    """
    Cross-match gaia id from coordinates in a dataframe
    """
    df = df.copy()

    # Some APERO database args
    params = ao.FAKE_PARAMS
    pconst = ao.Pconst()
    objdbm = None
    outdict = dict()

    # Loop through targets
    for _, row in df.iterrows():
        gaia_id = None
        ra = row['RA_DEG']
        dec = row['DEC_DEG']
        objname = row['OBJECT']
        mjd = row['MJDEND']

        # Setup object for this target
        astro_obj = ao.AstroObject(params, pconst, gaia_id, ra, dec, objdbm,
                                   objname, 0.0, 0.0, 0.0, 0.0, None)

        # Resolve tgarget with gaia id
        astro_obj.resolve_target(mjd)

        # Write to dict
        astro_obj.write_table(outdict)

    df['GAIADR2ID'] = outdict['GAIADR2ID']

    return df


def get_aliases_gaia(gaia_ids):
    """
    Get Simbad aliases corresponding to gaia ids.
    Args:
        gaia_ids (pd.Series): series of gaia_ids
    Returns:
        aliases (pd.Series): series of pipe-separated aliases for each gaia id.
    """
    alias_all = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for gaia_id in gaia_ids:
            tbl = Simbad.query_objectids('Gaia dr2 {}'.format(gaia_id))
            try:
                aliases = tbl.to_pandas()['ID'].str.decode('utf-8').tolist()
                aliases = '|'.join(aliases)
            except AttributeError:
                aliases = None

            alias_all.append(aliases)

    aliases = pd.Series(alias_all, name='ALIASES', index=gaia_ids.index)

    return aliases


def get_aliases_twomass(twomass):
    """
    Get Simbad aliases corresponding to gaia ids.
    Args:
        twomass (pd.Series): series of 2MASS IDs
    Returns:
        aliases (pd.Series): series of pipe-separated aliases for each ID
    """
    alias_all = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for tm_id in twomass:
            tbl = Simbad.query_objectids('2MASS J{}'.format(tm_id))
            try:
                aliases = tbl.to_pandas()['ID'].str.decode('utf-8').tolist()
                aliases = '|'.join(aliases)
            except AttributeError:
                aliases = None

            alias_all.append(aliases)

    aliases = pd.Series(alias_all, name='ALIASES', index=twomass.index)

    return aliases


def update_aliases(df, new_aliases=None, inplace=False):
    """
    For each object:
        1. Add new aliases if any
        2. Make sure unique pipe-separated list with OBJECT first in list
    Args:
        df (pd.DataFrame): dataframe with an ALIASES column
        new_aliases (Series, array, list): aliases to append for each object.
            If Series, make sure indexes match the dataframe
        inplace (bool): If true, edit df directly, otherwise return a copy
    """
    if not inplace:
        df = df.copy()

    # Make sure all strings
    df.ALIASES = df.ALIASES.fillna('')

    # Add new aliases
    if new_aliases is not None:
        new_aliases = new_aliases.fillna('')
        df.ALIASES = df.ALIASES.str.cat(new_aliases, sep='|')

    # Make sure 2MASS and Gaia included
    aliases = df.ALIASES.copy()
    mask_gaia = df.GAIADR2ID.notnull()
    aliases[mask_gaia] = aliases[mask_gaia].str.cat(
            'Gaia DR2 ' + df.GAIADR2ID[mask_gaia]
            )
    mask_2mass = df['TWOMASSID'].notnull()
    aliases[mask_2mass] = aliases[mask_2mass].str.cat(
            '2MASS J' + df['TWOMASSID'][mask_2mass]
            )

    # Make sure object is not current list, and that list is unique
    objs = df.OBJECT.tolist()
    aliases = df.ALIASES.apply(  # Remove object name and empty strings
            lambda row: [val for val in row.split('|')
                         if val not in objs + ['']
                         ]
            )
    aliases = aliases.apply(set).str.join('|')  # Unique and join with pipe

    # Add object name as first alias
    al_mask = aliases.str.len() > 0
    df.loc[al_mask, 'ALIASES'] = df.OBJECT.str.cat(aliases, sep='|')
    df.loc[~al_mask, 'ALIASES'] = df.OBJECT.copy()

    return None if inplace else df


def check_gaia_id(df, verbose=False):
    """
    Steps:
        1. String comparison with  formatting to match more objects. Some
           catalog prefixes are remapped (e.g. Gl->GJ) to improve matches.
        2. For TOIs, which are not always on Simbad, we also check the MAST
           database to verify if Gaia ID matches TOI or TIC
    A string comparison is done with a bit of formatting to match more objects.

    Args:
        df (pd.Dataframe): dataframe with object list
    Returns:
        df with CHECKED and FOUND columns updated according to search result
    """
    df = df.copy()

    # Get aliases and format (all caps no spaces), remove OBJECT name (first)
    aliases = df.loc[~df['CHECKED']]['ALIASES']
    aliases = aliases.str.upper().str.replace(' ', '').str.split('|').str[1:]
    aliases = aliases.tolist()

    # Format names for match
    names = df.loc[~df['CHECKED']]['OBJECT'].str.upper().str.replace(' ', '')
    names = names.tolist()

    # String-based checks and update
    match = check_str(names, aliases)
    df['FOUND'] = df['FOUND'].mask(~df['CHECKED'], match)  # Update matches
    if verbose:
        print('The following objects were not found with string matching.')
        print(df['OBJECT'].loc[~df['FOUND']])
    df['CHECKED'] = df['CHECKED'].mask(~df['CHECKED'], match)  # Update checks

    # TOI validation
    df_left = df.loc[~df['CHECKED']]  # Only target not found yet
    if df_left['OBJECT'].str.startswith('TOI').any():
        df_left = check_tess(df_left)
        df.loc[~df['CHECKED']] = df_left
        if verbose:
            print('The following objects were not found with TOI validation.')
            print(df['OBJECT'].loc[~df['FOUND']])
    elif verbose:
        print('Noi TOIs to look up.')

    df_left = df.loc[~df['CHECKED']]  # Only target not found yet
    df_left = check_simbad(df_left)
    df.loc[~df['CHECKED']] = df_left
    if verbose:
        print('The following objects were not found in SIMBAD.')
        print(df['OBJECT'].loc[~df['FOUND']])

    return df


def check_str(names, aliases):
    """
    Check if names match the simbad alias list.
    Args:
        names (list of str): object names in data headers
        aliases (list of str): aliases separated by '|' for each object.
    """
    # Basic string check
    match = []
    for name, alist in zip(names, aliases):
        # Some corner cases with catalog names
        if name[:2] == 'GL':
            name = 'GJ'+name[2:]
        elif name[:3] == 'BBW':
            name = 'BRAN'+name[3:]
        match.append(any(name in a for a in alist))

    return match


def check_tess(df):
    """
    Verify TESS objects with MAST archive
    Args:
        df (pd.Dataframe): Dataframe with objects not found yet
    Returns:
        df with match info updated
    """

    df = df.copy()

    # Load TOI/TIC list
    data = requests.get(TOI_URL)
    tbl = Table.read(data.text, format='ascii')
    df_toi = tbl.to_pandas()[['TIC ID', 'TOI']].astype(str)

    # Replace TOI names with TIC
    # Keep order consistent with dataframe
    toi_inds = df.index[df['OBJECT'].str.startswith('TOI')]
    tois = df['OBJECT'].loc[toi_inds]
    tois = tois.str.split('-').map(lambda x: x[1]) + '.01'
    tois = tois.tolist()
    toi_to_tic = df_toi.loc[df_toi['TOI'].isin(tois)]
    toi_to_tic = toi_to_tic.set_index('TOI', drop=False)
    toi_to_tic = toi_to_tic.reindex(tois).set_index(np.arange(len(tois)))
    tics = toi_to_tic['TIC ID']

    gaiaids = df['GAIADR2ID'].loc[toi_inds]

    # Check where TIC ID matches same as local Gaia ID
    if tics.size == 0:
        return df
    tic_gaia = Catalogs.query_criteria(
            catalog='TIC',
            objType='STAR',
            ID=tics.tolist()).to_pandas()[['ID', 'GAIA']]
    tic_gaia = tic_gaia.set_index('ID', drop=False)
    tic_gaia = tic_gaia.reindex(tics).set_index(np.arange(len(tics)))
    mask = tic_gaia['GAIA'].values == gaiaids.values

    # Record if match or not
    df['FOUND'].loc[toi_inds] = mask
    df['CHECKED'].loc[toi_inds] = mask

    # Where no match, replace by new ID
    df['GAIADR2ID'].loc[toi_inds] = df['GAIADR2ID'].loc[toi_inds].where(
                                                mask,
                                                other=tic_gaia['GAIA'])
    msg = 'Updated Gaia ID from TIC'
    df['COMMENTS'].loc[toi_inds] = df['COMMENTS'].loc[toi_inds].where(
                                                mask,
                                                other=msg)
    df['FOUND'].loc[toi_inds] = True
    df['CHECKED'].loc[toi_inds] = True

    return df


def check_simbad(df, verbose=False):
    """
    Check list of names in simbad and return gaia id if any.
    """
    df = df.copy()
    names = df['OBJECT']
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for i, name in enumerate(names):
            res = Simbad.query_objectids(name)
            count = 0
            try:
                for a in res['ID']:
                    if 'Gaia DR2' in a:
                        gaiaid = a.split(' ')[-1]
                        count += 1
            except TypeError:
                pass
            if count == 0:
                if verbose:
                    print('No GAIA ID found for {}'.format(name))
            elif count == 1:
                if verbose:
                    print('Found GAIA ID for {}'.format(name))
                df['GAIADR2ID'].iloc[i] = gaiaid
                df['COMMENTS'].iloc[i] = 'Updated Gaia ID from SIMBAD'
                df['FOUND'].iloc[i] = True
            elif count > 1:
                if verbose:
                    print('Multiple Gaia ID for {},'
                          ' returning no match'.format(name)
                          )

    return df


def vals_and_counts(df, kwd, threshold=None):
    """
    Get values and their counts in two separate arrays
    """

    vals = df[kwd]
    counts = vals.value_counts()
    vals = vals.dropna().unique()

    # reject unphysical values
    if threshold is not None:
        mask = np.logical_and(
                np.abs(vals) < threshold[1],
                np.abs(vals) > threshold[0],
                )
        vals = vals[mask]

    vals = vals[np.nonzero(vals)]
    counts = counts[vals].values

    return vals, counts


def pipelist(arr):
    """
    Return string with pipe-separated list
    """
    arr = np.array(arr)
    return '|'.join(list(arr.astype(str)))


def tryaddlast(alist, arr):
    """
    Try appending last index of array to list, if empty add nan
    """

    try:
        alist.append(arr[-1])
    except IndexError:
        alist.append(np.nan)


def twomass_from_gaia(df):
    """
    Get 2MASS IDs from Gaia DR2 ID and add in column of dataframe
    Args:
        df (pd.Dataframe): dataframe with object list
    Returns:
        df with updated 2mass IDs
    """
    gaia_ids = df['GAIADR2ID'].astype(str).tolist()
    # gaia_ids_end = df['GAIADR2ID'].astype(str).tolist()[100:]
    # gaia_ids.extend(gaia_ids_end)
    gaia_ids = '(' + ', '.join(gaia_ids) + ')'
    query = (
            'SELECT source_id, original_ext_source_id '
            'FROM gaiadr2.tmass_best_neighbour '
            'WHERE source_id IN {}'
            ).format(gaia_ids)
    gaia = TapPlus(url=GAIA_2MASS_URL)
    job = gaia.launch_job(query=query)
    tbl = job.get_results()
    df_2mass = tbl.to_pandas()

    df_2mass = df_2mass.astype(str).set_index('source_id', drop=False)
    df_2mass = df_2mass.reindex(df['GAIADR2ID'].astype(str).tolist())
    df_2mass = df_2mass['original_ext_source_id']

    # Weird encoding workaround
    df['TWOMASSID'] = df_2mass.str.strip("b\'\"").values

    return df


def twomass_from_simbad(df):
    """
    Check simbad aliases for missing 2mass ids
    """
    print(df.index[df['2MASSID'].isna()])
    print(df.ALIASES)
    for i in df.index[df['2MASSID'].isna()]:
        aliases = df['ALIASES'].loc[i].split('|')
        count = 0
        for alias in aliases:
            if '2MASS' in alias:
                twomass = alias
                count += 1
        if count == 1:
            df['2MASSID'].iloc[i] = twomass.split('J')[1]
            df['COMMENTS'].iloc[i] = (df['COMMENTS'].iloc[i]
                                      + ' 2MASS ID from SIMBAD'
                                      )

    return df
