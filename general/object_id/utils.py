"""
Some utility functions to play with astro databases
"""
import warnings

import requests
import numpy as np
from astropy.table import Table
from astroquery.simbad import Simbad
from astroquery.mast import Catalogs
from astroquery.utils.tap.core import TapPlus
import gspread_pandas as gspd

import berv_crossmatch as bcm

EXOFOP_URL = 'https://exofop.ipac.caltech.edu/'
TOI_URL = EXOFOP_URL+'tess/download_toi.php?sort=toi&output=csv'

# GAIA_2MASS_URL = 'https://gea.esac.esa.int/tap-server/tap'
GAIA_2MASS_URL = 'https://gaia.obspm.fr/tap-server/tap'

COLNAMES = ['OBJECT',
            'GAIADR2ID',
            '2MASSID',
            'RV',
            'RV_REF',
            'TEFF',
            'TEFF_REF',
            'ALIASES',
            'CHECKED',
            'FOUND',
            'MANUAL_EDIT',
            'COMMENTS']


def get_aliases_from_gaia(gaiaid):
    """Get Simbad Aliases for a Gaia ID.
    Args:
        GaiaID (str or int): Gaia ID to look up.
    Returns:
        Aliases (str): All names separated by '|'.
    """

    tbl = Simbad.query_objectids(f'Gaia dr2 {gaiaid}')
    try:
        aliases = '|'.join(tbl.to_pandas()['ID'].str.decode('utf-8').tolist())
    except AttributeError:
        aliases = None

    return aliases


def add_aliases(df):
    """
    Add alias column to dataframe of object based on Gaia ID
    Args:
        df (pd.Dataframe): dataframe with object list
    Returns:
        df with additional aliases column
    Note: There must be a GAIADR2ID column in the dataframe
    """

    for index, row in df.iterrows():
        if type(row['ALIASES']) == str:
            continue
        elif np.isnan(row['ALIASES']):
            simbad_ids = get_aliases_from_gaia(row['GAIADR2ID'])
            if simbad_ids is None:
                continue
            aliases = '|'.join([row['OBJECT'], simbad_ids])
            df.at[index, 'ALIASES'] = aliases
        else:
            raise ValueError('Unexpected alias type')

    return df


def check_id(df, replace=False):
    """
    Steps:
        1. String comparison with  formatting to match more objects. Some
           catalog prefixes are remapped (e.g. Gl->GJ) to improve matches.
        2. For TOIs, which are not always on Simbad, we also check the MAST
           database to verify if Gaia ID matches TOI or TIC
    A string comparison is done with a bit of formatting to match more objects.

    Args:
        df (pd.Dataframe): dataframe with object list
        replace (bool): replace Gaia ID with results from catalogs if found
    Returns:
        df with CHECKED and FOUND columns updated according to search result
    """

    df = df.copy()

    # Get alias list (without first, which is just object copied direclty)
    aliases = df.loc[~df['CHECKED']]['ALIASES']
    aliases = aliases.str.upper().str.replace(' ', '').str.split('|').str[1:]
    aliases = aliases.tolist()

    # Format names
    names = df.loc[~df['CHECKED']]['OBJECT'].str.upper().str.replace(' ', '')
    names = names.tolist()

    # String-based checks only
    match = check_str(names, aliases)
    df['FOUND'] = df['FOUND'].mask(~df['CHECKED'], match)  # Update matches
    print('The following objects were not found with string matching.')
    print(df['OBJECT'].loc[~df['FOUND']])

    # Update CHECKED only for objects found so far
    df['CHECKED'] = df['CHECKED'].mask(~df['CHECKED'], match)  # Update checks

    # TOI validation
    df_left = df.loc[~df['CHECKED']]  # Only target not found yet
    df_left = check_tess(df_left, replace=replace)
    df.loc[~df['CHECKED']] = df_left
    print('The following objects were not found with TOI validation.')
    print(df['OBJECT'].loc[~df['FOUND']])

    if replace:
        df_left = df.loc[~df['CHECKED']]  # Only target not found yet
        check_simbad(df_left)
        df.loc[~df['CHECKED']] = df_left
    print('The following objects were not found in SIMBAD (raw string only).')
    print(df['OBJECT'].loc[~df['FOUND']])

    # Mark all as checked
    df['CHECKED'] = df['CHECKED'].mask(~df['CHECKED'], True)

    return df


def check_str(names, aliases):
    """
    Check that names match the simbad alias list.
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


def check_tess(df, replace=False):
    """
    Verify TESS objects with MAST archive
    Args:
        df (pd.Dataframe): Dataframe with objects not found yet
        replace (bool): whether we should overwrite Gaia IDs of TICs with
            different match in catalog
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

    if replace:
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
                    print(f'No GAIA ID found for {name}')
            elif count == 1:
                if verbose:
                    print(f'Found GAIA ID for {name}')
                df['GAIADR2ID'].iloc[i] = gaiaid
                df['COMMENTS'].iloc[i] = 'Updated Gaia ID from SIMBAD'
                df['FOUND'].iloc[i] = True
            elif count > 1:
                if verbose:
                    print(f'Multiple Gaia ID for {name}, returning no match')

    return df


def update_full_sheet(sheet_id, df, ws=0):
    """
    Load a google sheet and overwrite all of its content with a new dataframe
    using gspread_pandas.
    Args:
        sheet (str): ID or URL of the google sheet
        df (pandas.Dataframe): dataframe to push to sheet
        ws (int or str): Name or index of the worksheet (tab)

    Returns:
        old_df (pandas.Dataframe): dataframe of old sheet content.
    """

    sh = gspd.spread.Spread(sheet_id)
    sh.open_sheet(ws)

    old_df = sh.sheet_to_df(index=0)  # Assume normal user so to indices
    sh.df_to_sheet(df, index=False, replace=True)

    return old_df


def get_full_sheet(sheet_id, ws=0):
    """
    Load a full Google Sheet with Gspread pandas
    Args:
        sheet_id (str): ID or URL of the google sheet
        ws (int or str): Name or index of the worksheet (tab)

    Returns:
        df (pandas.Dataframe): dataframe with sheet content
    """

    sh = gspd.spread.Spread(sheet_id)
    sh.open_sheet(ws)

    if sh.sheet.title in ['Info', 'Need Sorting']:
        # Make sure true/false handled properly in Info sheet
        try:
            df = sh.sheet_to_df(index=0, unformatted_columns=COLNAMES)
        except IndexError:
            df = sh.sheet_to_df(index=0)
    else:
        try:
            df = sh.sheet_to_df(index=0, unformatted_columns=['OBJECT'])
        except (KeyError, AttributeError, IndexError):
            df = sh.sheet_to_df(index=0)

    return df


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
    df['2MASSID'] = df_2mass.str.strip("b\'\"").values

    return df


def twomass_from_simbad(df):
    """
    Check simbad aliases for missing 2mass ids
    """
    for i in df.index[df['2MASSID'].isna()]:
        aliases = df['ALIASES'].iloc[i].split('|')
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


def gaia_crossmatch(df):
    """
    Cross-match gaia id from coordinates in a dataframe
    """
    df = df.copy()

    # Some APERO database args
    params = bcm.FAKE_PARAMS
    pconst = bcm.Pconst()
    objdbm = None
    outdict = dict()

    # Loop through targets
    for i, row in df.iterrows():
        gaia_id = None
        ra = row['RA_DEG']
        dec = row['DEC_DEG']
        objname = row['OBJECT']
        mjd = row['MJDEND']

        # Setup object for this target
        astro_obj = bcm.AstroObject(params, pconst, gaia_id, ra, dec, objdbm,
                                    objname, 0.0, 0.0, 0.0, 0.0, None)

        # Resolve tgarget with gaia id
        astro_obj.resolve_target(mjd)

        # Write to dict
        astro_obj.write_table(outdict)

    df['GAIADR2ID'] = outdict['GAIADR2ID']

    return df
