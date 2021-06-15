"""
Functions to facilitate maintenance of apero sheet.
Primarly designed for automated updates with update_sheet.py, but also useful
for interactive editing.

@author: vandalt
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import tqdm
from astropy.io import fits
from astropy.io.votable.tree import Table

import utils as ut

# =============================================================================
# Constants
# =============================================================================
# Pattern to search for raw files on disk
RAW_PATTERN = "/spirou/cfht_nights/common/raw/20*/*o.fits"
# Local file that stores info about previously parsed raw files
LOCAL_FILE = "object_info.csv"
# Directory with dfits outputs allowing to include objects from other databases
DFITS_DIR = "dfits_outputs"
OBJ_INFO_KEYS = [
    "OBJECT",
    "OBJRV",
    "OBJTEMP",
    "MJDEND",
    "RA_DEG",
    "DEC_DEG",
]

SHEET_ID = "1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8"

BAD_NAMES = [
    "sky.*_.*",
    ".*_sky.*",
    "sky",
    "SKY",
    "FakeTarget",
    "Test.*",
    ".*Test",
    "Engineering",
    "Calibration",
]
# Known rejected names
REJ_NAMES = [
    "Moon.*",
    "Neptune.*",
    "Saturn.*",
    "Venus.*",
    "17BE93-FT-1",
]


def get_object_info(
    fpattern, local_file=None, dfits_dir=None, keys=None, verbose=False
):
    """
    Fetch object information in headers of raw data from CFHT.
    If a local file is found, files already read will be removed for speed.
    Local txt files with dfits outputs can also be used to add objects from
    external sources
    Args:
        fpattern (str): Glob-compatible pattern for raw files
        local_file (str): Local file with observations already fetched.
                          (if this code has been run in the past)
        dfits_dir (str): directory with dfits output as txt files.
        keys    (list): Header keys to extract.
                Defaults to [OBJECT, OBJRV, OBJTEMP, MJDEND, RA_DEG, DEC_DEG].
        verbose   (bool): Additional information printed if True.
    """
    if keys is None:
        keys = OBJ_INFO_KEYS

    if local_file is not None:
        if os.path.isfile(local_file):
            df_loc = pd.read_csv(local_file)
            loc_files = df_loc.FILE.tolist()
        else:
            if verbose:
                print("No local file named {}".format(local_file))
            loc_files = []
    else:
        loc_files = []

    if dfits_dir is not None:
        dfits_frames = []
        for dfile in glob.glob(os.path.join(dfits_dir, "*.txt")):
            dfits_frames.append(
                Table.read(dfile, format="ascii.tab").to_pandas(use_nullable_int=False)
            )
        df_dfits = pd.concat(dfits_frames, ignore_index=True)
        df_dfits.FILE = df_dfits.FILE.apply(os.path.basename)
        df_dfits = df_dfits[["FILE"] + keys]
    else:
        df_dfits = pd.DataFrame([])

    # Find files to read
    full_list = glob.glob(fpattern)
    all_files = [os.path.basename(f) for f in full_list]
    loc_mask = np.isin(all_files, loc_files)
    iter_list = np.array(full_list)[~loc_mask].tolist()

    # Initialize dict
    outdict = dict()
    for k in keys:
        outdict[k] = []

    # Handle files separately because not inside header
    fnames = []

    # Loop through files
    # Slower than dfits | fitsort but avoids system calls
    if verbose:
        print("Fetching data from headers. This might take a few minutes.")
    for filepath in tqdm.tqdm(iter_list):
        hdr = fits.getheader(filepath)
        for k in keys:
            ut.hdr_to_dict(outdict, k, hdr)
        fnames.append(os.path.basename(filepath))

    # Create dataframe of new observationsfrom dict
    df = pd.DataFrame([])
    df["FILE"] = fnames
    for k in keys:
        df[k] = outdict[k]

    # Append to local df if any
    if len(loc_files) > 0:
        df = df_loc.append(df, ignore_index=True)
    if len(df_dfits) > 0:
        df = df.append(df_dfits, ignore_index=True)
    df = df.sort_values("FILE")
    df = df.drop_duplicates("FILE")  # Could have duplicates from externalf fits

    return df


def bad_and_rejected(names, sh=None, bad=None, rej=None, verbose=False):
    """
    Handle bad (calibration, sky, etc.) and rejected names.
    Rejected names are objects that we do not wish to include in the list for
    various reasons (e.g. solar system objects). They are stored in a separate
    tab of the Google Sheet.
    Args:
        names (pd.Series): series of names to update
        sh (gspread_pandas.Spread): Spreadsheet object to update with rejected
                                    names.
                                    Default: None, no remote update
        bad (list): list of names to completely remove, by default uses a
                    pre defined list.
        rej (list): list of rejected object names, by default uses a
                         pre defined list.
    Returns:
        good_names (pd.Series): Series of good names.
    """
    # Pre-process lists into regex list
    if bad is None:
        bad = BAD_NAMES
    if rej is None:
        rej = REJ_NAMES
    bad = "|".join(bad)  # pipe symbol as OR in regex
    rej = "|".join(rej)  # pipe symbol as OR in regex

    # Remove bad
    good_names = names[~names.str.fullmatch(bad)]

    # Remove rejected
    rej_names = good_names[good_names.str.fullmatch(rej)]

    # Update rejected names in google sheet
    if sh is not None:
        rej_names = pd.DataFrame(rej_names)
        df_rej = ut.get_full_sheet(sh, ws="Rejected")  # Remote names
        rej_all = "|".join(df_rej.OBJECT.tolist())
        rej_all = "|".join([rej, rej_all])
        rej_mask = rej_names.OBJECT.isin(df_rej.OBJECT)
        if not rej_mask.all():  # Update only  if new objects
            if verbose:
                print("Updating rejected sheet with new objects")
            df_rej = df_rej.append(rej_names[~rej_mask], ignore_index=True)
            df_rej = df_rej.sort_values("OBJECT")
            ut.update_full_sheet(sh, df_rej, ws="Rejected")
    else:
        rej_all = rej

    good_names = good_names[~good_names.str.fullmatch(rej)]

    return good_names


def check_gaia(names, df_sheet):
    """
    Check if names are duplicates from df_sheet based on Gaia DR2 ID.
    This does three things:
        1. Add matching IDs as aliases
        2. Add non-matched IDs as new object
        3. Keep objects without IDs for testing
    Args:
        names (pd.Series): series of object names
        df_sheet (pd.DataFrame): dataframe from Object ID Google Sheet
    Returns:
        not_found (pd.Series): names without ID
        df_sheet (pd.DataFrame): updated Google Sheet dataframe
    """
    # Don't edit reference to external array
    df_sheet = df_sheet.copy()

    # Remove spectral type prefix if there is one
    no_type_names = ut.strip_spectral_type(names)

    # Get gaia ID
    gaia_id = ut.get_gaia_simbad(no_type_names)

    # If matching ID in (non-null) sheet: duplicate, so add as name as alias
    dup_gaia_id = gaia_id[gaia_id.isin(df_sheet.GAIADR2ID.dropna())]
    dup_names = names[gaia_id.isin(df_sheet.GAIADR2ID.dropna())]
    inds = np.where(dup_gaia_id.values[:, None] == df_sheet.GAIADR2ID.values)[-1]
    new_al = df_sheet.loc[inds].ALIASES.str.cat(dup_names, sep="|")
    df_sheet.loc[inds, "ALIASES"] = new_al

    # If new ID in (non-null) sheet: new object, append to dataframe
    new_gaia_id = gaia_id[~gaia_id.isin(df_sheet.GAIADR2ID.dropna())]
    new_names = names[~gaia_id.isin(df_sheet.GAIADR2ID.dropna())]
    df_new = pd.DataFrame([new_names, new_gaia_id]).T
    df_new = df_new.dropna(subset=["GAIADR2ID"])  # Drop nan IDs
    df_new.MANUAL_EDIT = False
    if df_new.shape[0] > 0:
        # Append new and sort sheet (function does both)
        df_sheet = ut.sheet_append(df_new, df_sheet)

    # If no ID (nan) keep testing
    null_names = names[gaia_id.isnull()]

    return null_names, df_sheet


def gaia_match_position(df_pos, df_test, df_sheet):
    """
    Match Gaia ID of test objects, verify this match based on aliases found,
    and add valid matches to the sheet.
    Objects without ID or with failed match found objects are returned in a
    separate dataframe.

    Args:
        df_pos (pd.DataFrame): dataframe with position
        df_test (pd.DataFrame): dataframe with test objects
        df_sheet (pd.DataFrame): dataframe with objects in sheet
    Returns:
        df_test (pd.DataFrame): updated test dataframe with objects not found
        df_sheet (pd.DataFrame): updated Google Sheet dataframe
    """
    # Add GAIADR2ID by position cross-match
    df_pos = ut.gaia_crossmatch(df_pos)
    ord_ids = df_pos.set_index("OBJECT").loc[df_test.OBJECT].GAIADR2ID
    df_test["GAIADR2ID"] = ord_ids.values  # sorted above so safe to use values

    # Get aliases for each ID, make sure no nans
    aliases = ut.get_aliases_gaia(df_test.GAIADR2ID)

    # Update dataframe aliases with new ones
    df_test = ut.update_aliases(df_test, new_aliases=aliases)

    # Check Gaia ID and split test array
    df_test = ut.check_gaia_id(df_test)
    df_new = df_test[df_test.FOUND].copy()
    df_new.MANUAL_EDIT = False
    df_test = df_test[~df_test.FOUND].copy()

    # Append found objects to dataframe
    if df_new.shape[0] > 0:
        # Append new and sort sheet (function does both)
        df_sheet = ut.sheet_append(df_new, df_sheet)

    # Mark test dataframe as unchecked
    df_test.CHECKED = False

    return df_test, df_sheet


def str_match(df_test, df_sheet, min_len=2):
    """
    Find potential matches by string comparison. Since string comparison is
    prone to false positives, no updates are done in the Info sheet. Possible
    matches will appear as comments in the Need Sorting sheet.
    Args:
        df_test (pd.DataFrame): dataframe with test objects
        df_sheet (pd.DataFrame): dataframe with objects in sheet
    """
    df_test = df_test.copy()

    df = df_sheet.set_index("OBJECT")

    # Setup test names
    names = df_test.OBJECT
    test_names = ut.strip_spectral_type(names)
    test_names = names.str.upper().str.replace(" ", "").replace("_", "")
    test_names = test_names.str.replace("*", "")

    # Check if exact match
    matches_full = test_names.apply(
        lambda tn: [obj for obj in df.index if np.isin(tn, df.ALIASES[obj].split("|"))]
    )
    mask = matches_full.str.len() > 0
    comments = pd.Series([""] * test_names.size, index=test_names.index)
    comments[mask] += "Full: " + matches_full[mask].str.join(",")

    # Check if test name contains an alias
    esc_aliases = df.ALIASES.apply(
        lambda s: "|".join(re.escape(sub) for sub in s.split("|"))
    )
    matches_alias_in_test = test_names.apply(
        lambda tn: [obj for obj in df.index if bool(re.search(esc_aliases[obj], tn))]
    )
    mask = matches_alias_in_test.str.len() > 0
    comments[mask] += "  Alias in test: " + matches_alias_in_test[mask].str.join(",")

    # Check if an alias contains a test name
    matches_test_in_alias = test_names.apply(
        lambda tn: [
            obj for obj in df.index if tn in df.ALIASES[obj] and len(tn) > min_len
        ]
    )
    mask = matches_test_in_alias.str.len() > 0
    comments[mask] += "  Test in alias: " + matches_test_in_alias[mask].str.join(",")

    comments = comments.replace("", np.nan)
    print(comments)

    # Enter possible matches as comments
    df_test["COMM_NEW"] = comments
    print(df_test)
    comm_mask = df_test["COMM_NEW"].notnull()
    df_test.loc[comm_mask, "COMMENTS"] = df_test["COMM_NEW"][comm_mask]
    df_test = df_test.drop("COMM_NEW", axis=1)

    return df_test


def get_2mass_id(df_sheet):
    """
    Search 2MASS ID of objects
    """
    # Avoid crashing gaia query
    mask = df_sheet["TWOMASSID"].isnull() & df_sheet.GAIADR2ID.notnull()

    df_sheet[mask] = ut.twomass_from_gaia(df_sheet[mask])

    mask = df_sheet["TWOMASSID"].isnull()

    df_sheet[mask] = ut.twomass_from_simbad(df_sheet[mask])

    return df_sheet


def find_aliases_simbad(df_sheet):
    """
    Get simbad aliases for a sheet
    """

    # Get aliases from Gaia ID
    mask = df_sheet.ALIASES.str.split("|").str.len() < 4
    df_update = df_sheet[mask]
    aliases = ut.get_aliases_gaia(df_update.GAIADR2ID)
    df_update = ut.update_aliases(df_update, new_aliases=aliases)
    df_sheet[mask] = df_update

    # Get aliases from 2MASS ID
    mask = df_sheet.ALIASES.str.split("|").str.len() < 4
    df_update = df_sheet[mask]
    aliases = ut.get_aliases_twomass(df_update["TWOMASSID"])
    df_update = ut.update_aliases(df_update, new_aliases=aliases)
    df_sheet[mask] = df_update

    return df_sheet


def teff_and_rv(sh, df_sheet, df_loc):
    """
    Update TEFF and RV information.
    """

    # Initialize lists
    rv_list = []
    teff_list = []
    rv_all = []
    teff_all = []
    rv_all_count = []
    teff_all_count = []

    for remote_object in df_sheet["OBJECT"]:
        objrows = df_loc[df_loc["OBJECT"] == remote_object]
        rv, rv_count = ut.vals_and_counts(objrows, "OBJRV", threshold=(0, 2000))
        teff, teff_count = ut.vals_and_counts(
            objrows,
            "OBJTEMP",
            threshold=(1000, np.inf),
        )

        # Append latest value to "main" list
        ut.tryaddlast(rv_list, rv)
        ut.tryaddlast(teff_list, teff)

        # Keep pipe-separated string of all values and counts
        rv_all.append(ut.pipelist(rv))
        rv_all_count.append(ut.pipelist(rv_count))
        teff_all.append(ut.pipelist(teff))
        teff_all_count.append(ut.pipelist(teff_count))

    # Update RV and TEFF in Info worksheet
    update_mask = ~df_sheet["MANUAL_EDIT"]
    rv_teff_vals = np.array([rv_list, teff_list]).T
    rv_teff_vals = pd.DataFrame(rv_teff_vals, columns=["RV", "TEFF"])
    update_mask = update_mask.values[:, None] & rv_teff_vals.notnull()
    rv_teff_update = rv_teff_vals[update_mask]
    for col in ["RV", "TEFF"]:
        rv_teff_update[col + "_REF"] = "HEADER[CFHT]"  # Add reference too
        col2 = [col, col + "_REF"]
        df_sheet.loc[update_mask[col], col2] = rv_teff_update[col2]
    ut.update_full_sheet(sh, df_sheet, ws="Info")

    # Update the maintenance sheet (no mask for this one)
    df_maint = ut.get_full_sheet(sh, ws="Maintenance", index=1)
    df_maint = df_maint.reindex(df_sheet.OBJECT)
    df_maint = df_maint.fillna("")
    new_arr = np.array([rv_all, rv_all_count, teff_all, teff_all_count]).T
    df_new = pd.DataFrame(
        new_arr,
        columns=df_maint.columns,
        index=df_maint.index,
    )
    for col in df_maint.columns[1:]:
        try:
            df_maint[col] = df_maint[col].str.cat(df_new[col].astype(str), sep="|")
            df_maint[col] = (
                df_maint[col]
                .str.split("|")
                .apply(lambda vals: [str(float(val)) for val in vals if val != ""])
                .apply(pd.unique)
                .str.join("|")
            )
        except AttributeError:
            df_maint[col] = df_new[col]
    ut.update_full_sheet(sh, df_maint, ws="Maintenance", index=1)
