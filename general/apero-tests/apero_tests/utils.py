"""
General functions to use in apero tests
"""
import glob
import os
import warnings
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from apero.core import constants
from apero.core.constants import param_functions
from apero.core.constants.param_functions import ParamDict
from astropy.io import fits
from astropy.table import Table
from pandas import DataFrame, Series


def removext(name: str, ext: str = ".py") -> str:
    """
    Remove extension of a file/recipe name

    :param name: file/recipe name
    :type name: str
    :param ext: filetype extension, default is ".py"
    :type ext: str :returns: cleaned up name
    :rtype:  str
    """
    if not ext.startswith("."):
        ext = "." + ext

    while name.endswith(ext):
        name = name[:-3]

    return name


def get_nth_parent(path: str, order: int = 1):
    """
    Get nth order parent of path by applying os.dirname repeatedly.

    :param path: Input path
    :type path: str
    :param order: Number of directories up (number of times os.dirname is
                  applied).
                  Default is 1.
    :type order: int
    """
    parent_path = path
    for _ in range(order):
        parent_path = os.path.dirname(parent_path)

    return parent_path


def load_db_entries(db_id: str, instrument: str = "SPIROU") -> DataFrame:
    """
    Load an APERO database entries using the database 'ID'

    :param db_id: Database ID (e.g. 'CALIB' for calibdb)
    :type db_id: str
    :param instrument: Instrument to use when loading apero params
    :type instrument: str
    :return: Dataframe containing the database
    :rtype: DataFrame
    """
    # NOTE: For 0.7 version
    # from apero.core.core import drs_database
    # from apero.core import constants
    # params = constants.load()
    # calibdb = drs_database.CalibDatabase(params)
    # calibdb.load_db()
    # calib_table = calibdb.database.get('*', calibdb.database.tname,
    #                                    return_pandas=True)

    # TODO: Should we have check for db_id, does APERO have list of db names ?
    db_id = db_id.upper()

    params = constants.load(instrument)
    db_path = os.path.join(
        params[f"DRS_{db_id}_DB"], params[f"{db_id}_DB_NAME"]
    )

    colnames = params.listp(f"{db_id}_DB_COLS", dtype=str)
    db_arr = np.loadtxt(db_path, dtype=str, unpack=True)
    df = pd.DataFrame(dict(zip(colnames, db_arr)))

    # TODO: When sure that apero uses pandas >=1.0, use convert_dtypes here

    return df


def load_db_list(db_id: str, instrument: str = "SPIROU") -> DataFrame:
    """
    Load an APERO database using the database 'ID'

    :param db_id: Database ID (e.g. 'CALIB' for calibdb)
    :type db_id: str
    :param instrument: Instrument to use when loading apero params
    :type instrument: str
    :return: Dataframe containing the database
    :rtype: DataFrame
    """
    # TODO: Should we have check for db_id, does APERO have list of db names ?
    db_id = db_id.upper()

    params = constants.load(instrument)
    db_path = Path(params[f"DRS_{db_id}_DB"])

    # NOTE: Might not include only data files,
    # but need to get files independently of DB entries, so cannot just cross-check
    plist = list(db_path.glob("*.fits"))
    return [str(p.name) for p in plist]


def load_fits_df(pathlist: List[str]) -> DataFrame:
    """
    Load contents from list of fits tables in a dataframe, making sure to
    decode strings

    :param pathlist:
    :type pathlist: List[str]

    :returns:
    :rtype: DataFrame
    """

    # Pandas does not load fits files, so we convert using astropy tables
    paths = [p for p in pathlist if os.path.isfile(p)]
    df_list = [Table.read(f).to_pandas() for f in paths]
    df = pd.concat(df_list)

    # Decode astropy are encoded, this is different from pandas default
    # Pandas store string as 'object' so we keep only those columns
    for col in df.select_dtypes(include=object).columns:
        # Pandas also stores list as 'object', so check it's a 'bytes' type
        if isinstance(df[col].iloc[0], bytes):
            df[col] = df[col].str.decode("utf-8")

    return df


def make_full_index(
    real_index: DataFrame, missing_index: DataFrame
) -> DataFrame:

    real_index = real_index.copy()
    missing_index = missing_index.copy()

    # NOTE: Maybe index_key should be moved in more accessible scope
    index_key = "IN_INDEX"
    real_index[index_key] = True

    if len(missing_index) > 0:
        # NOTE: Need to add index_key here to avoid having len > 0
        missing_index[index_key] = False
        full_index = real_index.append(missing_index)
    else:
        full_index = real_index.copy()
    full_index = full_index.sort_values(["NIGHTNAME", "LAST_MODIFIED"])

    return full_index


def get_cdb_df(
    index_df: DataFrame,
    params: ParamDict,
    force: bool = False,
    cache_dir: str = None,
) -> DataFrame:
    """
    Get dataframe with all CDB files for each file in index. Index is file name
    """
    if cache_dir is not None:
        p = Path(cache_dir)
        p.mkdir(parents=True, exist_ok=True)
        cache_name = f"cdb_mjd_df_{params['INSTRUMENT']}.csv"
        cache_path = os.path.join(cache_dir, cache_name)
    else:
        cache_path = None

    index_df = index_df.copy()

    # TODO: Smarter mechanism for cached file
    if cache_path is not None and os.path.isfile(cache_path) and not force:
        cdb_mjd_df_cache = pd.read_csv(
            cache_path, index_col=[0, 1], header=[0, 1]
        )
    else:
        cdb_mjd_df_cache = None

    if cdb_mjd_df_cache is None:
        # We load the headers from two possible extensions
        # NOTE: Might not be necessary in 0.7
        headers = index_df.FULLPATH.apply(fits.getheader, ext=0)
        ext1_mask = headers.str.len() == 4
        headers_ext1 = index_df.FULLPATH[ext1_mask].apply(
            fits.getheader, ext=1
        )
        headers[ext1_mask] = headers_ext1
        headers_df = pd.DataFrame(headers.tolist())

        # We only want the CDB* keys. Other useful keys are already in index
        keys = [
            params[kw][0] for kw in list(params) if kw.startswith("KW_CDB")
        ]

        # Only apply if CDB* keys are in the headers
        if all([item in headers_df.columns for item in keys]):
            cdb_df = headers_df[keys]
        else:
            return index_df
        unique_cdb_files = get_unique_vals(cdb_df)
        sep = "" if params["DRS_CALIB_DB"].endswith(os.path.sep) else os.path.sep
        unique_cdb_files_fullpath = params["DRS_CALIB_DB"] + sep + unique_cdb_files
        mjds = unique_cdb_files_fullpath.apply(get_hkey, args=("MJDMID",))
        mjds.index = unique_cdb_files
        mjd_df = cdb_df.replace(mjds.to_dict())

        cdb_mjd_df = pd.concat(
            [cdb_df, mjd_df], axis=1, keys=["CALIB_FILE", "DELTA_MJD"]
        ).swaplevel(0, 1, 1)

        cdb_mjd_df.index = pd.MultiIndex.from_frame(
            index_df.reset_index()[["NIGHTNAME", "FILENAME"]]
        )

        mjd_ind_files = index_df.FULLPATH.apply(get_hkey, args=("MJDMID",))
        cdb_mjd_df.loc[:, (slice(None), "DELTA_MJD")] = (
            cdb_mjd_df.loc[:, (slice(None), "DELTA_MJD")].values
            - mjd_ind_files.values[:, None]
        )  # using numpy because pandas multi-index operation did not work well

        if cache_path is not None:
            cdb_mjd_df.to_csv(cache_path)
    else:
        cdb_mjd_df = cdb_mjd_df_cache.copy()

    return cdb_mjd_df


def get_hkey(fname, hkey):
    if not isinstance(fname, str) or not fname.endswith(".fits"):
        return np.nan
    try:
        hdf = fits.getheader(fname, ext=0)
        mjd = hdf[hkey]
    except KeyError:
        hdf = fits.getheader(fname, ext=1)
        mjd = hdf[hkey]
    except KeyError:
        warnings.warn(
            f"Could not find {hkey} in extension 0 or 1, using nan",
            RuntimeWarning,
        )
        mjd = np.nan

    return mjd


def get_unique_vals(df, keepna=False):
    vals_list = []
    for col in df.columns:
        colvals = df[col]
        if not keepna:
            colmask = ~colvals.isnull() | (colvals == "None")
            colvals = colvals[colmask]
        un_vals = pd.unique(colvals)
        vals_list.append(un_vals)
    series = pd.Series(np.concatenate(vals_list))

    return series


def global_index_check(full_index: DataFrame, full_log: DataFrame):
    """
    Run check on index to make sure that can be associated to recipes
    NOTE: Might be useless or **very** simplified in v0.7+
    """

    full_index = full_index.copy()

    # ???: Should these keys not be hardcoded or ok once in index and log ?
    # Some masks for all characteristics we want to check
    nan_pid_mask = full_index.KW_PID.isnull()
    blank_pid_mask = full_index.KW_PID == "--"
    no_pid_mask = nan_pid_mask | blank_pid_mask
    in_log_mask = full_index.KW_PID.isin(full_log.PID)
    full_index["IN_LOG"] = in_log_mask
    index_problem_mask = no_pid_mask | (~in_log_mask)

    # Identify PID types
    full_index["PID_TYPE"] = "PID"
    full_index.loc[nan_pid_mask, "PID_TYPE"] = "NaN"
    full_index.loc[blank_pid_mask, "PID_TYPE"] = "Blank"

    # TODO: Do the checks/output here
    # TODO: This can be displayed in a bokeh table with some specific columns ?

    # TODO: When have way of knowing which recipe, return also based on recipe
    # (not just log).
    # THOMAS-FROM-A-FEW-MONTHS-LATER: not sure what I was talking about, 0.7 maybe ?
    return full_index[~index_problem_mask], full_index[index_problem_mask]


def load_log_df(
    output_parent: str,
    log_fname: str = "log.fits",
    return_missing: bool = False,
) -> DataFrame:
    """
    Load all log.fits files in single dataframe

    :param output_parent: Parent path of output directories
    :type output_parent: str
    :param log_fname: Log filename
    :type log_fname: str
    :param return_missing: Return list of missing log files
    :type return_missing: bool
    :return: Dataframe of log files
    :rtype: DataFrame
    """
    # TODO: In 0.7, this will change to loading a DB
    # ???: Keep info of output_parent in df ?

    allpaths = [
        os.path.join(output_parent, d, log_fname)
        for d in os.listdir(output_parent)
    ]
    log_df = load_fits_df(allpaths)

    # Use DIRECTORY as index and keep only relevant entries
    # whith self.recipe or extract recipes called
    log_df = log_df.set_index(["DIRECTORY"])

    if return_missing:
        # Dirs without log files
        missing_logs = [p for p in allpaths if not os.path.isfile(p)]
        return log_df, missing_logs
    else:
        return log_df


def load_index_df(
    output_parent: str,
    index_fname: str = "index.fits",
    return_missing: bool = False,
) -> DataFrame:
    """
    Load all index.fits files in a single dataframe

    :param output_parent: Parent path of output directories
    :type output_parent: str
    :param index_fname: index filename
    :type index_fname: str
    :param return_missing: Return list of missing index files
    :type return_missing: bool
    :return: Dataframe of index files
    :rtype: DataFrame
    """
    # TODO: In 0.7, this will change to loading a DB

    # ???: Keep info of output_parent in df ?
    # Get all index.fits in a dataframe
    allpaths = [
        os.path.join(output_parent, d, index_fname)
        for d in os.listdir(output_parent)
    ]
    ind_df = load_fits_df(allpaths)

    # Add full paths to dataframe
    parent_path = os.path.dirname(os.path.dirname(allpaths[0]))
    sep = os.path.sep
    ind_df["FULLPATH"] = (
        parent_path + sep + ind_df.NIGHTNAME + sep + ind_df.FILENAME
    )

    # Use NIGHTNAME as index
    ind_df = ind_df.set_index(["NIGHTNAME"])

    if return_missing:
        # Dirs without index.fits
        missing_inds = [p for p in allpaths if not os.path.isfile(p)]
        return ind_df, missing_inds
    else:
        return ind_df


def get_names_no_index(params: ParamDict) -> List[str]:
    """
    Get file (base)names that should not be in index dataframe.

    :param params: APERO params
    :type params: ParamDict
    :return: List of file names derived from APERO
    :rtype: List[str]
    """

    # Get paths that contain filenames we want to exclude
    calib_reset_path = param_functions.get_relative_folder(
        params["DRS_PACKAGE"], params["DRS_RESET_CALIBDB_PATH"]
    )
    tellu_reset_path = param_functions.get_relative_folder(
        params["DRS_PACKAGE"], params["DRS_RESET_TELLUDB_PATH"]
    )
    runs_reset_path = param_functions.get_relative_folder(
        params["DRS_PACKAGE"], params["DRS_RESET_RUN_PATH"]
    )

    # We know that log and index are also not data
    # NOTE: index/log will be in db in 0.7+
    exclude_list = ["index.fits", "log.fits"]
    exclude_list.extend(os.listdir(calib_reset_path))
    exclude_list.extend(os.listdir(tellu_reset_path))
    exclude_list.extend(os.listdir(runs_reset_path))

    return exclude_list


def get_output_files(
    output_parent: str, exclude_fname: Optional[Union[List[str], str]] = None
) -> Series:
    """
    Load all output files in a pandas series

    :param output_parent: Parent output directory to scan
    :type output_parent: str
    :param exclude_fname: File names to exclude (none by default).
    :type exclude_fname: Optional[Union[List[str], str]]
    :return: List of full path to output files
    :rtype: Series
    """
    if isinstance(exclude_fname, str):
        exclude_fname = [exclude_fname]
    all_paths = glob.glob(os.path.join(output_parent, "*", "*.fits"))
    out_series = pd.Series(all_paths)
    if exclude_fname is not None:
        exclude_mask = out_series.apply(os.path.basename).isin(exclude_fname)
        out_series = out_series[~exclude_mask]

    return out_series


def missing_index_headers(
    ind_df: DataFrame,
    output_files: Optional[pd.Series] = None,
    instrument: str = "SPIROU",
    force: bool = False,
    cache_dir: str = None,
    cache_suffix: str = None,
) -> Tuple[DataFrame, Series]:
    """
    Find files with missing index and store header info in a dataframe.
    The output is compatible with APERO index.

    :param ind_df: Dataframe with all index content
    :type ind_df: DataFrame
    :param output_files: List of output files,
                         automatically loaded if None (default)
    :type output_files: Optional[List[str]]
    :param instrument: Instrument name for APERO. Default is SPIROU.
    :type instrument: str
    :param force: Force reading all files even if cache is available,
                  defaults to False
    :type force: bool, optional
    :param cache_dir: Directory where the cache file should be found,
                      defaults to None
    :type cache_dir: str, optional
    :param cache_suffix: Suffix to add right before the cache file extension,
                         defaults to None
    :type cache_suffix: str, optional
    :return: Dataframe with info of missing index files,
             Series of files not found at all
    :rtype: Tuple[DataFrame, Series]
    """

    # TODO: Could have better naming/ID method s.t. repeated runs update OK
    if cache_dir is not None:
        p = Path(cache_dir)
        p.mkdir(parents=True, exist_ok=True)
        cache_name = f"missing_index_df_{instrument}{cache_suffix}.csv"
        cache_path = os.path.join(cache_dir, cache_name)
        nf_cache = f"not_found_{instrument}{cache_suffix}.csv"
        cache_not_found = os.path.join(cache_dir, nf_cache)
    else:
        cache_path = None

    # If not output files are given, we load them from disk,
    # using apero to discard some filenames
    params = constants.load(instrument)
    exclude_fname = get_names_no_index(
        params
    )  # Exclude log.fits and stuff like that

    if output_files is None:
        parent_dir = get_nth_parent(ind_df.FULLPATH.iloc[0], order=2)
        output_files = get_output_files(
            parent_dir, exclude_fname=exclude_fname
        )
    else:
        exclude_mask = output_files.apply(os.path.basename).isin(exclude_fname)
        output_files = output_files[~exclude_mask]

    # Debug files are not needed
    debug_mask = output_files.apply(os.path.basename).str.startswith("DEBUG")
    output_files = output_files[~debug_mask]

    # We reset the index to match created dataframe by default
    index_mask = output_files.apply(os.path.basename).isin(ind_df.FILENAME)
    out_not_in_index = output_files[~index_mask].reset_index(drop=True)

    # If there is a cached file, we load it and return it directly
    # TODO: Add a mechanism to use the cache properly, by combining and
    # re-saving, not just skipping reading any file is there is cache
    if (
        cache_path is not None
        and os.path.isfile(cache_path)
        and os.path.isfile(cache_not_found)
        and not force
    ):
        missing_ind_df_cache = pd.read_csv(cache_path, index_col=0)
        not_found_series_cache = pd.read_csv(
            cache_not_found, index_col=0, squeeze=True
        )
    else:
        missing_ind_df_cache = None

    if missing_ind_df_cache is None:
        # Some files (images) have file in ext=0, others (tables) in ext=1
        # NOTE: This may change in future versions ?
        # TODO: have this in two place: move to function ?

        # If no missing files
        if len(out_not_in_index) == 0:
            return out_not_in_index, pd.Series([])

        # We treat missing files on a per-basename basis
        out_basename = out_not_in_index.apply(os.path.basename)

        # We iterate unique values to find their obs night
        dup_mask = out_basename.duplicated()
        unique_base = out_basename[~dup_mask]
        unique_full = out_not_in_index[~dup_mask]

        # Get the header for each missing unique basename
        uheaders = unique_full.apply(fits.getheader, ext=0)
        ext1_mask = uheaders.str.len() == 4
        uheaders_ext1 = unique_full[ext1_mask].apply(fits.getheader, ext=1)
        uheaders[ext1_mask] = uheaders_ext1

        # Get the dates, then combine with basename to get path ends
        dkey = params["KW_DATE_OBS"][0]
        og_dates = uheaders.apply(lambda x: x[dkey])

        # Of the expected original dates, get the ones we can find
        og_ends = og_dates + os.path.sep + unique_base
        not_in_index_ends = (
            out_not_in_index.str.split(os.path.sep)
            .str[-2:]
            .str.join(os.path.sep)
        )
        found_mask = og_ends.isin(not_in_index_ends)
        prefix_path = get_nth_parent(unique_full.iloc[0], order=2)
        og_found = prefix_path + os.path.sep + og_ends[found_mask]
        og_missing = unique_base[~found_mask]

        headers = og_found.apply(fits.getheader, ext=0)
        ext1_mask = headers.str.len() == 4
        headers_ext1 = og_found[ext1_mask].apply(fits.getheader, ext=1)
        headers[ext1_mask] = headers_ext1

        # Get header keys corresponding to index columns
        pconstant = constants.pload(instrument)
        index_cols = pconstant.OUTPUT_FILE_HEADER_KEYS()
        keys = [params[col][0] for col in index_cols]

        # Get df in index format (tolist expands header values automatically)
        missing_headers_df = DataFrame(headers.tolist())
        if len(missing_headers_df) > 0:
            missing_ind_df = DataFrame(
                missing_headers_df[keys].values, columns=index_cols
            )
        else:
            missing_ind_df = DataFrame([], columns=index_cols)

        # Add fields that are not in the headers
        missing_ind_df["FILENAME"] = og_found.apply(os.path.basename)
        missing_ind_df["NIGHTNAME"] = og_found.apply(os.path.dirname).apply(
            os.path.basename
        )
        missing_ind_df["LAST_MODIFIED"] = og_found.apply(
            os.path.getmtime
        ).astype(
            str
        )  # APERO stores these times as string, so we convert them here
        missing_ind_df["FULLPATH"] = og_found

        missing_ind_df = missing_ind_df.set_index(["NIGHTNAME"])
        missing_ind_df = missing_ind_df[ind_df.columns]

        if cache_path is not None:
            missing_ind_df.to_csv(cache_path)
            og_missing.to_csv(cache_not_found)
    else:
        missing_ind_df = missing_ind_df_cache.copy()
        og_missing = not_found_series_cache.copy()

    return missing_ind_df, og_missing


def is_series(o: Any) -> bool:

    return isinstance(o, Series)
