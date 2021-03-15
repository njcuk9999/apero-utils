"""
General functions to use in apero tests
"""
import glob
import os

import pandas as pd
from astropy.io import fits
from astropy.table import Table

from typing import List, Union, Optional

from pandas import DataFrame

from apero.core import constants


def get_nth_parent(path: str, order: int = 1):
    """
    Get nth order parent of path by applying os.dirname repeatedly.

    :param path: Input path
    :type path: str
    :param order: Number of directories up (number of times os.dirname is applied).
                  Default is 1.
    :type order: int
    """
    parent_path = path
    for _ in range(order):
        parent_path = os.path.dirname(parent_path)

    return parent_path


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


def load_log_df(
    output_parent: str, log_fname: str = "log.fits", return_missing: bool = False
) -> DataFrame:
    """
    Load all index.fits files in single dataframe

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
        os.path.join(output_parent, d, log_fname) for d in os.listdir(output_parent)
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
    output_parent: str, index_fname: str = "index.fits", return_missing: bool = False
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
        os.path.join(output_parent, d, index_fname) for d in os.listdir(output_parent)
    ]
    ind_df = load_fits_df(allpaths)

    # Add full paths to dataframe
    parent_path = os.path.dirname(os.path.dirname(allpaths[0]))
    sep = os.path.sep
    ind_df["FULLPATH"] = parent_path + sep + ind_df.NIGHTNAME + sep + ind_df.FILENAME

    # Use NIGHTNAME as index
    ind_df = ind_df.set_index(["NIGHTNAME"])

    if return_missing:
        # Dirs without index.fits
        missing_inds = [p for p in allpaths if not os.path.isfile(p)]
        return ind_df, missing_inds
    else:
        return ind_df


def get_output_files(
    output_parent: str, exclude_fname: Optional[Union[List[str], str]] = None
) -> pd.Series:
    """

    :param output_parent: Parent output directory to scan
    :type output_parent: str
    :param exclude_fname: File name prefixes to exclude
                          (i.e. file name without extension).
                          By default, index and log files are ignored.
    :type exclude_fname: Optional[Union[List[str], str]]
    :return: List of full path to output files
    :rtype: List[str]
    """
    if exclude_fname is None:
        exclude_fname = ["index", "log"]
    if isinstance(exclude_fname, str):
        exclude_fname = [exclude_fname]
    exclude_pattern = "|".join(exclude_fname)
    all_paths = glob.glob(
        os.path.join(output_parent, "*", f"[!{exclude_pattern}]*.fits")
    )
    out_series = pd.Series(all_paths)

    return out_series


def missing_index_headers(
    ind_df: DataFrame,
    output_files: Optional[pd.Series] = None,
    instrument: str = "SPIROU",
) -> DataFrame:
    """
    Find files with missing index and store header info in headers

    :param ind_df: Dataframe with all index content
    :type ind_df: DataFrame
    :param output_files: List of output files,
                         automatically loaded if None (default)
    :type output_files: Optional[List[str]]
    :param instrument: Instrument name for APERO. Default is SPIROU.
    :type instrument: str
    :return: Dataframe with info of missing index files
    :rtype: DataFrame
    """
    if output_files is None:
        parent_dir = get_nth_parent(ind_df.FULLPATH.iloc[0], order=2)
        output_files = get_output_files(parent_dir)

    # We reset the index to match created dataframe by default
    index_mask = output_files.isin(ind_df.FULLPATH)
    out_not_in_index = output_files[~index_mask].reset_index(drop=True)

    # Some files (images) have file in ext=0, others (tables) in ext=1
    # NOTE: This may change in future versions ?
    headers = out_not_in_index.apply(fits.getheader, ext=0)
    ext1_mask = headers == 4
    headers_ext1 = out_not_in_index[ext1_mask].apply(fits.getheader, ext=1)
    headers[ext1_mask] = headers_ext1

    # Get header keys corresponding to index columns
    pconstant = constants.pload(instrument)
    params = constants.load(instrument)
    index_cols = pconstant.OUTPUT_FILE_HEADER_KEYS()
    keys = [params[col][0] for col in index_cols]

    # Get dataframe in index format (tolist expands header values automatically)
    missing_headers_df = DataFrame(headers.tolist())
    missing_ind_df = DataFrame(missing_headers_df[keys].values, columns=index_cols)

    # Add fields that are not in the headers
    missing_ind_df["FILENAME"] = out_not_in_index.apply(os.path.basename)
    missing_ind_df["NIGHTNAME"] = out_not_in_index.apply(os.path.dirname).apply(
        os.path.basename
    )
    missing_ind_df["LAST_MODIFIED"] = out_not_in_index.apply(os.path.getmtime).astype(
        str
    )  # APERO stores these time as string
    missing_ind_df["FULLPATH"] = out_not_in_index

    missing_ind_df = missing_ind_df.set_index(["NIGHTNAME"])
    missing_ind_df = missing_ind_df[ind_df.columns]

    return missing_ind_df
