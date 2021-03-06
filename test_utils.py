"""
General functions to use in apero tests
"""
import glob
import os

import pandas as pd
from astropy.table import Table

from typing import List, Union, Optional


def load_fits_df(pathlist: List[str]) -> pd.DataFrame:
    """
    Load contents from list of fits tables in a dataframe, making sure to
    decode strings

    :param pathlist:
    :type pathlist: List[str]

    :returns:
    :rtype: pd.DataFrame
    """
    # Get all logs in a dataframe
    paths = [p for p in pathlist if os.path.isfile(p)]  # Check all real files
    df_list = [Table.read(f).to_pandas() for f in paths]  # Check List of dfs
    df = pd.concat(df_list)  # List

    # Decode strings (fits from astropy are encoded)
    for col in df.columns:
        # try/except for non-string columns
        try:
            df[col] = df[col].str.decode("utf-8")
        except AttributeError:
            pass

    return df


def load_index_df(output_parent: str,
                  index_fname: str = "index.fits") -> pd.DataFrame:
    """
    Load all index.fits files in a single dataframe

    :param output_parent: Parent path of output directories
    :type output_parent: str
    :param index_fname: index filename
    :type index_fname: str
    :return: Dataframe
    :rtype: pd.DataFrame
    """
    # TODO: In 0.7, this will change to loading a DB

    # TODO: Keep info of output_parent in df
    # Get all index.fits in a dataframe
    allpaths = [
        os.path.join(output_parent, d, index_fname)
        for d in os.listdir(output_parent)
    ]
    ind_df = load_fits_df(allpaths)

    # Use NIGHTNAME as index
    ind_df = ind_df.set_index(["NIGHTNAME"])

    # Report paths without index.fits at all
    missing_inds = [p for p in allpaths if not os.path.isfile(p)]

    return ind_df, missing_inds


def list_output_files(
        output_parent: str,
        exclude_fname: Union[List[str], str] = ["index", "log"]) -> List[str]:
    """

    :param output_parent: Parent output directory to scan
    :type output_parent: str
    :param exclude_fname: File name prefixes to exclude
                          (i.e. file name without extension)
    :type exclude_fname: Union[List[str], str]
    :return: List of full path to output files
    :rtype: List[str]
    """
    if isinstance(exclude_fname, str):
        exclude_fname = [exclude_fname]
    exclude_pattern = '|'.join(exclude_fname)
    all_paths = glob.glob(
        os.path.join(output_parent, "*", f"[!{exclude_pattern}]*.fits"))

    return all_paths

def find_missing_index(ind_df: pd.DataFrame, output_files: Optional[List[str]] = None) -> pd.DataFrame:
    """

    :param ind_df: Dataframe with all index content
    :type ind_df: pd.DataFrame
    :param output_files: List of output files,
                         automatically loaded if None (default)
    :type output_files: Optional[List[str]]
    :return: Dataframe with info of missing index files
    :rtype: pd.DataFrame
    """
    if output_files is None:
        # TODO: Use ind_df output parent info to load this
        raise TypeError("Nonetype output_files not yet supported, please provide a path")