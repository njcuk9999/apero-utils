"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import os
import warnings
import glob
import numpy as np
from datetime import datetime
from typing import List, Optional, Tuple, Union

import pandas as pd
from astropy.io import fits
from jinja2 import Environment, FileSystemLoader, select_autoescape
from pandas import DataFrame

from apero.core import constants
from apero.core.core.drs_recipe import DrsRecipe

import test_utils as ut

# from . import OUTDIR, TEMPLATEDIR

# TODO: Maybe we can move this so it's accessible to all files, maybe parameter
# in drs later
PARENTDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEMPLATEDIR = os.path.join(PARENTDIR, "templates")
OUTDIR = os.path.join(PARENTDIR, "out")


def removext(name: str, ext: str = ".py") -> str:
    """
    Remove extension of a file/recipe name

    :param name: file/recipe name
    :type name: str
    :param ext: filetype extension, default is ".py"
    :type ext: str
    :returns: cleaned up name
    :rtype:  str
    """
    if not ext.startswith("."):
        ext = "." + ext

    while name.endswith(ext):
        name = name[:-3]

    return name


class DrsTest:
    def __init__(
        self,
        instrument: Optional[str] = None,
        drs_recipe: Optional[DrsRecipe] = None,
        setup: Optional[str] = None,
        testnum: int = 1,
        all_log_df: Optional[DataFrame] = None,
        all_index_df: Optional[DataFrame] = None,
    ):
        """

        :param instrument: Instrument of the tested recipe
        :type instrument: Optional[str]
        :param drs_recipe: Recipe under test
        :type drs_recipe: Optional[DrsRecipe]
        :param setup: Path APERO setup
        :type setup: Optional[str]
        :param testnum: Number ID N of the test (Nth test for recipe),
                        default is 1
        :type testnum: int
        """
        # Get setup path
        if setup is None:
            self.setup = os.environ["DRS_UCONFIG"]
        else:
            self.setup = setup

        # Set date at start of test
        self.date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # test date

        # Set default parameter values, can overwrite later
        self.recipe = None  # Recipe object for test
        self.name = "Unknown Test for unknown Recipe"  # Full test name
        self.test_id = "unknown_test"  # Short test ID
        self.instrument = instrument  # Instrument may be set as kwarg
        self.recipe_name = None
        self.params = None  # DRS parameters
        self.ismaster = False  # Is it a master recipe?
        self.fibers = []  # TODO: Implement further only if needed in methods
        self.output_hkeys = []  # Output header keys
        self.calibdb_keys = []  # Calibdb keys
        self.calib_hkeys = []  # All calib header keys
        self.calls_extract = False  # TODO: Remove if not needed
        self.input_path = ""  # Path to input dir
        self.output_path = ""  # Path to output dir
        self.log_df = None
        self.ind_df = None

        # Overwrite some (most!) parameters with recipe info automatically
        if drs_recipe is not None:
            # The recipe object and test ID
            self.recipe = drs_recipe
            if self.instrument is not None:
                warnings.warn(
                    "Overwriting kwarg instrument with recipe info", RuntimeWarning
                )
            self.instrument = self.recipe.instrument
            self.recipe_name = removext(self.recipe.name, ext=".py")
            self.test_id = self.recipe_name + f"_test{testnum}"
            self.name = f"Test for recipe {self.recipe_name}"

            # General params
            self.params = self.recipe.drs_params
            self.ismaster = self.recipe.master

            # Path to input and output directories
            self.dirpaths = self.get_dir_path()
            self.input_path = self.dirpaths[self.recipe.inputdir]
            self.output_path = self.dirpaths[self.recipe.outputdir]

            # Get keys for outputs
            self.output_dict = self.recipe.outputs
            self.output_hkeys = list(
                map(
                    lambda o: o.required_header_keys["KW_OUTPUT"],
                    list(self.output_dict.values()),
                )
            )
            self.calibdb_keys = list(
                map(lambda o: o.get_dbkey(), list(self.output_dict.values()))
            )
            self.calib_hkeys = [
                self.recipe.drs_params[kw][0]
                for kw in list(self.recipe.drs_params)
                if kw.startswith("KW_CDB")
            ]

            self.log_df = self.load_log_df(all_log_df=all_log_df)
            self.ind_df = self.load_ind_df(all_ind_df=all_index_df)

            # Load calibdb DF (master entries and calibdb files)
            # TODO: Do we need to restructure like log and ind here ?
            self.master_calib_df = self._load_master_calibdb()

    # =========================================================================
    # Functions to load output information
    # =========================================================================
    def load_log_df(
        self, all_log_df: Optional[DataFrame] = None, force: bool = False
    ) -> Union[DataFrame, None]:
        """Get log content in a dataframe

        Parse all log files and return a dataframe with only entries that have
        the current recipe in LOGFILE.

        :param all_log_df: Dataframe with log info for multiple recipes
        :type all_log_df: DataFrame
        :param force: Force to reload, by default (False) returns self.log_df if exists
        :type force: bool
        :returns: dataframe of log content and list of missing log files
        :rtype: Tuple[pd.DataFrame, list]
        """
        # NOTE: missing_log was removed, list of all missing logs can be generated with
        #  util load_log_df function by setting `return_missing` to True

        if not force and self.log_df is not None:
            return self.log_df

        if self.recipe_name is None:
            warnings.warn(
                "Cannot load log dataframe if no recipe is set, returning None",
                RuntimeWarning,
            )
            return None

        if all_log_df is None:
            all_log_df = ut.load_log_df(self.output_path)

        # Important for extract recipes: keep only the recipe under test
        log_df = all_log_df[all_log_df.LOGFILE.str.contains(self.recipe_name)]

        return log_df

    def load_ind_df(
        self, all_ind_df: Optional[DataFrame] = None, force: bool = True
    ) -> Union[DataFrame, None]:
        """Get index contents in a dataframe

        Parse all index.fits files and return a dataframe.

        :param all_ind_df: Dataframe with all index entries for multiple recipes
        :type all_ind_df: Optional[DataFrame]
        :param force: Force to reload, by default (False) returns self.log_df if exists
        :type force: bool
        :returns: dataframe of index content and list of missing index files
        :rtype: Tuple[pd.DataFrame, list]
        """
        # NOTE: missing_ind was removed, list of all missing inds can be generated with
        #  util load_ind_df function by setting `return_missing` to True

        if not force and self.ind_df is not None:
            return self.ind_df

        if self.recipe_name is None:
            # TODO: Make sure this is true when done
            warnings.warn(
                "Cannot load index dataframe if no recipe is set, returning None",
                RuntimeWarning,
            )
            return None

        if all_ind_df is None:
            all_ind_df = ut.load_index_df(self.output_path)

        # TODO: Likely need to filter recipes here
        # ???: Is extract in keywords or need extra ?
        ind_df = all_ind_df.copy()

        return ind_df

    def _load_master_calibdb(self):
        """
        Load Master calibdb

        :return: Dataframe with master calibdb entries
        :rtype: pd.DataFrame
        """
        # PORT: For 0.7 version
        # from apero.core.core import drs_database
        # from apero.core import constants
        # params = constants.load()
        # calibdb = drs_database.CalibDatabase(params)
        # calibdb.load_db()
        # calib_table = calibdb.database.get('*', calibdb.database.tname,
        #                                    return_pandas=True)

        calibdb_file_path = os.path.join(
            self.params["DRS_CALIB_DB"], self.params["CALIB_DB_NAME"]
        )

        colnames = self.params.listp("CALIB_DB_COLS", dtype=str)
        values = np.loadtxt(calibdb_file_path, unpack=True, dtype=str)
        calib_df = pd.DataFrame(dict(zip(colnames, values)))

        return calib_df

    def find_missing_index(self):
        """
        Get df of files not in index.fits but found on disk
        """

        # TODO: Make sure that loading files here is optimized/moved if bottleneck
        # Search for files not in the index.fits
        # Paths of the files inside the index.fits
        indpaths = (
            self.output_path
            + self.ind_df.index.get_level_values("NIGHTNAME")
            + os.path.sep
            + self.ind_df.FILENAME
        )
        # Paths for all .fits files on disk
        filepaths = [
            d
            for d in glob.glob(self.output_path + "*/*.fits")
            if not (
                os.path.basename(d).startswith("index")
                or os.path.basename(d).startswith("log")
            )
        ]
        # Paths for files not in index.fits
        missing_paths = np.setdiff1d(filepaths, indpaths.tolist())
        # Load headers
        headers0 = pd.Series(missing_paths).apply(fits.getheader, ext=0)
        index_ext1 = np.where(pd.Index(headers0.str.len()).get_loc(4))[0]
        headers1 = pd.Series(missing_paths[index_ext1], index=index_ext1).apply(
            fits.getheader, ext=1
        )
        headers = headers0
        headers.update(headers1)
        header_df = pd.DataFrame(headers.tolist(), index=headers.index)

        pconstant = constants.pload(self.instrument)
        # list of index.fits columns
        index_cols = pconstant.OUTPUT_FILE_HEADER_KEYS()
        # Header keywords associate with index.fits columns
        keys = [self.params[col][0] for col in index_cols]

        missing_ind_df = pd.DataFrame(header_df[keys].values, columns=index_cols)

        # Populate the three columns not in the header
        filename = []
        nightname = []
        mtime = []
        for i in range(len(missing_ind_df)):
            sep_index = (
                missing_paths[i].replace(self.output_path, "").index(os.path.sep)
            )
            filename.append(
                missing_paths[i].replace(self.output_path, "")[(sep_index + 1) :]
            )
            nightname.append(missing_paths[i].replace(self.output_path, "")[:sep_index])
            mtime.append(os.path.getmtime(missing_paths[i]))
        missing_ind_df.insert(loc=0, column="LAST_MODIFIED", value=mtime)
        missing_ind_df.insert(loc=0, column="NIGHTNAME", value=nightname)
        missing_ind_df.insert(loc=0, column="FILENAME", value=filename)
        missing_ind_df = missing_ind_df.set_index(["NIGHTNAME"])

        return missing_ind_df

    # =========================================================================
    # Properties of the test
    # =========================================================================

    # =========================================================================
    # Utility functions
    # =========================================================================
    def get_dir_path(self):
        """
        Get dictionary mapping keywords to data directories
        """
        dirpaths = dict()
        if self.params is not None:
            dirpaths["raw"] = self.params["DRS_DATA_RAW"]
            dirpaths["tmp"] = self.params["DRS_DATA_WORKING"]
            dirpaths["red"] = self.params["DRS_DATA_REDUC"]
            dirpaths["reduced"] = self.params["DRS_DATA_REDUC"]
        else:
            dirpaths["raw"] = ""
            dirpaths["tmp"] = ""
            dirpaths["red"] = ""
            dirpaths["reduced"] = ""

        return dirpaths

    # =========================================================================
    # Function to run tests and write their output
    # =========================================================================
    def run_test(self):
        """Run the test given a proper recipe test script"""
        # TODO: General runtest method for external calls
        raise NotImplementedError("The run_test method has not yet been implemented")

    def gen_html(self, html_dict: dict):
        """Generate HTML summary from jinja2 template.

        :param html_dict: dict with all values used in html template
        :type html_dict: dict
        """

        # Jinja2 env
        env = Environment(
            loader=FileSystemLoader(TEMPLATEDIR),
            autoescape=select_autoescape(["html", "xml"]),
        )

        # Create template for test
        template = env.get_template(".".join([self.test_id, "html"]))

        html_text = template.render(html_dict)

        output_path = os.path.join(
            OUTDIR, self.test_id, ".".join([self.test_id, "html"])
        )
        with open(output_path, "w") as f:
            f.write(html_text)
