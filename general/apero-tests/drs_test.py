"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import glob
import os
import warnings
from datetime import datetime
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from apero.core import constants
from apero.core.core.drs_argument import DrsArgument
from apero.core.core.drs_file import DrsFitsFile
from apero.core.core.drs_recipe import DrsRecipe
from astropy.io import fits
from jinja2 import Environment, FileSystemLoader, select_autoescape
from pandas import DataFrame

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
        self.fileargs = None  # Recipe args and kwargs that are files
        self.drs_files = None  # Recipe args and kwargs that are files
        self.name = "Unknown Test for unknown Recipe"  # Full test name
        self.test_id = "unknown_test"  # Short test ID
        self.instrument = instrument  # Instrument may be set as kwarg
        self.recipe_name = None
        self.params = None  # DRS parameters
        self.ismaster = False  # Is it a master recipe?
        self.output_hkeys = []  # Output header keys
        self.calibdb_keys = []  # Calibdb keys
        self.calib_hkeys = []  # All calib header keys
        self.input_path = ""  # Path to input dir
        self.output_path = ""  # Path to output dir
        self.log_df = None
        self.ind_df = None

        # Overwrite some (most!) parameters with recipe info automatically
        if drs_recipe is not None:
            # The recipe object and test ID
            self.recipe = DrsRecipe()
            self.recipe.copy(drs_recipe)
            self.fileargs = self.get_fileargs()
            self.drs_files = [f for farg in self.fileargs for f in farg.files]
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
            self.calibdb_keys = self.get_dbkeys("calibration")
            self.calib_hkeys = [
                self.recipe.drs_params[kw][0]
                for kw in list(self.recipe.drs_params)
                if kw.startswith("KW_CDB")
            ]  # Previous calibrations

            # TODO: Need to filter ind_df with outputs
            self.log_df = self.load_log_df(all_log_df=all_log_df)
            self.ind_df = self.load_ind_df(all_ind_df=all_index_df)

            # Load calibdb DF (master entries and calibdb files)
            # TODO: Will need to load previous calibs
            # NOTE: Seems will need to parse **all** headers (better outside I guess)
            self.master_calib_df = self._load_master_calibdb()

    # =========================================================================
    # Functions to load output information
    # =========================================================================
    def get_dbkeys(self, dbname: str) -> List[str]:
        """
        Get db keys for a given database (any valid APERO dbname)

        :param dbname: Database name (any valid APERO dbname)
        :type dbname: str
        :return: List of db keys for the recipe
        :rtype: List[str]
        """
        if self.drs_files is None:
            raise TypeError("Cannot load dbkeys if drs_files is 'None'.")

        keylist = []
        for rfile in self.drs_files:
            if rfile.dbname is None or rfile.dbkey is None:
                continue

            if rfile.dbname == dbname:
                if rfile.fibers is None:
                    # Some files have just one fiber, so dbkey works out of the box
                    keylist.append(rfile.get_dbkey())
                elif isinstance(rfile.fibers, list):
                    # if multiple fibers, need to set one by one (copy for safety)
                    for fiber in rfile.fibers:
                        rfile2 = rfile.newcopy()
                        rfile2.fiber = fiber
                        keylist.append(rfile2.get_dbkey())
                else:
                    raise TypeError(
                        f"Got an unexpected fibers attribute for file {rfile}"
                    )

        return keylist

    def get_fileargs(self) -> List[DrsArgument]:
        fileargs = [
            self.recipe.args[key]
            for key in self.recipe.args
            if self.recipe.args[key].dtype in ["file", "files"]
        ]
        fileargs.extend(
            [
                self.recipe.kwargs[key]
                for key in self.recipe.kwargs
                if self.recipe.kwargs[key].dtype in ["file", "files"]
            ]
        )

        return fileargs

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
        # NOTE: For 0.7 version
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
