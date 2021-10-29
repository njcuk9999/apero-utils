"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import os
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

import pandas as pd
from apero.core import constants
from apero.core.core.drs_argument import DrsArgument
from apero.core.core.drs_file import DrsFitsFile
from apero.core.core.drs_recipe import DrsRecipe
from jinja2 import Environment, FileSystemLoader, select_autoescape
from pandas import DataFrame

import apero_tests.subtest as st
import apero_tests.utils as ut

# from . import OUTDIR, TEMPLATEDIR

# TODO: Maybe we can move this so it's accessible to all files, maybe parameter
# in drs later
PARENTDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEMPLATEDIR = os.path.join(PARENTDIR, "templates")
OUTDIR = os.path.join(PARENTDIR, "out")
CACHEDIR = os.path.join(PARENTDIR, "cache")

# TODO: Maybe link to full docs, not github
RECIPE_DOCS_LINK = "https://github.com/njcuk9999/apero-drs#8-APERO-Recipes"


class DrsTest:
    def __init__(
        self,
        instrument: Optional[str] = None,
        drs_recipe: Optional[DrsRecipe] = None,
        setup: Optional[str] = None,
        testnum: int = 1,
        pp_flag: bool = False,
        all_log_df: Optional[DataFrame] = None,
        all_index_df: Optional[DataFrame] = None,
        all_master_calib_df: Optional[DataFrame] = None,
        all_tellu_df: Optional[DataFrame] = None,
        all_cdb_used_df: Optional[DataFrame] = None,
    ):
        """
        Test class that contains information about test results for a given
        recipe.

        The all_* dataframes are useful if many recipes are checked and users
        want to parse the database only once.

        :param instrument: Instrument of the tested recipe
        :type instrument: Optional[str]
        :param drs_recipe: Recipe under test
        :type drs_recipe: Optional[DrsRecipe]
        :param setup: Path APERO setup
        :type setup: Optional[str]
        :param testnum: Number ID N of the test (Nth test for recipe),
                        default is 1
        :type testnum: int
        :pp_flag: TEMPORARY flag for PP tests when different
        :param all_log_df: DataFrame with logs for all recipes
        :type  all_log_df: Optional[DataFrame]
        :param all_index_df: DataFrame with index for all recipes
        :type all_index_df: Optional[DataFrame]
        :param all_master_calib_df: DataFrame with master calib info for all
                                    recipes
        :type all_master_calib_df: Optional[DataFrame]
        :param all_tellu_df: DataFrame with tellu info for all recipes
        :type all_tellu_df: Optional[DataFrame]
        :param all_cdb_used_df: DataFrame with info about used calibs for all
                                recipes
        :type all_cdb_used_df: Optional[DataFrame]
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
        self.arg_drs_files = None  # Recipe args and kwargs that are files
        self.name = "Unknown Test for unknown Recipe"  # Full test name
        self.test_id = "unknown_test"  # Short test ID
        self.instrument = instrument  # Instrument may be set as kwarg
        self.recipe_name = None
        self.params = None  # DRS parameters
        self.pconst = None  # DRS constants
        self.ismaster = False  # Is it a master recipe?
        self.output_dict = None
        self.output_drs_files = None
        self.output_hkeys = []  # Output header keys
        self.calibdb_keys = []  # Calibdb keys
        self.input_path = ""  # Path to input dir
        self.output_path = ""  # Path to output dir
        self.html_path = ""  # Path to HTML report
        self.log_df = None
        self.ind_df = None
        self.calib_df = None
        self.tellu_df = None
        self.cdb_used_df = None
        self.subtest_list = None

        # Overwrite some (most!) parameters with recipe info automatically
        if drs_recipe is not None:
            # The recipe object and test ID
            self.recipe = DrsRecipe()
            self.recipe.copy(drs_recipe)
            self.fileargs = self.get_fileargs()
            self.arg_drs_files = [
                f for farg in self.fileargs for f in farg.files
            ]
            if self.instrument is not None:
                warnings.warn(
                    "Overwriting kwarg instrument with recipe info",
                    RuntimeWarning,
                )
            self.instrument = self.recipe.instrument
            self.recipe_name = ut.removext(self.recipe.name, ext=".py")
            self.test_id = self.recipe_name + f"_test{testnum}"
            self.name = f"Test for recipe {self.recipe_name}"

            # General params
            self.params = self.recipe.drs_params
            self.pconst = constants.pload(self.instrument)
            self.ismaster = self.recipe.master

            # PP flag
            self.pp_flag = pp_flag

            # Path to input and output directories
            self.dirpaths = self.get_dir_path()
            self.input_path = self.dirpaths[self.recipe.inputdir]
            self.output_path = self.dirpaths[self.recipe.outputdir]
            self.calibdb_path = os.path.join(
                self.params["DRS_CALIB_DB"], self.params["CALIB_DB_NAME"]
            )

            # Path to HTML report
            self.html_path = Path(
                OUTDIR, self.test_id, ".".join([self.test_id, "html"])
            )

            # Get keys for outputs
            self.output_dict = self.recipe.outputs
            self.output_drs_files = [
                o
                for o in list(self.output_dict.values())
                if isinstance(o, DrsFitsFile)
            ]

            # TODO: Handle all recipes without KW_OUTPUT in empty
            # required_header_keys
            if not self.pp_flag:
                self.output_hkeys = list(
                    map(
                        lambda o: o.required_header_keys["KW_OUTPUT"],
                        self.output_drs_files,
                    )
                )

            # NOTE: index is filtered with log to keep only relevant entries.
            # A global test is used to report outputs that match no logs
            self.log_df = self.load_log_df(all_log_df=all_log_df)
            self.ind_df = self.load_ind_df(all_ind_df=all_index_df)

            # Add QC informatio to index (useful for calibdb comparisons)
            pid_qc_mapping = self.log_df.groupby("PID").PASSED_ALL_QC.all()
            qc_mask = pid_qc_mapping.loc[self.ind_df.KW_PID]
            self.ind_df["QC_PASSED"] = qc_mask.values

            self.calibdb_keys = self.get_dbkeys("calibration")
            self.telludb_keys = self.get_dbkeys("telluric")
            if not self.pp_flag:
                self.calibdb_to_output = self.get_db_output_mapping(
                    "calibration"
                )

                # This addds CALIB_KEY to index df for future comparisons
                self.add_calib_to_ind()

            # Load calibdb DF (master entries and calibdb files)
            self.calib_df = self.load_calib_df(
                all_calib_df=all_master_calib_df
            )

            if not self.pp_flag:
                self.cdb_used_df = self.load_cdb_used_df(
                    all_cdb_used_df=all_cdb_used_df
                )

            # Load telludb DF
            self.tellu_df = self.load_tellu_df(all_tellu_df=all_tellu_df)

            # Generate list of default subtests
            self.set_subtests()

    # =========================================================================
    # Functions to load output information
    # =========================================================================
    def get_dbkeys(self, dbname: str) -> List[str]:
        """
        Get output db keys for a given database (any valid APERO dbname)

        :param dbname: Database name (any valid APERO dbname)
        :type dbname: str
        :return: List of db keys for the recipe
        :rtype: List[str]
        """
        if self.output_drs_files is None:
            raise TypeError("Cannot load dbkeys if drs_files is 'None'.")

        keylist = []
        for rfile in self.output_drs_files:
            if rfile.dbname is None or rfile.dbkey is None:
                continue

            if rfile.dbname == dbname:
                if rfile.fibers is None:
                    # Some files have just one fiber, dbkey works directly
                    keylist.append(rfile.get_dbkey())
                elif isinstance(rfile.fibers, list):
                    # If multiple fibers, need to set one by one
                    for fiber in rfile.fibers:
                        rfile2 = rfile.newcopy()  # copy for safety
                        rfile2.fiber = fiber
                        keylist.append(rfile2.get_dbkey())
                else:
                    raise TypeError(
                        f"Got an unexpected fibers attribute for file {rfile}"
                    )

        return keylist

    def get_db_output_mapping(self, dbname: str) -> Dict[str, str]:
        """
        Get mapping from output db keys to their corresponding KW_OUTPUT

        :param dbname: Database name (any valid APERO dbname)
        :type dbname: str
        :return: Mapping from dbkey to KW_OUTPUT
        :rtype: Dict[str, str]
        """
        if self.output_drs_files is None:
            raise TypeError("Cannot load dbkeys if drs_files is 'None'.")

        keydict = dict()
        for rfile in self.output_drs_files:
            out_key = rfile.required_header_keys["KW_OUTPUT"]
            if rfile.dbname is None or rfile.dbkey is None:
                continue

            if rfile.dbname == dbname:
                if rfile.fibers is None:
                    # Some files have just one fiber, dbkey works directly
                    db_key = rfile.get_dbkey()
                    keydict[db_key] = out_key
                elif isinstance(rfile.fibers, list):
                    # If multiple fibers, need to set one by one
                    for fiber in rfile.fibers:
                        rfile2 = rfile.newcopy()  # copy for safety
                        rfile2.fiber = fiber
                        db_key = rfile2.get_dbkey()
                        keydict[db_key] = out_key
                else:
                    raise TypeError(
                        f"Got an unexpected fibers attribute for file {rfile}"
                    )

        return keydict

    def get_fileargs(self) -> List[DrsArgument]:
        """
        Get list of argument with type 'file[s]' for the recipe
        """
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
        :param force: Force to reload, by default (False) returns self.log_df
                      if exists
        :type force: bool :returns: dataframe of log content and list of
                           missing log files
        :rtype: Tuple[pd.DataFrame, list]
        """
        # NOTE: missing_log was removed.
        #  Use utils.load_log_df with `return_missing=True` for that.

        # Don't reload if told not to
        if not force and self.log_df is not None:
            return self.log_df

        if self.recipe_name is None:
            raise ValueError("Cannot load log dataframe if no recipe is set")

        if all_log_df is None:
            all_log_df = ut.load_log_df(self.output_path)

        # Important for extract recipes: keep only the recipe under test
        log_df = all_log_df[all_log_df.LOGFILE.str.contains(self.recipe_name)]

        if self.pp_flag:
            log_df["ARGS"] = log_df.ARGS.str.replace("persi_", "")
            log_df["ODOMETER"] = log_df.ARGS.str.split(".fits").str[0].str[-8:]

        return log_df

    def load_ind_df(
        self, all_ind_df: Optional[DataFrame] = None, force: bool = True
    ) -> Union[DataFrame, None]:
        """Get index contents in a dataframe

        Parse all index.fits files and return a dataframe.

        :param all_ind_df: Dataframe with all index entries for multiple
                           recipes
        :type all_ind_df: Optional[DataFrame]
        :param force: Force to reload, by default (False) returns self.log_df
                      if exists
        :type force: bool
        :returns: dataframe of index content and list of missing index files
        :rtype: Tuple[pd.DataFrame, list]
        """
        # NOTE: missing_ind was removed.
        #  Use utils.load_ind_df with `return_missing=True` for that.

        # Don't reload if told not to
        if not force and self.ind_df is not None:
            return self.ind_df

        # TODO: Make sure this is true when done
        if self.recipe_name is None:
            raise ValueError("Cannot load index dataframe if no recipe is set")

        if all_ind_df is None:
            all_ind_df = ut.load_index_df(self.output_path)

        if self.log_df is None:
            raise AttributeError(
                "A non-null log_df is required to filter the index"
            )

        # We use the log to filter the index, also keep only output hkeys
        # NOTE: Might have an indepent way that does not use log in v0.7
        ind_df = all_ind_df[all_ind_df.KW_PID.isin(self.log_df.PID)]

        if not self.pp_flag:
            ind_df = ind_df[ind_df.KW_OUTPUT.isin(self.output_hkeys)]

        return ind_df

    def load_calib_df(
        self, all_calib_df: Optional[DataFrame] = None, force: bool = True
    ) -> Union[DataFrame, None]:
        """
        Load Master calibdb

        :return: Dataframe with master calibdb entries
        :rtype: pd.DataFrame
        """

        if not force and self.calib_df is not None:
            return self.ind_df

        if all_calib_df is None:
            all_calib_df = ut.load_db("CALIB")

        if self.calibdb_keys is None:
            msg = (
                "Cannot load index dataframe if no calibdb keys are set,"
                " returning None"
            )
            warnings.warn(msg, RuntimeWarning)
            return None

        # NOTE: Copy is used here so that future assignemnt never affect
        #       original (suppressing SettingWithCopyWarning)
        calib_df = all_calib_df[
            all_calib_df["key"].isin(self.calibdb_keys)
        ].copy()

        # ???: Maybe this won't be required with pandas >= 1.0 and
        # convert_dtypes (see utils function)
        calib_df["master"] = calib_df["master"].astype(int)

        return calib_df

    def load_tellu_df(
        self, all_tellu_df: Optional[DataFrame] = None, force: bool = True
    ) -> Union[DataFrame, None]:
        """
        Load telludb

        :return: Dataframe with telludb entries
        :rtype: pd.DataFrame
        """
        if not force and self.tellu_df is not None:
            return self.ind_df

        if all_tellu_df is None:
            all_tellu_df = ut.load_db("TELLU")

        if self.telludb_keys is None:
            msg = (
                "Cannot load index dataframe if no telludb keys are set,"
                " returning None"
            )
            warnings.warn(msg, RuntimeWarning)
            return None

        tellu_df = all_tellu_df[all_tellu_df["key"].isin(self.telludb_keys)]

        return tellu_df

    # def load_db_output(self, db_path: str):
    #     full_db_list = os.listdir(db_path)

    def load_cdb_used_df(
        self, all_cdb_used_df: Optional[DataFrame] = None, force: bool = True
    ) -> Union[DataFrame, None]:
        if not force and self.cdb_used_df is not None:
            return self.cdb_used_df

        if self.ind_df is None:
            msg = (
                "Cannot load previous calibs if index dataframe is not set."
                " Returning None"
            )
            warnings.warn(msg, RuntimeWarning)

        if all_cdb_used_df is None:
            all_cdb_used_df = ut.get_cdb_df(self.ind_df, self.params)

        # We keep only calib files whose index is in the APERO index dataframe
        if len(self.ind_df) > 0:
            ind_night_file = self.ind_df.reset_index()[
                ["NIGHTNAME", "FILENAME"]
            ]
            keep_ind = pd.MultiIndex.from_frame(ind_night_file)
            cdb_used_df = all_cdb_used_df.loc[keep_ind]
        else:
            cdb_used_df = all_cdb_used_df.iloc[0:0]

        return cdb_used_df

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

    def add_calib_to_ind(self):
        if "loc" in self.name:
            import ipdb; ipdb.set_trace()

        # Get all possible fibers
        sci_fibers, ref_fiber = self.pconst.FIBER_KINDS()
        all_fibers = sci_fibers + [ref_fiber]

        fiber_series = "_" + self.ind_df.KW_FIBER
        fiber_series = fiber_series.where(
            fiber_series.isin([f"_{fib}" for fib in all_fibers]), ""
        )
        out_fiber_series = self.ind_df.KW_OUTPUT + fiber_series

        # Map output keywords to calibdb keys (including fiber as suffix)
        out_to_calib = dict()
        for ck, ok in self.calibdb_to_output.items():

        out_to_calib = {
            (
                v
                if v in out_fiber_series.values
                else v + f"_{k.split('_')[-1]}"
            ): k
            for k, v in self.calibdb_to_output.items()
        }

        self.ind_df["CALIB_KEY"] = out_fiber_series.map(out_to_calib)

    # =========================================================================
    # Function to run tests and write their output
    # =========================================================================
    def set_subtests(self):
        """
        Create a all subtestes for this recipe and put them in a list
        """
        subtest_list = []

        # Count all raw files
        if self.pp_flag:
            subtest_list.append(st.CountRawTest(self.input_path))

        # Count log entries
        subtest_list.append(st.CountLogTest(self.log_df))

        # Count all outputs
        subtest_list.append(
            st.CountOutTest(self.ind_df, self.output_path, self.output_hkeys)
        )

        # Count unique outputs
        subtest_list.append(
            st.CountOutTest(
                self.ind_df, self.output_path, self.output_hkeys, unique=True
            )
        )

        # TODO Compare output and log results

        # QC failed count
        subtest_list.append(st.CountQCTest(self.log_df, self.html_path))

        # QC Plot
        subtest_list.append(
            st.PlotQCTest(self.log_df, self.html_path, self.recipe_name)
        )

        # Not ended count
        st_ended = st.CountEndedTest(self.log_df, self.html_path)
        subtest_list.append(st_ended)

        if not self.pp_flag:
            st_index_calib = st.CountIndexCalib(self.ind_df)
            subtest_list.append(st_index_calib)

            # calibDB entries count
            st_cdb_entries = st.CountCalibEntries(
                self.calib_df, self.calibdb_keys
            )
            subtest_list.append(st_cdb_entries)

            # TODO: Compare calibdb checks
            subtest_list.append(
                st.ComparisonTest(st_index_calib, st_cdb_entries)
            )

        # TODO: Calibdb

        # TODO: Telludb

        # telluDB output count
        subtest_list.append(st.CountTelluEntries(self.tellu_df))

        self.subtest_list = subtest_list

    def run_test(self):
        """
        Run all subtests and generate HTML summary for the recipe
        """

        # Run all subtests one by one
        final_list = []
        for i, subtest in enumerate(self.subtest_list):
            # Ensure unique subtest IDs within one DrsTest
            # This is important for subtests that have their own directory
            # or file. Too keep informative IDs when possible, we add number
            # to existing ID
            subtest.id = subtest.id + str(i)
            # Some might have outputs directly assiged and no run -> try/except
            try:
                subtest.run()
            except NotImplementedError:
                pass

            # No need to keep "None" results. Means the subtest has no info
            if subtest.result is not None:
                final_list.append(subtest)
        subtest_list = final_list

        html_dict = {
            # Summary header info
            "name": self.name,
            "setup": self.setup,
            "instrument": self.instrument,
            "recipe": self.recipe_name,
            "date": self.date,
            "output_path": self.output_path,
            "output_list": self.output_hkeys,
            "calibdb_list": self.calibdb_keys,
            "calibdb_path": self.calibdb_path,
            # TODO: Get links per recipe automatically
            "docs_link": RECIPE_DOCS_LINK,
            # Checks
            "subtest_list": subtest_list,
        }

        self.gen_html(html_dict)

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
        env.tests["series"] = ut.is_series

        # Create template for test
        # template = env.get_template(".".join([self.test_id, "html"]))
        template = env.get_template("test_report.html")

        html_text = template.render(html_dict)

        with open(self.html_path, "w") as f:
            f.write(html_text)
