"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import os
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union

import apero_tests.subtest as st
import apero_tests.utils as ut
import pandas as pd
from apero.core import constants
from apero.core.core.drs_argument import DrsArgument
from apero.core.core.drs_base_classes import BinaryDict
from apero.core.core.drs_file import DrsFitsFile
from apero.core.instruments.spirou.recipe_definitions import recipes
from apero.core.utils.drs_recipe import DrsRecipe
from apero.tools.module.processing.drs_processing import skip_clean_arguments
from jinja2 import Environment, FileSystemLoader, select_autoescape
from pandas import DataFrame

PARENTDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEMPLATEDIR = os.path.join(PARENTDIR, "templates")
OUTDIR = os.path.join(PARENTDIR, "out")
CACHEDIR = os.path.join(PARENTDIR, "cache")

RECIPE_DICT = {x.name: x for x in recipes}


# TODO: Link to documentation instead of Github
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
        :pp_flag: True if test is for preprocessing recipe
        :type pp_flag: bool
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
        self.output_drs_files = None
        self.output_hkeys = []  # Output header keys
        self.calibdb_keys = []  # Calibdb keys
        self.input_path = ""  # Path to input dir
        self.output_path = ""  # Path to output dir
        self.calibdb_path = ""
        self.html_path = ""  # Path to HTML report
        self.log_df = None
        self.ind_df = None
        self.calib_df = None
        self.tellu_df = None
        self.cdb_used_df = None
        self.subtest_list = None

        # PP flag
        self.pp_flag = pp_flag

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
            # TODO: Confirm this attribute with Neil
            self.params = self.recipe.params
            self.pconst = constants.pload(self.instrument)
            self.ismaster = self.recipe.master
            # HACK: Should use DRS to define that. Remove when fixed in apero-drs
            if "wave_master" in self.recipe_name:
                self.ismaster = True

            # Path to input and output directories
            self.dirpaths = self.get_dir_path()
            # TODO: Confirm these attribute with Neil
            self.input_path = self.dirpaths[self.recipe.in_block_str]
            self.output_path = self.dirpaths[self.recipe.out_block_str]
            self.calibdb_path = self.params["DRS_CALIB_DB"]
            self.calibdb_seed_path = os.path.join(
                self.params["DRS_DATA_ASSETS"], "databases", "reset.calib.csv"
            )

            # Path to HTML report
            self.html_path = Path(
                OUTDIR, self.test_id, ".".join([self.test_id, "html"])
            )

            # Get output files for this recipe
            self.output_drs_files = self.get_recipe_outputs()

            self.output_hkeys = [
                o.required_header_keys["KW_OUTPUT"]
                for o in self.output_drs_files
                if "KW_OUTPUT" in o.required_header_keys
            ]

            # If out KW_OUTPUT, it will likely be "--" in index. Useful for PP
            if len(self.output_hkeys) == 0:
                # ???: Is there something in DRS that gives this null value str
                # Could use that instead of hardcoding
                self.output_hkeys = ["--"]

            # NOTE: index is filtered with log to keep only relevant entries.
            # A global test is used to report outputs that match no logs
            # (see test_definitions)
            self.log_df = self.load_log_df(all_log_df=all_log_df)
            self.ulog_df = self.log_df.drop_duplicates(
                subset=["RUNSTRING", "RECIPE", "SUBLEVEL", "LEVELCRIT", "PID"]
            ).copy()
            self.ind_df = self.load_ind_df(all_ind_df=all_index_df)

            self.model_ind_df = self.generate_model_index()

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

            self.calibdb_list = ut.load_db_list(
                "CALIB", instrument=self.instrument
            )

            # Load calibdb DF (master entries and calibdb files)
            self.calib_df = self.load_calib_df(
                all_calib_df=all_master_calib_df
            )

            # Needed for "used calibs" test
            self.all_calib_df = all_master_calib_df

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

        if self.recipe_name is None:
            raise ValueError("Cannot load index dataframe if no recipe is set")

        if all_ind_df is None:
            all_ind_df = ut.load_index_df(self.output_path)

        if self.log_df is None:
            raise AttributeError(
                "A non-null log_df is required to filter the index"
            )

        # We use the log to filter the index, also keep only output hkeys
        # FUTURE: Might have an indepent way that does not use log in v0.7
        ind_df = all_ind_df[all_ind_df.KW_PID.isin(self.log_df.PID)]

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
            return self.calib_df

        if all_calib_df is None:
            all_calib_df = ut.load_db_entries("CALIB")

        if self.calibdb_keys is None:
            msg = (
                "Cannot load calibdb dataframe if no calibdb keys are set,"
                " returning None"
            )
            warnings.warn(msg, RuntimeWarning)
            return None

        # NOTE: Copy is used here so that future assignemnt never affect
        #       original (suppressing SettingWithCopyWarning)
        # TODO: Check if need to use fiber as well when comparing with calibdb_keys
        # TODO: Make sure works with or without fibers
        calib_df = all_calib_df[
            all_calib_df["KEYNAME"].isin(self.calibdb_keys)
        ].copy()

        # Don't keep seed values
        calib_seed = pd.read_csv(self.calibdb_seed_path, skipinitialspace=True)
        calib_df = calib_df[~calib_df.UHASH.isin(calib_seed.UHASH)]

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
            all_tellu_df = ut.load_db_entries("TELLU")

        if self.telludb_keys is None:
            msg = (
                "Cannot load index dataframe if no telludb keys are set,"
                " returning None"
            )
            warnings.warn(msg, RuntimeWarning)
            return None

        # TODO: Check if need to use fiber as well when comparing with calibdb_keys
        # TODO: Make sure works with or without fibers
        tellu_df = all_tellu_df[
            all_tellu_df["KEYNAME"].isin(self.telludb_keys)
        ]

        return tellu_df

    def load_cdb_used_df(
        self, all_cdb_used_df: Optional[DataFrame] = None, force: bool = True
    ) -> Union[DataFrame, None]:
        """
        Load dataframe with used calibrations for each files. Calls
        utils.get_cdb_df() if all_cdb_used_df not passed.

        Dataframe has 2 level column
        (CDB type on first leve, calib file and DELTA_MJD on second)

        :param all_cdb_used_df: Dataframe with all used calibrations for
                                reduction, defaults to None
        :type all_cdb_used_df: Optional[DataFrame], optional
        :param force: Re-load even if already set as attribute,
                      defaults to True
        :type force: bool, optional
        :return: DataFrame with used calibration info
        :rtype: Union[DataFrame, None]
        """
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
            ind_night_file = self.ind_df.reset_index()[["OBS_DIR", "FILENAME"]]
            keep_ind = pd.MultiIndex.from_frame(ind_night_file)
            cdb_used_df = all_cdb_used_df.loc[keep_ind]
        else:
            cdb_used_df = all_cdb_used_df.iloc[0:0]

        return cdb_used_df

    # =========================================================================
    # Utility functions
    # =========================================================================
    def get_recipe_outputs(
        self, recipe: Optional[Union[DrsRecipe, str]] = None
    ):
        recipe = recipe or self.recipe

        # Try to get recipe from name if we have a string
        if isinstance(recipe, str):
            try:
                recipe = RECIPE_DICT[recipe]
            except AttributeError:
                recipe = RECIPE_DICT[recipe + ".py"]

        odict = recipe.outputs
        outfiles = [
            o for o in list(odict.values()) if isinstance(o, DrsFitsFile)
        ]

        return outfiles

    def get_dir_path(self) -> Dict:
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
        """
        Add CALIB_KEY with corresponding calibdb key for each output in index
        """

        # Get all possible fibers
        sci_fibers, ref_fiber = self.pconst.FIBER_KINDS()
        all_fibers = sci_fibers + [ref_fiber]

        fiber_series = "_" + self.ind_df.KW_FIBER
        fiber_series = fiber_series.where(
            fiber_series.isin([f"_{fib}" for fib in all_fibers]), ""
        )
        # out_fiber_series = self.ind_df.KW_OUTPUT + fiber_series
        out_fiber_series = self.ind_df.KW_OUTPUT.copy()

        # TODO: Check the two columns instead of merging here
        # TODO: out_fiber_series has LOC_ORDERP_AB and out_to_calib maps from LOC_ORDERP only
        # because now the key is is a separate column
        # - First, just don't use fibers and see if counts match
        # - We used fibers before, so should we do the extra work to add them back ?
        # - Actually fibers are needed to make count right (could have twice as much for example)
        out_to_calib = {
            (
                v + f"_{k.split('_')[-1]}"
                if k.split("_")[-1] in all_fibers
                else v
            ): k
            for k, v in self.calibdb_to_output.items()
        }

        # If map OK, this will do
        self.ind_df["CALIB_KEY"] = out_fiber_series.map(out_to_calib)

    def generate_model_index(self) -> pd.DataFrame:
        """
        Generate expected index based on the log dataframe
        """

        # Use unique log, exact same row generatees exact same file (duh)
        log_df = self.ulog_df.copy()

        # Drop rows with same runstring (the generate the same output)
        # Before that, remove some args from runstring if needed
        log_df.RUNSTRING = log_df.RUNSTRING.apply(skip_clean_arguments)
        log_df = log_df[
            ~log_df.duplicated(
                subset=["RUNSTRING", "RECIPE", "SUBLEVEL", "LEVELCRIT"]
            )
        ]

        # Drop runs with "--master" because they re-generate the same files
        # NOTE: Might be different in 0.7
        # TODO: Use apero function to remove unwanted arguments regardless of if recipe is master
        master_mask = log_df.RUNSTRING.str.contains("--master")
        log_df = log_df[~master_mask] if not self.ismaster else log_df
        log_df = log_df.copy()

        # Add binary flags from recipe
        # This ensures that all keys are used.
        # NOTE: The default value None won't work for boolean masks
        # but default value depends on what we check. So need to set these before
        # perform action (see quicklook block for example)
        # TODO: Set default values automatically, in DRS or in a separate default dict ?
        flag_names = list(log_df["FLAGSTR"].str.split("|").explode().unique())
        log_df[
            flag_names
        ] = None  # None by default in case some are undefined for some rows
        for i, (_, row) in enumerate(log_df.iterrows()):
            flags = BinaryDict()
            flags.add_keys(row["FLAGSTR"].split("|"))
            flags.encode(row["FLAGNUM"])
            log_df.iloc[i, log_df.columns.isin(flag_names)] = [
                flags[k] if k in flags else None for k in flag_names
            ]

        if "QUICKLOOK" in log_df.columns:
            # Don't drop files without quicklook flag (set to false for later in loop)
            log_df["QUICKLOOK"] = log_df["QUICKLOOK"].fillna(False)

        # Get all recipes in the log
        recipes = list(log_df.RECIPE.unique())
        log_recipe_dict = dict()

        # Store recipe info in a dict
        for rname in recipes:
            recipe_obj = RECIPE_DICT[rname + ".py"]
            out_files = self.get_recipe_outputs(recipe=recipe_obj)
            fibers = [f.fibers for f in out_files]

            rdict = {
                "name": rname,
                "recipe_obj": RECIPE_DICT[rname + ".py"],
                "out_files": self.get_recipe_outputs(recipe=recipe_obj),
                "fibers": fibers,
            }

            log_recipe_dict[rname] = rdict

        def _get_dprtype(
            ofile: Union[DrsFitsFile, List[DrsFitsFile], None]
        ) -> str:
            """
            Get DPR Type of file recursively

            :param ofile: Output DRS File for which we want DPR_TYPE
            :type ofile: Union[DrsFitsFile, List[DrsFitsFile], None]
            """
            if isinstance(ofile.intype, DrsFitsFile):
                # Go up intype until get DPRTYPE. Should end with None or value
                if "KW_DPRTYPE" in ofile.intype.required_header_keys:
                    return ofile.intype.required_header_keys["KW_DPRTYPE"]
                else:
                    return _get_dprtype(ofile.intype)
            # If a list or none, just return None
            # FUTURE: Probably not great, maybe update in 0.7 ? (Talk to Neil)
            else:
                return None

        # Most recipes have multiple levels that correspond to a specific loop
        # LEAKM has two loops in the same levels, and we care only about fiber for the
        # leakm files
        if "leak_master" in self.recipe_name:
            dflist = []
            for _, pdf in log_df.groupby(["PID", "LEVEL"], sort=False):
                first_crit = pdf["LEVELCRIT"].str.split("=").str[0]
                # Want same number of each
                if first_crit.nunique() == 1:
                    dflist.append(pdf)
                else:
                    # NOTE: Alternative way to get last loop in level, kept in case we need
                    # it in the future
                    # used_crit = pdf.groupby(first_crit)["SUBLEVEL"].min().idxmax()
                    used_crit = "fiber"

                    dflist.append(pdf[first_crit == used_crit])

            log_df = pd.concat(dflist)

        rows = []
        row_inds = []
        for i, (log_night, log_row) in enumerate(log_df.iterrows()):

            recipe_dict = log_recipe_dict[log_row.RECIPE]

            # For recipe under test, will need to check calls to other recipes
            # NOTE: Important assumption that not 2 or more chained calls
            # NOTE: Shoud we use RECIPE_KIND here?
            # The assumption should be OK (for now) according to Neil
            if log_row.RECIPE == self.recipe_name:
                other_dict = {
                    k: v
                    for k, v in log_recipe_dict.items()
                    if k != log_row.RECIPE
                }
                other_files_list = [
                    r["out_files"] for r in other_dict.values()
                ]
                other_files = [f.name for fl in other_files_list for f in fl]
            else:
                other_files = []

            if (
                log_row["LEVELCRIT"] is not None
                and "fiber" in log_row["LEVELCRIT"]
            ):
                # Get fiber from level_crit to avoid repeating entries if
                # recipe has more than one fiber in output files but this log
                # entry only generates one fiber
                # NOTE: Fiber must be a list of lists
                # (list of fibers for each output type)
                fibers = [
                    [
                        log_row["LEVELCRIT"]
                        .split("fiber")[-1]
                        .split("=")[1]
                        .strip()
                    ]
                ] * len(recipe_dict["out_files"])
            else:
                # If no fiber from level crit, use values given by recipe
                fibers = recipe_dict["fibers"]

            # Generate index entry for each DRS output file of recipe
            for j, ofile in enumerate(recipe_dict["out_files"]):

                # We exclude DEBUG from input
                if ofile.outclass.debug:
                    continue

                # Handle fiber list being None
                if "SCIFIBER" in log_df.columns:
                    if log_row["SCIFIBER"]:
                        fiber_list = ["AB"]
                    else:
                        fiber_list = ["C"]
                else:
                    fiber_list = fibers[j] or ["--"]

                # Set KW_OUTPUT
                kw_output = (
                    ofile.name
                    if "KW_OUTPUT" in ofile.required_header_keys
                    else "--"
                )

                # Skip files that are not outputs of the recipe
                if kw_output not in self.output_hkeys:
                    continue

                # If an output file is in this recipe and in another from the
                # same log (only useful for "main" recipe of this test)
                # For e.g., an extract row will always have other_files = []
                # (see def above)
                if kw_output in other_files:
                    continue

                # HACK: Ideally this would not rely on QL_ prefix to skip ext files I think
                if (
                    "QUICKLOOK" in log_row.keys()
                    and not log_row["QUICKLOOK"]
                    and kw_output.startswith("QL_")
                ):
                    continue
                elif (
                    "QUICKLOOK" in log_row.keys()
                    and log_row["QUICKLOOK"]
                    and not kw_output.startswith("QL_")
                ):
                    continue

                # HACK: Ideally this would not rely on EXT_ prefix to skip ext files I think
                if (
                    self.ismaster
                    and not ofile.outclass.master
                    and not ofile.name.startswith("EXT_")
                ):
                    if "wave_master" not in self.recipe_name:
                        print(
                            f"UNEXPECTED file {ofile} for recipe {self.recipe_name}"
                        )
                    continue

                dprtype = _get_dprtype(ofile)
                # HACK: Will have a better way to do this in 0.7 (talk to Neil)
                if dprtype is not None and log_row.ARGS.endswith("]"):
                    # If we can match a DPR type for the call in log and it's
                    # not the right one, we skip
                    if not log_row.ARGS.endswith(f"[{dprtype}]"):
                        continue

                for k, fiber in enumerate(fiber_list):
                    fake_fname = log_night + f"_{i}_{j}_{k}"
                    kw_pid = log_row.PID
                    kw_fiber = fiber

                    row_dict = {
                        "FILENAME": fake_fname,  # Required for counting
                        "KW_OUTPUT": kw_output,
                        "KW_PID": kw_pid,
                        "KW_FIBER": kw_fiber,
                    }

                    rows.append(row_dict)
                    row_inds.append(log_night)

        return pd.DataFrame(rows, index=row_inds)

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
        subtest_list.append(st.CountLogTest(self.log_df, self.html_path, master=self.ismaster))

        # Count all outputs
        subtest_list.append(
            st.CountOutTest(self.ind_df, self.output_path, self.output_hkeys)
        )

        # Count unique outputs
        st_unique_outputs = st.CountOutTest(
            self.ind_df, self.output_path, self.output_hkeys, unique=True
        )
        subtest_list.append(st_unique_outputs)

        # Count unique outputs in model index
        st_unique_model_outputs = st.CountOutTest(
            self.model_ind_df,
            self.output_path,
            self.output_hkeys,
            unique=True,
            description="# of unique outputs in model index df",
        )
        subtest_list.append(st_unique_model_outputs)

        subtest_list.append(
            st.ComparisonTest(st_unique_outputs, st_unique_model_outputs)
        )

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
                self.calib_df, self.calibdb_keys, master=self.ismaster
            )
            subtest_list.append(st_cdb_entries)

            subtest_list.append(
                st.ComparisonTest(st_index_calib, st_cdb_entries)
            )

            subtest_list.append(
                st.CheckIndexCalibFiles(
                    self.ind_df, self.calibdb_list, self.html_path
                )
            )

            subtest_list.append(
                st.CheckCalibEntriesFiles(
                    self.calib_df, self.calibdb_list, self.html_path
                )
            )

            subtest_list.append(
                st.CheckUsedCalibs(
                    self.cdb_used_df, self.html_path, self.all_calib_df
                )
            )

        # telluDB output count
        # TODO: Uncomment when tellu df filtering is fixed
        # subtest_list.append(st.CountTelluEntries(self.tellu_df))

        self.subtest_list = subtest_list

    def run_test(self):
        """
        Run all subtests and generate HTML summary for the recipe
        """

        # Keep this here so that if subtests access it via attribute it is set
        if self.html_path == "":
            # Path to HTML report
            self.html_path = Path(
                OUTDIR, self.test_id, ".".join([self.test_id, "html"])
            )

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
