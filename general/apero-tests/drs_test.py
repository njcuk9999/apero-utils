"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import os
from datetime import datetime
from typing import List, Optional, Tuple

import pandas as pd
from astropy.table import Table
from jinja2 import Environment, FileSystemLoader, select_autoescape

# from . import OUTDIR, TEMPLATEDIR
PARENTDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
TEMPLATEDIR = os.path.join(PARENTDIR, "templates")
OUTDIR = os.path.join(PARENTDIR, "out")


def removext(name: str, ext: str = ".py"):
    """Remove extension of a file/recipe name

    :param name: file/recipe name
    :type name: str

    :returns: cleaned up name
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
        drs_recipe: Optional[str] = None,
        setup: Optional[str] = None,
        testnum: int = 1,
    ):

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
        self.params = None  # DRS parameters
        self.ismaster = False  # Is it a master recipe?
        self.fibers = []  # NOTE: not sure still needed
        self.output_hkeys = []  # Output header keys
        self.calibdb_keys = []  # Calibdb keys
        self.calib_hkeys = []  # All calib header keys
        self.calls_extract = False  # NOTE: not sure still needed
        self.input_path = ""  # Path to input dir
        self.output_path = ""  # Path to output dir

        # Overwrite some parameters with recipe info automatically
        if drs_recipe is not None:
            # The recipe object and test ID
            self.recipe = drs_recipe
            if self.instrument is not None:
                print("WARNING: Overwriding kwarg instrument with recipe info")
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

        # TODO: Load log, outputs and calibdb in a nice way with APERO
        self.log_df, _ = self._load_log()

    # =========================================================================
    # Functions to load output information
    # =========================================================================
    def _load_log(self) -> Tuple[pd.DataFrame, List[str]]:
        """Get log content in a dataframe

        Parse all log files and return a dataframe with only entries that have
        the current recipe in LOGFILE.

        :returns: dataframe of log content and list of missing log files
        :rtype: Tuple[pd.DataFrame, list]
        """
        # Get all logs in a dataframe
        allpaths = [
            os.path.join(self.output_path, ld, "log.fits")
            for ld in os.listdir(self.output_path)
        ]
        paths = [p for p in allpaths if os.path.isfile(p)]
        dfs = [Table.read(f).to_pandas() for f in paths]
        log_df = pd.concat(dfs)

        # Decode strings (fits from astropy are encoded)
        for col in log_df.columns:
            # try/except for non-string columns
            try:
                log_df[col] = log_df[col].str.decode("utf-8")
            except AttributeError:
                pass

        # Use DIRECTORY as index and keep only relevant entries
        # whith self.recipe or extract recipes called
        log_df = log_df.set_index(["DIRECTORY"])

        # Important for extract recipes: keep only the recipe under test
        log_df = log_df[log_df.LOGFILE.str.contains(self.recipe)]

        # Paths without files
        missing_logs = [p for p in allpaths if not os.path.isfile(p)]

        return log_df, missing_logs

    # =========================================================================
    # Utility functions
    # =========================================================================
    def get_dir_path(self):
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
