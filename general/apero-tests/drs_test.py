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

# TODO: Maybe we can move this so it's accessible to all files, maybe parameter
# in drs later
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


def load_fits_df(pathlist: List[str]) -> pd.DataFrame:
    """Load contents from list of fits tables in a dataframe

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
        self.fibers = []  # NOTE: not sure still needed. (output file def.?)
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
                print("WARNING: Overwriting kwarg instrument with recipe info")
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

        # Load log and index in a DataFrame
            self.log_df, _ = self._load_log()
            self.ind_df, _ = self._load_index()

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
            os.path.join(self.output_path, d, "log.fits")
            for d in os.listdir(self.output_path)
        ]
        log_df = load_fits_df(allpaths)

        # Use DIRECTORY as index and keep only relevant entries
        # whith self.recipe or extract recipes called
        log_df = log_df.set_index(["DIRECTORY"])

        # Important for extract recipes: keep only the recipe under test
        log_df = log_df[log_df.LOGFILE.str.contains(self.recipe_name)]

        # Paths without files
        missing_logs = [p for p in allpaths if not os.path.isfile(p)]

        return log_df, missing_logs

    def _load_index(self) -> Tuple[pd.DataFrame, List[str]]:
        """Get index contents in a dataframe

        Parse all index.fits files and return a dataframe.

        :returns: dataframe of index content and list of missing index files
        :rtype: Tuple[pd.DataFrame, list]
        """
        # Get all index.fits in a dataframe
        allpaths = [
            os.path.join(self.output_path, d, "index.fits")
            for d in os.listdir(self.output_path)
        ]
        ind_df = load_fits_df(allpaths)

        # Use NIGHTNAME as index and keep only relevant entries
        # whith self.recipe or extract recipes called
        ind_df = ind_df.set_index(["NIGHTNAME"])


        # Paths without files
        missing_inds = [p for p in allpaths if not os.path.isfile(p)]

        return ind_df, missing_logs


    def _load_header(self) -> pd.DataFrame: 
        # WIP

        # Search for all .fits files on disk
        filepaths = [d for d in glob.glob(self.params['DRS_DATA_REDUC']+'*/*.fits') 
                if not (os.path.basename(d).startswith('index') or 
                os.path.basename(d).startswith('log'))]



        pconstant = constants.pload(self.params['INSTRUMENT'])
        # list of index.fits columns
        col_index = pconstant.OUTPUT_FILE_HEADER_KEYS()

        # From apero drs

        for output in outputs:
        # get absfilename
        absoutput = os.path.join(opath, output)

        if not os.path.exists(absoutput):
            mtime = np.nan
        else:
            mtime = os.path.getmtime(absoutput)

        # get filename
        if 'FILENAME' not in col_index:
            col_index['FILENAME'] = [output]
            col_index['NIGHTNAME'] = [nightname]
            col_index['LAST_MODIFIED'] = [mtime]
        else:
            col_index['FILENAME'].append(output)
            col_index['NIGHTNAME'].append(nightname)
            col_index['LAST_MODIFIED'].append(mtime)

        # loop around index columns and add outputs to istore
        for icol in icolumns:
            # get value from outputs
            if icol not in outputs[output]:
                value = 'None'
            else:
                value = outputs[output][icol]
            # push in to istore
            if icol not in istore:
                istore[icol] = [value]
            else:
                istore[icol].append(value)



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
