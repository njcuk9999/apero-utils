"""
Test abstract class and factory for APERO.

@author: vandal
"""
import os
import glob
from datetime import datetime
from abc import ABC, abstractmethod
from typing import Optional, Union, List

import pandas as pd
from astropy.io import fits

import apero_tests_func as atf
from apero.core import constants


class Test(ABC):
    """Test."""

    def __init__(self, inst: str = 'SPIROU', setup: Optional[str] = None):
        """__init__.

        :param inst: Instrument name used in APERO (Default is SPIROU).
        :type inst: str
        :param setup: Path to apero setup. Using DRS_UCONFIG if None.
        :type setup: Optional[str]
        """
        self._params = constants.load(inst)

        # setup path
        if setup is None:
            self._setup = os.environ['DRS_UCONFIG']
        else:
            self._setup = setup

        self._instrument = self._params['INSTRUMENT']  # instrument

        self._date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # test date

    # =========================================================================
    # Properties
    # =========================================================================
    @property
    def setup(self):
        """setup."""
        return self._setup

    @property
    def instrument(self):
        """instrument."""
        return self._instrument

    @property
    def date(self):
        """date."""
        return self._date

    # =========================================================================
    # Abstract methods common to all tests
    # =========================================================================
    @property
    @abstractmethod
    def name(self):
        """Test full unique name."""

    @property
    @abstractmethod
    def recipe(self):
        """Recipe checked by the tests"""

    # @property
    # @abstractmethod
    # def outdict(self):
    #     """outdict."""

    @abstractmethod
    def runtest(self):
        """Method called to run the tests on APERO products"""

    # @abstractmethod
    # def htmlgen(self):
    #     """htmlgen."""


class CalibTest(Test):
    """CalibTest."""

    def __init__(self,
                 inst: str = 'SPIROU',
                 setup: Optional[str] = None,
                 logdir: Union[str, list] = 'night'):
        """__init__.

        :param inst:
        :type inst: str
        :param setup:
        :type setup: Optional[str]
        :param logdir:
        :type logdir: str
        """

        super().__init__(inst=inst, setup=setup)

        # Paths without end separator if any
        sep = os.path.sep
        self._reduced_path = self._params['DRS_DATA_REDUC'].rstrip(sep)
        self._calibdb_path = self._params['DRS_CALIB_DB'].rstrip(sep)

        # List of all reduced nights
        self._reduced_nights = atf.list_nights(self.reduced_path)

        # Where logs can be found
        if logdir == 'night':
            self._logdirs = self.reduced_nights.copy()
        elif isinstance(logdir, str):
            self._logdirs = [logdir]
        elif isinstance(logdir, list):
            self._logdirs = logdir
        else:
            raise TypeError('logdir must be a list or a string')

    @property
    def reduced_path(self) -> str:
        """Path to reduced directory

        :return: reduced_path
        :rtype: str
        """
        return self._reduced_path

    @property
    def calibdb_path(self) -> str:
        """Path do calibDB directory.

        :return: calibdb_path
        :rtype: str
        """
        return self._calibdb_path

    @property
    def reduced_nights(self) -> List[str]:
        """List of reduced night directories

        :return: reduced_nights
        :rtype: List[str]
        """
        return self._reduced_nights

    @property
    def log_df(self) -> pd.DataFrame:
        """Dataframe with log information

        :return: log_df
        :rtype: pd.DataFrame
        """
        return self._log_df

    @property
    def logdirs(self) -> List[str]:
        """List of directories where log files are expected for this test

        Currently supported:
            - single dir (e.g. other or master ight)
            - One log dir per night (typically the same as night)

        :return: logdirs
        :rtype: List[str]
        """
        return self._logdirs

    @property
    @abstractmethod
    def output_files(self) -> pd.Series:
        """output_files.

        :return: output_files
        :rtype: pd.Series
        """

    @property
    @abstractmethod
    def recipe(self) -> str:
        """Name of the reof the recipe"""

    def _gen_output_files(self) -> pd.Series:
        """Get output files
        Organize output files by output type and night using two indexing
        levels in pandas

        :return: Output series of all output files. Each file has two indices
                 (output pattern, night date)
        :rtype: pd.Series
        """

        inds = pd.MultiIndex.from_product(  # multi-index output and night
                [self.output_list, self.reduced_nights],
                names=['output', 'night'],
                )
        sep = os.path.sep
        paths = pd.Series(self.reduced_path  # Series of path per night+output
                          + sep + inds.get_level_values('night')
                          + sep + inds.get_level_values('output'),
                          index=inds
                          )
        files = paths.apply(glob.glob)  # Get file list for each pattern
        files = files.explode()         # Pop lists in individual cells

        return files

    def _gen_log_df(self):
        """_gen_log_df.

        Generate a dataframe with log for each directory, with two index
        levels:
            - DIRECTORY
            - RECIPE
        """
        # Get all logs in a dataframe
        allpaths = [os.path.join(self.reduced_path, ld, 'log.fits')
                    for ld in self.logdirs]
        paths = [p for p in allpaths if os.path.isfile(p)]
        dfs = [pd.DataFrame(fits.getdata(f)) for f in paths]
        log_df = pd.concat(dfs).set_index(['DIRECTORY'])
        log_df = log_df.loc[log_df.RECIPE == self.recipe]
        self._log_df = log_df

        # Store missing logs in another list
        self._missing_logs = [p for p in allpaths if not os.path.isfile(p)]

    def do_stop(self):
        """Do stop comparison for two inputs"""

    def get_outputs(self):
        """Get list of output files for each directory in reduced_nights"""

    def check_logs(self):
        """Check logs in directories, either by night, master, or other"""

    # =========================================================================
    # Abstract methods common to all calib tests
    # =========================================================================
    @property
    @abstractmethod
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: List[str]
        """

    @property
    @abstractmethod
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: List[str]
        """
