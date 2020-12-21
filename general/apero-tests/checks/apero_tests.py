"""
Test abstract class and factory for APERO.

@author: vandal
"""
import os
import glob
from datetime import datetime
from abc import ABC, abstractmethod
from typing import Optional, Union

import numpy as np
import pandas as pd

import apero_tests_func as atf
from logobjects import Log
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
    def outdict(self):
        """outdict."""

    @abstractmethod
    def runtest(self):
        """Method called to run the tests on APERO products"""

    @abstractmethod
    def htmlgen(self):
        """htmlgen."""


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

        :return: __init__
        """

        super().__init__(inst=inst, setup=setup)

        # Paths without end separator if any
        sep = os.path.sep
        self._reduced_path = self._params['DRS_DATA_REDUC'].rstrip(sep)
        self._calib_db_path = self._params['DRS_CALIB_REDUC'].rstrip(sep)

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

        # Log dataframe
        self._log_df = self._get_log_df()

    @property
    def reduced_path(self) -> str:
        """Path to reduced directory

        :return: reduced_path
        :rtype: str
        """
        return self._reduced_path

    @property
    def calib_db_path(self) -> str:
        """Path do calibDB directory.

        :return: calib_db_path
        :rtype: str
        """
        return self._calib_db_path

    @property
    def reduced_nights(self) -> list[str]:
        """List of reduced night directories

        :return: reduced_nights
        :rtype: list[str]
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
    def logdirs(self) -> list[str]:
        """List of directories where log files are expected for this test

        Currently supported:
            - single dir (e.g. other or master ight)
            - One log dir per night (typically the same as night)

        :return: logdirs
        :rtype: list[str]
        """
        return self._logdirs

    def _get_output_df(self) -> pd.DataFrame:
        """Get output df
        Generate dataframe of output nights, and a full column with all files.

        :return: Output dataframe with output patterns as indices and
                 nights as columns
        :rtype: pd.DataFrame
        """
        df = pd.DataFrame([],
                          index=self.output_list,
                          columns=self.reduced_nights)
        df = df.apply(
                lambda x: self.reduced_path + os.path.sep + x.name
                                            + os.path.sep + x.index)
        df = df.applymap(glob.glob)  # Get list of files per night
        df['full'] = df.values.tolist()
        df.full = df.full.apply(np.concatenate)  # Full list of files per row

        return df

    def _get_log_df(self) -> pd.Dataframe:
        """_get_log_df.

        Generate DataFrame of logobjects.Log with log infomation.
        Initialized with only one column containing Log objects for each night.

        :return: log_df
        :rtype: pd.DataFrame
        """
        # TODO: Will need safety check for the missing logs...
        inds = pd.Index(self.logdirs)
        logpaths = (self.reduced_path
                    + os.path.sep + inds
                    + os.path.sep + 'log.fits')
        logs = pd.Series(logpaths, index=inds)
        logs = logs.apply(Log)
        log_df = pd.Dataframe(logs, index=inds, columns=['log'])

        return log_df

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
    def output_list(self) -> list[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """

    @property
    @abstractmethod
    def calibdb_list(self) -> list[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
    @property
    @abstractmethod
    def output_df(self) -> pd.DataFrame:
        """Dataframe with output files per night and cumulated

        :return: output_df
        :rtype: pd.DataFrame
        """
        return self._output_df
