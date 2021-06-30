from typing import Any, List, Optional, Sequence

import pandas as pd
from pandas.core.frame import DataFrame

import apero_tests.utils as ut


class SubTest:
    # TODO: result should probably more restricted,
    # probably Union[Mapping, Sequence]
    # TODO: Restrict details type/formatting
    def __init__(self,
                 subtest_id: Optional[str] = None,
                 description: Optional[str] = None,
                 result: Any = None,
                 comments: Optional[str] = None,
                 details: Optional[str] = None):
        """
        Individual SubTest to test an APERO DRS recipe.

        To make a subtest "silent", set its result to None so that DrsTest
        ignores it.

        :param subtest_id: Subtest ID, set to lowercase class name if None,
                           defaults to None.
        :type subtest_id: Optional[str], optional
        :param description: Subtest description, defaults to None
        :type description: Optional[str], optional
        :param result: Result of the subtest, defaults to None
        :type result: Any, optional
        :param comments: Extra comments about subtest, defaults to None
        :type comments: Optional[str], optional
        :param details: Additional details about subtest (second html page),
                        defaults to None
        :type details: Any, optional
        """

        if subtest_id is None:
            self.id = type(self).__name__.lower()
        else:
            self.id = subtest_id
        if description is None:
            self.description = "Generic Drs Test Entry"
        else:
            self.description = description
        self.result = result
        self.comments = comments
        self.details = details

    def run(self):
        """
        Method that runs the subtests, then modifying its attributes.

        :raises NotImplementedError: By default, not implemented for the
                                     base SubTest class
        """
        # Using __name__ in case a children class also has no `run()`
        msg = (f"The run method is not implemented for {type(self).__name__}."
               " Create a children class or assign result directly.")
        raise NotImplementedError(msg)


class CountRawTest(SubTest):
    def __init__(self,
                 input_path: str):
        """
        Subtest that counts the number of raw files on disk.

        The result is the count.

        :param input_path: Path to raw file directory
        """
        self.input_path = input_path

        super().__init__(description=f"# of raw files in {self.input_path}")

    def run(self):

        exclude_list = ["index.fits", "log.fits"]
        raw_series = ut.get_output_files(self.input_path, exclude_fname=exclude_list)
        self.result = raw_series.size


class CountLogTest(SubTest):
    def __init__(self,
                 log_df: DataFrame,
                 master_flag: str = "--master",
                 group_kwds: Optional[Sequence[str]] = None):
        """
        Subtest that counts the number of calls in the log for a recipe.

        The result is a pandas Series with the counts for each unique
        combination of keys specified by `group_kwds`.

        :param log_df: DataFrame with log entries for a recipe.
        :type log_df: DataFrame
        :param master_flag: Flag that denotes a master call in the log
                            RUNSTRING, defaults to "--master"
        :type master_flag: str, optional
        :param group_kwds: Columns used to group the log DF before counting.
                           The default ["RECIPE", "SUBLEVEL", "LEVEL_CRIT"]
                           is used when set to None.
        :type group_kwds: Optional[Sequence[str]], optional
        """

        super().__init__(description="# of calls in logs")

        self.log_df = log_df
        self.master_flag = master_flag
        if group_kwds is not None:
            self.group_kwds = group_kwds
        else:
            self.group_kwds = ["RECIPE", "SUBLEVEL", "LEVEL_CRIT"]

    def run(self):

        master_mask = self.log_df.RUNSTRING.str.contains(self.master_flag)
        master_df = self.log_df[master_mask]
        master_count = master_df.groupby(self.group_kwds).count().KIND
        tot_count = self.log_df.groupby(self.group_kwds).count().KIND
        log_count = tot_count.sub(master_count, fill_value=0).astype(int)

        self.result = log_count


class CountOutTest(SubTest):
    def __init__(self,
                 ind_df: DataFrame,
                 output_path: str,
                 output_hkeys: List[str],
                 unique: bool = False,
                 group_kwds: Optional[Sequence[str]] = None):
        """
        Subtest that counts the number of output files for a recipe

        :param ind_df: Index dataframe from a DrsTest object
        :type ind_df: DataFrame
        :param output_path: Path to the output directory of the recipe
        :type output_path: str
        :param output_hkeys: Header keys for each output type, used for
        consistency
        :type output_hkeys: List[str]
        :param unique: Return unique count if true, defaults to False
        :type unique: bool, optional
        :param group_kwds: Keys used to group the index before counting,
                           The default ["KW_OUTPUT", "KW_DPRTYPE", "KW_FIBER"]
                           is used if None.
        :type group_kwds: Optional[Sequence[str]], optional
        """

        self.unique = unique
        self.output_path = output_path

        if self.unique:
            description = f"# of unique outputs in {self.output_path}"
        else:
            description = f"# of outputs in {self.output_path}"

        super().__init__(description=description)

        self.ind_df = ind_df
        self.output_hkeys = output_hkeys
        if group_kwds is not None:
            self.group_kwds = group_kwds
        else:
            # TODO: Should we use drs params to get keys here or is literal ok?
            self.group_kwds = ["KW_OUTPUT", "KW_DPRTYPE", "KW_FIBER"]

    def run(self):
        if not self.ind_df.empty:
            if self.unique:
                self.result = self.ind_df.groupby(
                    self.group_kwds).FILENAME.nunique()
            else:
                self.result = self.ind_df.groupby(
                    self.group_kwds).FILENAME.count()
        else:
            self.result = pd.Series(0,
                                    index=self.output_hkeys,
                                    name="FILENAME")


class CountQCTest(SubTest):
    def __init__(self, log_df: DataFrame, test_html_path: str):
        super().__init__(
            description="# of log entries that failed one or more QC")

        self.log_df = log_df
        self.test_html_path = test_html_path

    def run(self):

        log_qc_failed = self.log_df[~self.log_df.PASSED_ALL_QC]
        self.result = len(log_qc_failed)

        odometer_flag = 'ODOMETER' in log_qc_failed.columns

        if self.result > 0 and odometer_flag:
            self.comments = 'One or more recipe have failed QC.'
            self
            log_reset = log_qc_failed.reset_index()
            data_dict_check_qc = {
                'Night': log_reset.DIRECTORY.values,
                'Odometer': log_qc_failed.ODOMETER.values,
                'QC_STRING': log_qc_failed.QC_STRING.values,
            }
            self.details = ut.inspect_table(
                self.test_html_path, self.id, data_dict_check_qc,
                'Odometers that Failed One or More Quality Control')

        if self.result > 0 and not odometer_flag:
            self.comments = 'One or more recipe have failed QC.'
            self
            log_reset = log_qc_failed.reset_index()
            data_dict_check_qc = {
                'Night': log_reset.DIRECTORY.values,
                'QC_STRING': log_qc_failed.QC_STRING.values,
            }
            self.details = ut.inspect_table(
                self.test_html_path, self.id, data_dict_check_qc,
                'Nights that Failed Quality Control')


class PlotQCTest(SubTest):
    def __init__(self, log_df: DataFrame, test_html_path: str,
                 recipe_name: str):
        super().__init__(
            description="Plot of the various QCs as a function of odometer, order or time")

        self.log_df = log_df
        self.test_html_path = test_html_path
        self.recipe_name = recipe_name

    def run(self):
        qc_names = self.log_df.QC_NAMES.str.split(r'\|\|', expand=True)
        qc_names = qc_names.iloc[0]  # Keep only one row
        qc_values = self.log_df.QC_VALUES.str.split(r'\|\|', expand=True)
        qc_values.columns = qc_names
        try:
            # NOTE: .convert_dtypes will do in pd versions >= 1.0.0
            for col in qc_values.columns:
                try:
                    qc_values[col] = qc_values[col].astype(float)
                except ValueError:
                    del qc_values[col]
            if len(qc_values.columns) == 0:
                self.result = None
                return
        except KeyError:
            pass

        odometer_flag = 'ODOMETER' in self.log_df.columns
        order_flag = len(qc_names) == 49

        if odometer_flag:
            data_dict_qc_plot = {'Night': qc_values.index.tolist()}
            data_dict_qc_plot['Odometer'] = self.log_df.ODOMETER.tolist()
            for key, series in qc_values.iteritems():
                data_dict_qc_plot[key] = series.tolist()

        elif order_flag:
            data_dict_qc_plot = {'Order': list(range(1, 50)),
                                   qc_names[0]: qc_values.values[0]}

        else:
            data_dict_qc_plot = {'Night': qc_values.index.tolist()}
            for key, series in qc_values.iteritems():
                data_dict_qc_plot[key] = series.tolist()

        self.result = "See the Details column"
        self.details = ut.inspect_plot(
            self.test_html_path, self.id, data_dict_qc_plot,
            f'{self.recipe_name}.py Quality Control')


class CountEndedTest(SubTest):
    # TODO: Very similar to QC, could probably put most of code together
    # and change a few args
    def __init__(self, log_df: DataFrame, test_html_path: str):
        super().__init__(description="# of log entries that failed to finish")

        self.log_df = log_df
        self.test_html_path = test_html_path

    def run(self):

        log_ended_false = self.log_df[~self.log_df.ENDED]
        self.result = len(log_ended_false)

        if self.result > 0:
            self.comments = 'One or more recipe have failed QC.'
            self
            log_reset = log_ended_false.reset_index()
            data_dict_check_ended = {
                'Night': log_reset.DIRECTORY.values,
                'ERRORS': log_ended_false.ERRORS.values,
                'LOGFILE': log_ended_false.LOGFILE.values,
            }
            self.details = ut.inspect_table(self.test_html_path, self.id,
                                            data_dict_check_ended,
                                            'Nights that Failed to finish')


class CountCalibEntries(SubTest):
    def __init__(self, calib_df: DataFrame, calib_keys: List[str]):

        super().__init__(description="# of entries in calib DB")

        self.calib_df = calib_df
        self.calib_keys = calib_keys

    def run(self):
        zero_count = pd.Series(0, index=self.calib_keys)

        if not self.calib_df.empty:
            master_mask = self.calib_df["master"].astype(bool)
            if master_mask.sum() == 0:
                master_calib_count = zero_count.copy()
            else:
                master_calib_df = self.calib_df[master_mask]
                master_calib_count = master_calib_df.groupby("key").size()
                self.comments = f"Additional entries with master == 1: {master_calib_count}"

            calib_count = self.calib_df.groupby("key").size()
            calib_count -= master_calib_count

        else:
            self.result = None
            return

        self.result = calib_count

class CountTelluEntries(SubTest):
    def __init__(self, tellu_df: DataFrame):

        super().__init__(description="# of entries in tellu DB")

        self.tellu_df = tellu_df

    def run(self):
        if not self.tellu_df.empty:
            tellu_count = self.tellu_df.groupby("key").size()

        else:
            self.result = None
            return

        self.result = tellu_count
