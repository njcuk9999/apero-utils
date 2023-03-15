import numbers
from typing import Any, List, Optional, Sequence

import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame

import apero_tests.utils as ut
from apero_tests.display import delta_mjd_plot, inspect_plot, inspect_table

COMPARISONS = {
    "eq": np.equal,
    "gt": np.greater,
    "lt": np.less,
}


class SubTest:
    def __init__(
        self,
        subtest_id: Optional[str] = None,
        description: Optional[str] = None,
        result: Any = None,
        comments: Optional[str] = None,
        details: Optional[str] = None,
        color: Optional[str] = None,
    ):
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
        :type details: Optional[str], optional
        :param color: Color to use for the cell background in test output,
                      useful to show if passed or fail when applicable
        :type color: Optional[str], optional
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
        self.color = color

    def run(self):
        """
        Method that runs the subtests, then modifying its attributes.

        :raises NotImplementedError: By default, not implemented for the
                                     base SubTest class
        """
        # Using __name__ in case a children class also has no `run()`
        msg = (
            f"The run method is not implemented for {type(self).__name__}."
            " Create a child class or assign result directly."
        )
        raise NotImplementedError(msg)


class CountRawTest(SubTest):
    def __init__(self, input_path: str):
        """
        Subtest that counts the number of raw files on disk.

        The result is the count.

        :param input_path: Path to raw file directory
        """
        self.input_path = input_path

        super().__init__(description=f"# of raw files in {self.input_path}")

    def run(self):

        exclude_list = ["index.fits", "log.fits"]
        raw_series = ut.get_output_files(
            self.input_path, exclude_fname=exclude_list
        )
        self.result = raw_series.size


class CountLogTest(SubTest):
    def __init__(
        self,
        log_df: DataFrame,
        test_html_path: str,
        ref_flag: str = "--ref",
        group_kwds: Optional[Sequence[str]] = None,
    ):
        """
        Subtest that counts the number of calls in the log for a recipe.

        The result is a pandas Series with the counts for each unique
        combination of keys specified by `group_kwds`.

        :param log_df: DataFrame with log entries for a recipe.
        :type log_df: DataFrame
        :param ref_flag: Flag that denotes a ref call in the log
                            RUNSTRING, defaults to "--ref"
        :type ref_flag: str, optional
        :param group_kwds: Columns used to group the log DF before counting.
                           The default ["RECIPE", "SUBLEVEL", "LEVELCRIT"]
                           is used when set to None.
        :type group_kwds: Optional[Sequence[str]], optional
        """

        super().__init__(description="# of calls in logs")

        self.log_df = log_df
        self.ref_flag = ref_flag
        self.test_html_path = test_html_path
        if group_kwds is not None:
            self.group_kwds = group_kwds
        else:
            self.group_kwds = ["RECIPE", "SUBLEVEL", "LEVELCRIT"]

    def run(self):

        ref_mask = self.log_df.RUNSTRING.str.contains(self.ref_flag)
        ref_df = self.log_df[ref_mask]
        ref_count = ref_df.groupby(self.group_kwds).count().RECIPE_KIND
        tot_count = self.log_df.groupby(self.group_kwds).count().RECIPE_KIND
        log_count = tot_count.sub(ref_count, fill_value=0).astype(int)

        dup_subset = ["RUNSTRING", "RECIPE", "SUBLEVEL", "LEVELCRIT", "PID"]
        ulog_df = self.log_df.drop_duplicates(subset=dup_subset).copy()

        if ulog_df.equals(self.log_df):
            self.color = "Lime"
        else:
            self.comments = "Some entries are duplicated"
            self.color = "Yellow"
            ulog_df = self.log_df
            # Keep one duplicated entry
            # Keep one duplicated entry
            dup_log_df = self.log_df[
                self.log_df.duplicated(subset=dup_subset)
            ].drop_duplicates(subset=dup_subset)

            dup_dict = dup_log_df.to_dict(orient="list")

            self.details = inspect_table(
                self.test_html_path,
                self.id,
                dup_dict,
                "Duplicated log entries",
            )

        self.result = log_count


class ComparisonTest(SubTest):
    def __init__(
        self,
        subtest1: SubTest,
        subtest2: SubTest,
        passed: str = "eq",
        conditional: Optional[str] = None,
        subtest_id: Optional[str] = None,
        description: Optional[str] = None,
        result: Any = None,
        comments: Optional[str] = None,
        details: Optional[str] = None,
        color: Optional[str] = None,
    ):
        """
        Subtest that compares two other subtests.

        :param subtest1: First subtest to compare
        :type subtest1: SubTest
        :param subtest2: Second subtest to compare
        :type subtest2: SubTest
        :param passed: Condition for passed test Must be eq, gt or lt,
                       defaults to 'eq'
        :type passed: str, optional
        :param conditional: Condition for a conditional pass, defaults to None
        :type conditional: Optional[str], optional
        :param subtest_id: ID of this subtest, defaults to None
        :type subtest_id: Optional[str], optional
        :param description: Description of the test, defaults to None
        :type description: Optional[str], optional
        :param result: Test result, defaults to None
        :type result: Any, optional
        :param comments: Comments to add in report, defaults to None
        :type comments: Optional[str], optional
        :param details: Details about result (usually link to other page),
                        defaults to None
        :type details: Optional[str], optional
        :param color: Color of the test (Lime, Yellow, Red), defaults to None
        :type color: Optional[str], optional
        """

        if passed not in COMPARISONS and passed is not None:
            raise ValueError(
                f"Supported passed values are {list(COMPARISONS.keys())}."
            )
        if conditional not in COMPARISONS and conditional is not None:
            raise ValueError(
                f"Supported conditional values are {list(COMPARISONS.keys())}."
            )

        if description is None:
            description = f"Comparison: {subtest1.id} {passed} {subtest2.id} ?"

        super().__init__(
            subtest_id=subtest_id,
            description=description,
            result=result,
            comments=comments,
            details=details,
            color=color,
        )

        self.passed = passed
        self.conditional = conditional
        self.subtest1 = subtest1
        self.subtest2 = subtest2

    def check_types(self):
        result_types = (numbers.Number, pd.Series)
        for st in [self.subtest1, self.subtest2]:
            if st.result is None:
                raise TypeError(
                    "Test results should not be None. Make sure to run tests "
                    "before comparing them (this can be done by reordering "
                    "the subtest list)."
                )
            elif not isinstance(st.result, result_types):
                msg = (
                    "Comparison tests work only with the following result"
                    f"types: {result_types}"
                )
                raise TypeError(msg)

    def prepare_series(self):
        r1, r2 = self.subtest1.result, self.subtest2.result

        if isinstance(r1, pd.Series) and isinstance(r2, pd.Series):
            dtype = r1.dtype
            r1, r2 = r1.align(r2, fill_value=0)
            r1 = r1.astype(dtype)
            r2 = r2.astype(dtype)

        return r1, r2

    def run(self):

        self.check_types()

        r1, r2 = self.prepare_series()

        # Check if main condition passed
        check_passed = COMPARISONS[self.passed](r1, r2)

        self.result = check_passed

        # If conditinal passed is defined (yellow output)
        if self.conditional is not None:
            check_conditional = COMPARISONS[self.conditional](r1, r2)
        else:
            check_conditional = False

        # Set color based on status
        if np.all(self.result):
            self.color = "Lime"
        elif np.all(check_conditional):
            self.color = "Yellow"
        else:
            self.color = "Red"


class CountOutTest(SubTest):
    def __init__(
        self,
        ind_df: DataFrame,
        output_path: str,
        output_hkeys: List[str],
        unique: bool = False,
        description: Optional[str] = None,
        group_kwds: Optional[Sequence[str]] = None,
    ):
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
                           The default ["KW_OUTPUT", "KW_FIBER"], unless KW_FIBER is all NA
                           is used if None.
        :type group_kwds: Optional[Sequence[str]], optional
        """

        self.unique = unique
        self.output_path = output_path

        if description is None:
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
            # ???: Should we use drs params to get keys here or is literal ok?
            # NOTE: We could remove KW_DPRTYPE in future if not needed
            # self.group_kwds = ["KW_OUTPUT", "KW_DPRTYPE", "KW_FIBER"]
            self.group_kwds = ["KW_OUTPUT"]
            if not self.ind_df["KW_FIBER"].isna().all():
                self.group_kwds.append("KW_FIBER")

    def run(self):
        if not self.ind_df.empty:
            if self.unique:
                self.result = self.ind_df.groupby(
                    self.group_kwds
                ).FILENAME.nunique()
            else:
                self.result = self.ind_df.groupby(
                    self.group_kwds
                ).FILENAME.count()
        else:
            self.result = pd.Series(
                0,
                index=pd.Index(self.output_hkeys, name="KW_OUTPUT"),
                name="FILENAME",
            )


class CountQCTest(SubTest):
    def __init__(self, log_df: DataFrame, test_html_path: str):
        """
        Count files that failed one or more QC

        :param log_df: Log dataframe
        :type log_df: DataFrame
        :param test_html_path: HTML report path of parent test
        :type test_html_path: str
        """
        super().__init__(
            description="# of log entries that failed one or more QC"
        )

        self.log_df = log_df
        self.test_html_path = test_html_path

    def run(self):

        log_qc_failed = self.log_df[~self.log_df.PASSED_ALL_QC.astype(bool)]
        self.result = len(log_qc_failed)

        odometer_flag = "ODOMETER" in log_qc_failed.columns

        if self.result > 0 and odometer_flag:
            self.comments = "One or more recipe have failed QC."
            self
            log_reset = log_qc_failed.reset_index()
            data_dict_check_qc = {
                "Night": log_reset.OBS_DIR.values,
                "Odometer": log_qc_failed.ODOMETER.values,
                "QC_STRING": log_qc_failed.QC_STRING.values,
            }
            self.details = inspect_table(
                self.test_html_path,
                self.id,
                data_dict_check_qc,
                "Odometers that Failed One or More Quality Control",
            )

        if self.result > 0 and not odometer_flag:
            self.comments = "One or more recipe have failed QC."
            log_reset = log_qc_failed.reset_index()
            data_dict_check_qc = {
                "Night": log_reset.OBS_DIR.values,
                "QC_STRING": log_qc_failed.QC_STRING.values,
            }
            self.details = inspect_table(
                self.test_html_path,
                self.id,
                data_dict_check_qc,
                "Nights that Failed Quality Control",
            )


class PlotQCTest(SubTest):
    def __init__(
        self, log_df: DataFrame, test_html_path: str, recipe_name: str
    ):
        """
        Plot QC results for test

        :param log_df: Log dataframe
        :type log_df: DataFrame
        :param test_html_path: Path to HTML report of parent test
        :type test_html_path: str
        :param recipe_name: Name of the recipe under test
        :type recipe_name: str
        """
        super().__init__(
            description=(
                "Plot of the various QCs as a function of odometer,"
                " order or time"
            )
        )

        self.log_df = log_df
        self.test_html_path = test_html_path
        self.recipe_name = recipe_name

    def run(self):
        qc_names = self.log_df.QC_NAMES.str.split(r"\|\|", expand=True)
        qc_names = qc_names.iloc[0]  # Keep only one row
        qc_values = self.log_df.QC_VALUES.str.split(r"\|\|", expand=True)
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

        odometer_flag = "ODOMETER" in self.log_df.columns
        order_flag = len(qc_names) == 49

        if odometer_flag:
            data_dict_qc_plot = {"Night": qc_values.index.tolist()}
            data_dict_qc_plot["Odometer"] = self.log_df.ODOMETER.tolist()
            for key, series in qc_values.items():
                data_dict_qc_plot[key] = series.tolist()

        elif order_flag:
            data_dict_qc_plot = {
                "Order": list(range(1, 50)),
                qc_names[0]: qc_values.values[0],
            }

        else:
            data_dict_qc_plot = {"Night": qc_values.index.tolist()}
            for key, series in qc_values.items():
                # HACK: Just making the code run, probably a a bad result and need to fix underlying cause in tests setup
                if key is None:
                    key = "None"
                data_dict_qc_plot[key] = series.tolist()

        if "Night" in data_dict_qc_plot and "other" in data_dict_qc_plot["Night"]:
            other_list = [n for n in data_dict_qc_plot["Night"] if n == "other"]
            data_dict_qc_plot["Night"] = [n for n in data_dict_qc_plot["Night"] if n != "other"]
            self.comments = f"Dropped {len(other_list)} rows with OBS_DIR = other."

        self.result = "See the Details column"
        self.details = inspect_plot(
            self.test_html_path,
            self.id,
            data_dict_qc_plot,
            f"{self.recipe_name}.py Quality Control",
        )


class CountEndedTest(SubTest):
    def __init__(self, log_df: DataFrame, test_html_path: str):
        """
        Count test that failed to finish.

        :param log_df: Log dataframe
        :type log_df: DataFrame
        :param test_html_path: Path to HTML report of parent test
        :type test_html_path: str
        """
        super().__init__(description="# of log entries that failed to finish")

        self.log_df = log_df
        self.test_html_path = test_html_path

    def run(self):

        log_ended_false = self.log_df[~self.log_df.ENDED.astype(bool)]
        self.result = len(log_ended_false)

        if self.result > 0:
            self.comments = "One or more recipe have failed to finish."
            self
            log_reset = log_ended_false.reset_index()
            data_dict_check_ended = {
                "Night": log_reset.OBS_DIR.values,
                "ERRORS": log_ended_false.ERRORS.values,
                "LOGFILE": log_ended_false.LOGFILE.values,
            }
            self.details = inspect_table(
                self.test_html_path,
                self.id,
                data_dict_check_ended,
                "Nights that Failed to finish",
            )


class CountIndexCalib(SubTest):
    def __init__(
        self,
        ind_df: DataFrame,
        calib_key: str = "CALIB_KEY",
        qc_key: str = "QC_PASSED",
    ):
        """
        Count calibration files in index.

        :param ind_df: Index dataframe
        :type ind_df: DataFrame
        :param calib_key: Index key that has calibdb keys in index df,
                          defaults to "CALIB_KEY"
        :type calib_key: str, optional
        :param qc_key: Index df column saying if QC has passed for this file,
                       defaults to "QC_PASSED"
        :type qc_key: str, optional
        """
        super().__init__(
            description="# of index entries for each calibdb key",
        )
        # Keep only index entries that passed QC
        self.ind_df = ind_df[ind_df[qc_key]]
        self.key = calib_key

    def run(self):
        self.result = self.ind_df.groupby(self.key).FILENAME.nunique()
        if len(self.result) == 0:
            self.result = 0


class CheckIndexCalibFiles(SubTest):
    def __init__(
        self,
        ind_df: DataFrame,
        calib_list: List[str],
        test_html_path: str,
        calib_key: str = "CALIB_KEY",
        qc_key: str = "QC_PASSED",
    ):
        """
        Check that calib files in index appear in calibdb directory

        :param ind_df: Index dataframe
        :type ind_df: DataFrame
        :param calib_list: List of files in calibdb dir
        :type calib_list: List[str]
        :param calib_key: Index key that has calibdb keys in index df,
                          defaults to "CALIB_KEY"
        :type calib_key: str, optional
        :param qc_key: Index df column saying if QC has passed for this file,
                       defaults to "QC_PASSED"
        :type qc_key: str, optional
        """

        super().__init__(
            description="Check that calib files in index are in calibdb dir"
        )

        # NOTE: For 0.7, there will probably be a way to do more detailed, key-by-key
        # comparison, but in 0.6 this would require reading many^TM fits files to get
        # all keys

        self.ind_df = ind_df[ind_df[qc_key]].dropna(subset=[calib_key])
        self.calib_list = calib_list
        self.test_html_path = test_html_path

    def run(self):

        in_calib_mask = self.ind_df.FILENAME.isin(self.calib_list)

        self.result = np.all(in_calib_mask)

        self.details = inspect_table(
            self.test_html_path,
            self.id,
            self.ind_df[~in_calib_mask].to_dict(orient="list"),
            "Missing files in calibdb directory",
        )

        # Set color based on status
        if self.result:
            self.color = "Lime"
        else:
            self.color = "Red"


class CheckCalibEntriesFiles(SubTest):
    def __init__(
        self,
        calib_df: DataFrame,
        calib_list: List[str],
        test_html_path: str,
    ):
        """
        Check that calibdb entries are in calibdb dir

        :param calib_df: Dataframe with calibdb info
        :type calib_df: DataFrame
        :param calib_list: List of files in calibdb
        :type calib_list: List[str]
        """

        super().__init__(
            description="Check that calibdb entries are in calibdb dir"
        )

        # FUTURE: 0.7 will probably have way to do more detailed, key-by-key
        # comparison, but in 0.6 this would require reading many^TM fits files
        # to get all keys

        self.calib_df = calib_df
        self.calib_list = calib_list
        self.test_html_path = test_html_path

    def run(self):

        in_dir_mask = self.calib_df.FILENAME.isin(self.calib_list)

        self.result = np.all(in_dir_mask)

        self.details = inspect_table(
            self.test_html_path,
            self.id,
            self.calib_df[~in_dir_mask].to_dict(orient="list"),
            "Missing files in calibdb directory",
        )

        # Set color based on status
        if self.result:
            self.color = "Lime"
        else:
            self.color = "Red"


class CountCalibEntries(SubTest):
    def __init__(
        self, calib_df: DataFrame, calib_keys: List[str], ref: bool = False
    ):
        """
        Count entries in calibdb

        :param calib_df: Dataframe with calibdb info
        :type calib_df: DataFrame
        :param calib_keys: Calib db keys for current test
        :type calib_keys: List[str]
        :param ref: True if ref recipe, defaults to False
        :type ref: bool, optional
        """

        super().__init__(description="# of entries in calib DB")

        self.calib_df = calib_df
        self.calib_keys = list(set(calib_keys))
        self.isref = ref

    def run(self):
        zero_count = pd.Series(0, index=self.calib_keys)

        if not self.calib_df.empty:
            # TODO: Check with Neil that OK to use nunique here
            # calib_count = self.calib_df.groupby("KEYNAME").size()
            calib_count = self.calib_df.groupby("KEYNAME").nunique().FILENAME
            if not self.isref:
                ref_mask = self.calib_df["REFCAL"].astype(bool)
                if ref_mask.sum() == 0:
                    ref_calib_count = zero_count.copy()
                else:
                    ref_calib_df = self.calib_df[ref_mask]
                    ref_calib_count = ref_calib_df.groupby(
                        "KEYNAME"
                    ).size()
                    msg = "Additional entries with ref == 1"
                    self.comments = f"{msg}: {ref_calib_count}"

                calib_count -= ref_calib_count

        else:
            self.result = 0
            return

        self.result = calib_count


class CheckUsedCalibs(SubTest):
    def __init__(
        self,
        used_calib_df: DataFrame,
        test_html_path: str,
        calib_df: DataFrame,
        delta_mjd_max: float = 7.0,
    ):
        """
        Check if used calibrations were close in time.

        :param used_calib_df: Dataframe of used calibration files and time
                              diffs
        :type used_calib_df: DataFrame
        :param test_html_path: Path of HTML report of parent test.
        :type test_html_path: str
        :param calib_df: Dataframe with calibdb info *for all recipes*
        :type calib_df: DataFrame
        :param delta_mjd_max: Flag above this time difference (in days),
                              defaults to 7.0
        :type delta_mjd_max: float, optional
        """
        super().__init__(
            description=f"# of sci files with |delta MJD| from calib > {delta_mjd_max}"
        )

        # Drop calib files not relevant to sci file under test with dropna
        mask = ~used_calib_df.isna().all(axis=0).groupby(level=0).any()
        self.used_calib_df = used_calib_df[mask.index[mask]]
        self.calib_df = calib_df
        self.test_html_path = test_html_path
        self.dmax = 7.0

    def run(self):

        if len(self.used_calib_df.columns) > 0:
            calib_with_fname = self.calib_df.set_index("FILENAME")
            fdf = self.used_calib_df.xs(
                "CALIB_FILE", axis=1, level=1
            ).reset_index(drop=True)

            stacked_fdf = fdf.stack().reset_index(0, drop=True)
            ind_for_db = stacked_fdf[stacked_fdf.isin(self.calib_df.FILENAME)]
            calib_ind_used = calib_with_fname.loc[ind_for_db]
            drop_map = (
                calib_ind_used.REFCAL.astype(int).groupby("FILENAME").all()
            )
            drop_mask = fdf.replace(drop_map)
            drop_mask_all = drop_mask.all()
            # Always keep CDBWAVE because means that wavelength solution
            # in header is old, even if it is indeed a ref==1 file
            if "CDBWAVE" in drop_mask_all.index:
                drop_mask_all["CDBWAVE"] = False
            drop_inds = drop_mask_all.index[~drop_mask_all]
            # TODO: For apero_mk_template, end up with empty used_calib_df -> is this OK
            self.used_calib_df = self.used_calib_df[drop_inds]
            if len(self.used_calib_df.columns) == 0:
                n_fail = 0
            else:
                delta_mjd_df = self.used_calib_df.xs(
                    "DELTA_MJD", axis=1, level=1
                ).abs()
                check_per_ctype = delta_mjd_df > self.dmax
                bool_mask = check_per_ctype.any(axis=1)
                n_fail = bool_mask.sum()
        else:
            n_fail = 0

        self.result = n_fail

        if n_fail == 0:
            self.color = "Lime"
            return
        else:
            self.color = "Red"

        inspect_df = self.used_calib_df[bool_mask]

        nights = inspect_df.index.get_level_values("OBS_DIR")
        if "other" in nights:
            nother = nights.value_counts()["other"]
            inspect_df = inspect_df.drop("other")
            self.comments = f"Dropped {nother} rows with OBS_DIR = other."

        self.details = delta_mjd_plot(
            self.test_html_path,
            self.id,
            inspect_df,
            "Time difference between calib files and science file (hover for file name)",
        )


class CountTelluEntries(SubTest):
    def __init__(self, tellu_df: DataFrame):
        """
        Count entries in telluDB

        :param tellu_df: Dataframe with tellu db
        :type tellu_df: DataFrame
        """

        super().__init__(description="# of entries in tellu DB")

        self.tellu_df = tellu_df

    def run(self):
        if not self.tellu_df.empty:
            tellu_count = self.tellu_df.groupby("KEYNAME").size()

        else:
            self.result = 0
            return

        self.result = tellu_count
