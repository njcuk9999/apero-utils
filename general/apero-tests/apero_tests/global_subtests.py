import warnings

import pandas as pd

from apero_tests.display import inspect_table, inspect_hkey_plot
from apero_tests.drs_test import DrsTest
from apero_tests.subtest import SubTest


class CustomBadpixTest(SubTest):
    """
    Quick example of how to write a very simple test.
    See below or in subtest.py for more complete ones.
    """

    def __init__(self, input_msg: str):
        """
        Subtest that does nothing, example for how we can define custom tests.
        It is 100% sure to pass and has a smiley so I'm leaving it here for
        now, pour le moral :)
        """
        self.msg = input_msg

        super().__init__(description="Test to show something")

    def run(self):

        self.result = self.msg + " :)"

        self.color = "Lime"


class GlobalIndexCheck(SubTest):
    def __init__(self, global_bad_index: pd.DataFrame, parent_test: DrsTest):
        """
        Check that index files have a PID matched in log, otherwise report

        :param global_bad_index: Full bad index (all files without match)
        :type global_bad_index: pd.DataFrame
        :param parent_test: Parent test calling the subtest. Needed for HTML
                            report path, tno using path directly in case it has
                            not been generated when creating the subtest.
        :type parent_test: DrsTest
        """
        super().__init__(
            description="Files in index that have no PID matched in log"
        )

        self.global_bad_index = global_bad_index
        self.parent_test = parent_test

        self.group_kwds = ["PID_TYPE", "IN_INDEX", "IN_LOG"]
        self.inspect_kwds = ["KW_PID", "KW_OUTPUT", "KW_DPRTYPE", "KW_FIBER"]
        self.count_column = "FILENAME"

    def run(self):

        if len(self.global_bad_index) == 0:
            self.result = 0
            self.color = "Lime"
            return

        try:
            global_bad_index_summary = self.global_bad_index.groupby(
                self.group_kwds, dropna=False
            )[self.count_column].count()
        except TypeError:
            pd_msg = (
                "Your pandas version does not support NaN grouping. "
                "Some entries might be missing from the index summary"
            )
            warnings.warn(pd_msg, RuntimeWarning)
            global_bad_index_summary = self.global_bad_index.groupby(
                self.group_kwds
            )[self.count_column].count()

        bad_index_reset = self.global_bad_index.reset_index()
        data_dict_index_pid = {
            "Night": bad_index_reset.OBS_DIR,
            "File": bad_index_reset.FILENAME,
        }
        for k in self.group_kwds + self.inspect_kwds:
            data_dict_index_pid[k] = bad_index_reset[k]

        self.result = global_bad_index_summary
        self.color = "Red"
        self.comments = (
            "Some files are in the index but have no PID match in log"
        )
        self.details = inspect_table(
            self.parent_test.html_path,
            self.id,
            data_dict_index_pid,
            title="Files that are in index but have no PID match in log",
        )


class CheckMissingAdded(SubTest):
    def __init__(self, missing_ind_df: pd.DataFrame, parent_test: DrsTest):
        """
        Check missing files that were added to the index

        :param missing_ind_df: Dataframe of missing index file
        :type missing_ind_df: pd.DataFrame
        :param parent_test: Parent test calling the subtest.
        :type parent_test: DrsTest
        """
        super().__init__(description="Files on disk that were not in index")

        self.inspect_kwds = [
            "KW_OUTPUT",
            "KW_DPRTYPE",
            "KW_FIBER",
            "KW_PID",
        ]
        self.missing_ind_df = missing_ind_df
        self.parent_test = parent_test

    def run(self):
        n_missing_added = len(self.missing_ind_df)

        self.result = n_missing_added

        if n_missing_added > 0:
            df_reset = self.missing_ind_df.reset_index()
            data_dict = {
                "Night": df_reset.OBS_DIR,
                "File": df_reset.FILENAME,
            }
            for k in self.inspect_kwds:
                data_dict[k] = data_dict[k]

            self.details = inspect_table(
                self.parent_test.html_path,
                self.id,
                data_dict,
                title="Files that were on disk and not in index (added to index df)",
            )
            self.color = "Yello"
        else:
            self.color = "Lime"


class CheckNotFound(SubTest):
    def __init__(self, not_found: pd.Series, parent_test: DrsTest):
        """
        Check if some files were not found on disk. For example, a copied
        calib files whose original file is missing

        :param not_found: Series of missing files
        :type not_found: pd.Series
        :param parent_test: DRS Test that calls the subtest.
        :type parent_test: DrsTest
        """
        super().__init__(
            description="Copied files with original file not found"
        )

        self.not_found = not_found
        self.parent_test = parent_test

    def run(self):
        n_not_found = len(self.not_found)

        self.result = n_not_found

        if n_not_found > 0:
            data_dict = {"File": self.not_found.values}

            self.details = inspect_table(
                self.parent_test.html_path,
                self.id,
                data_dict,
                title="Copied files with original file not found",
            )
        else:
            self.color = "Lime"


class CustomWaveTest(SubTest):
    def __init__(
        self, hkey_df: pd.DataFrame, parent_test: DrsTest):
        """
        Plot the FP cavity properties as a function of time

        :param hkey_df: Log dataframe
        :type hkey_df: DataFrame
        :param parent_test: Parent test calling the subtest. Needed for HTML
                            report path, not using path directly in case it has
                            not been generated when creating the subtest.
        :type parent_test: DrsTest
        """
        super().__init__(
            description=(
                "Plot the FP cavity properties as a function of time"
            )
        )

        self.hkey_df = hkey_df
        self.parent_test = parent_test

    def run(self):
        
        # We only want the WCAV* keys.
        keys = [kw for kw in self.hkey_df.columns if kw.startswith("WCAV")]

        self.result = "See the Details column"
        self.details = inspect_hkey_plot(
            self.parent_test.html_path,
            self.id,
            self.hkey_df[keys].drop('other'),
            f"{self.parent_test.recipe_name}.py Quality Control",
        )


class CustomFluxTest(SubTest):
    def __init__(
        self, hkey_df: pd.DataFrame, parent_test: DrsTest):
        """
        Test that the flux in AB is roughly equal to flux in A + B

        :param hkey_df: Log dataframe
        :type hkey_df: DataFrame
        :param parent_test: Parent test calling the subtest. Needed for HTML
                            report path, not using path directly in case it has
                            not been generated when creating the subtest.
        :type parent_test: DrsTest
        """
        super().__init__(description="Check that Flux AB = Flux A + Flux B")

        self.hkey_df = hkey_df
        self.parent_test = parent_test

    def run(self):

        # We only want the EXTSN* keys.
        keys = [kw for kw in self.hkey_df.columns if kw.startswith("EXTSN")]
        self.hkey_df = self.hkey_df[keys].drop('other')

        # Find 'good' file type (objects, e2ds, e2dsff, tcorr)
        file_list = self.hkey_df.index.get_level_values(1).tolist()
        good = [n for n, x in enumerate(file_list) if ('o_pp' in x) and 
            ('AB.fits' in x or 'A.fits' in x or 'B.fits' in x) and 
            ('e2ds' in x or 'tcorr' in x)]
        self.hkey_df = self.hkey_df.iloc[good]

        self.result = "See the Details column"
        self.details = inspect_hkey_plot(
            self.parent_test.html_path,
            self.id,
            self.hkey_df,
            f"{self.parent_test.recipe_name}.py Quality Control",
        )
