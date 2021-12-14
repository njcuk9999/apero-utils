from apero_tests.subtest import SubTest
import pandas as pd
import warnings

from apero_tests.drs_test import DrsTest
from apero_tests.display import inspect_table


class CustomBadpixTest(SubTest):
    """Quick example"""

    def __init__(self, input_msg: str):
        """
        Subtest that does nothing, example for how we can define custom tests
        """
        self.msg = input_msg

        super().__init__(description="Test to show something")

    def run(self):

        self.result = self.msg + " :)"

        self.color = "Lime"


# TODO: Move this to separate function or subtest (except returned things)
# Full summary of index with things to flags
class GlobalIndexCheck(SubTest):
    def __init__(self, global_bad_index: pd.DataFrame, parent_test: DrsTest):
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
            "Night": bad_index_reset.NIGHTNAME,
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
    def __init__(self, missing_ind_df, parent_test: DrsTest):
        super().__init__(
            description="Files on disk that were not in index"
        )

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
                "Night": df_reset.NIGHTNAME,
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
        else:
            self.color = "Lime"


class CheckNotFound(SubTest):
    def __init__(self, not_found: pd.Series, parent_test: DrsTest):
        super().__init__(
            description="Copied files with original file not found"
        )

        self.not_found = not_found
        self.parent_test = parent_test

    def run(self):
        n_not_found = len(self.not_found)

        self.result = n_not_found

        if n_not_found > 0:
            data_dict = {
                "File": self.not_found.values
            }

            self.details = inspect_table(
                self.parent_test.html_path,
                self.id,
                data_dict,
                title="Copied files with original file not found",
            )
        else:
            self.color = "Lime"
