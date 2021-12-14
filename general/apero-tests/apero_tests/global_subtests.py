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

        super().__init__(description=f"Test to show something")

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
        for k in self.group_kwds:
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
