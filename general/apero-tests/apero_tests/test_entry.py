from typing import Optional, Sequence

from pandas.core.frame import DataFrame


class TestEntry:
    def __init__(self):

        self.id = None
        self.description = "Generic Drs Test Entry"
        self.result = None
        self.comments = None
        self.details = None

    @property
    def id(self):
        pass

    @property
    def description(self):
        pass

    @property
    def result(self):
        pass

    @property
    def comments(self):
        pass

    @property
    def details(self):
        pass


class CountLogEntry(TestEntry):
    def __init__(self,
                 log_df: DataFrame,
                 master_flag: str = "--master",
                 group_kwds: Optional[Sequence[str]] = None):

        self.log_df = log_df
        self.master_flag = master_flag
        if group_kwds is not None:
            self.group_kwds = group_kwds
        else:
            self.group_kwds = ["RECIPE", "SUBLEVEL", "LEVEL_CRIT"]

    def run(self):

        master_mask = self.log_df.RUNSTRING.str.contains("--master")
        master_df = self.log_df[master_mask]
        master_count = master_df.groupby(self.group_kwds).count().KIND
        tot_count = self.log_df.groupby(self.group_kwds).count().KIND

        log_count = tot_count.sub(master_count, fill_value=0).astype(int)

        self.result = log_count
