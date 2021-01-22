"""
Check if dark master calib worked fine.

Tests performed:
    check1: how many recipes were run? (cal_dark_master_{instrument} in
            other/log.fits)
    check2: how many of each output do we have?
            output1: {ODOMETER_CODE}_pp_dark_master.fits
    check3: how many of each unique output do we have?
            output1: {ODOMETER_CODE}_pp_dark_master.fits
    stop1: check3 == check1?
    check4: using the log.fits how many entry failed one or more QC?
            Which nights? Which QC?
    check5: using the log.fits how many entry failed to finish? Which nights?
            Why (using the ERRORS and LOGFILE columns)?
    check6: how many entry DARKM in master_calib_{INSTRUMENT}.txt?
    check7: for each calib entry how many are in the calibDB?
    stop2: check7 == check6?

@author: charles
"""
from typing import Optional, Union, List

import pandas as pd

from .tests import Test, CalibTest
from . import utils as ut


class DarkMTest(CalibTest):
    """DarkMTest."""

    def __init__(self,
                 inst: str = 'SPIROU',
                 setup: Optional[str] = None,
                 logdir: Union[str, list] = 'other'):
        """__init__.

        :param inst:
        :type inst: str
        :param setup:
        :type setup: Optional[str]
        :param logdir:
        :type logdir: Union[str, list]
        """
        super().__init__(inst=inst, setup=setup, logdir=logdir)

    # =========================================================================
    # Specific properties
    # =========================================================================
    @property
    def name(self) -> str:
        """name.

        :rtype: str
        """
        return 'Dark Master Recipe Test #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        return 'darkmaster_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return ['*_pp_dark_master.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['DARKM']

    @property
    def previous_calibs(self) -> None:
        """previous_calibs.

        :rtype: None
        """

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return f'cal_dark_master_{self.instrument.lower()}'

    @property
    def ismaster(self) -> bool:
        """Is the test for a master recipe.

        :rtype: bool
        """
        return True

    @property
    def fibers(self) -> None:
        """fibers.
        No fibers for darkmaster
        """

    # =========================================================================
    # Overwritten parent methods
    # =========================================================================
    @property
    def output_num_entry(self) -> pd.Series:
        """Total number of entries per output in calibdb.

        For a master recipe, we don't need to subtract master entries.

        :rtype: pd.Series
        """
        return self.tot_num_entry

    @property
    def output_missing(self) -> pd.DataFrame:
        """Overwrite output_missing for specific case of dark master.

        Missing outputs defined as every night without a dark master calib

        :rtype: pd.DataFrame
        """
        # NOTE: should change in APERO v0.7
        output_frame = self.output_files.index.to_frame().reset_index(drop=True)
        night_comb = [(p, n)
                      for p in self.output_list
                      for n in self.reduced_nights]
        all_nights = pd.DataFrame(night_comb, columns=['PATTERN', 'DIRECTORY'])
        output_missing = all_nights[~all_nights.isin(output_frame).all(axis=1)]

        return output_missing

    def check_duplicates(self) -> pd.DataFrame:
        """check_duplicates.

        Duplicate defined as every night with more than one dark master calib

        :rtype: pd.DataFrame
        """
        # NOTE: should change in APERO v0.7
        dup = self.output_num_night[self.output_num_night > 1]
        dup.name = 'COUNT'
        dup = dup.reset_index()

        return dup

    # =========================================================================
    # Running the test
    # =========================================================================
    def runtest(self):
        """runtest."""

        # =====================================================================
        # Inspect outputs
        # =====================================================================
        # Checking for duplicates
        dup = self.check_duplicates()

        # QC/ENDED
        comments_check4, inspect_check4 = self.check_qc(ncheck=4)
        comments_check5, inspect_check5 = self.check_ended(ncheck=5)

        dict_stop1 = self.stop_output_log(dup)

        # =====================================================================
        # Inspect calibDB
        # =====================================================================
        # Check if some files have duplicaltes in db (empty if none)
        # Use tot_num_entry for dark MASTER recipe (i.e. not remove masters)
        calib_dup_mask = self.master_calib_df.duplicated(keep=False)
        # master_mask = self.master_calib_df.MASTER  # Master recip
        calib_dup = self.master_calib_df[calib_dup_mask]
        calib_dup = calib_dup.set_index('FILE', append=True)
        ind_names = calib_dup.index.names  # get file and key to count
        calib_dup['COUNT'] = calib_dup.groupby(ind_names).size()
        calib_dup = calib_dup.reset_index('FILE')
        calib_dup = calib_dup.drop_duplicates()

        dict_stop2 = self.stop_calibdb(calib_dup)

        html_dict = {
                # Summary header info
                'name': self.name,
                'setup': self.setup,
                'instrument': self.instrument,
                'recipe': self.recipe,
                'reduced_path': self.reduced_path,
                'date': self.date,
                'output_list': self.output_list,
                'calibdb_list': self.calibdb_list,
                'calibdb_path': self.calibdb_path,

                # Check 1: number of calls in logfits
                'recipe_num_logfits': self.log_recipe.tot_num,  # Master recipe

                # Check 2 number of outputs
                'output_num_total': self.output_num_total,

                # Check 3 number of unique outputs
                'output_num_unique': self.output_num_unique,

                # Stop 1
                'dict_stop1': dict_stop1,

                # Check 4: QC failed
                'log_qc_failed': self.log_all.qc_failed,
                'comments_check4': comments_check4,
                'inspect_check4': inspect_check4,

                # Check 5: not ended
                'log_ended_false': self.log_all.ended_false,
                'comments_check5': comments_check5,
                'inspect_check5': inspect_check5,

                # Check 6: calibdb entries
                'output_num_entry': self.tot_num_entry,  # Master recipe
                # 'comments_check6': comments_check6,

                # Check 7: calibdb outputs
                'output_num_calibdb': self.output_num_calibdb,
                'output_dict': self.calib_output_dict,

                # Stop 2: calib output == entries
                'dict_stop2': dict_stop2,
                }

        self.gen_html(html_dict)


if __name__ == '__main__':
    test = DarkMTest()
    test.runtest()
