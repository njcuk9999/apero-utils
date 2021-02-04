"""
Check if bad pixel calib worked fine.

Tests preformed
    check1: how many recipes were run (cal_badpix_{instrument} in log.fits)?
            how many in the master directory?
    check2: how many of each output do we have?
            #output1: {ODOMETER_CODE}_pp_badpixel.fits
            #output2: {ODOMETER_CODE}_pp_bmap.fits
    check3: how many of each unique output do we have?
            #output1: {ODOMETER_CODE}_pp_badpixel.fits
            #output2: {ODOMETER_CODE}_pp_bmap.fits
    stop1: check3 == check1?
    check4: using the log.fits how many entry failed one or more QC?
            Which nights? Which QC?
    check5: using the log.fits how many entry failed to finish? Which nights?
            Why (using the ERRORS and LOGFILE columns)?
    check6: how many entry BADPIX and BKGRDMAP in
            master_calib_{INSTRUMENT}.txt?
    check7: for each calib entry how many are in the calibDB?
    stop2: check7 == check6?

@author: charles
"""
from typing import List, Optional, Union

from .tests import CalibTest


class BadPixTest(CalibTest):
    """BadPixTest."""

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
        return 'Bad Pixel Correction Recipe Test #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        # TODO: Use IDs derived from recipe info (name, master)
        return 'badpixel_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """
        return ['*_pp_badpixel.fits', '*_pp_bmap.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['BADPIX', 'BKGRDMAP']

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
        return f'cal_badpix_{self.instrument.lower()}'

    @property
    def ismaster(self) -> bool:
        """Is the test for a master recipe.

        :rtype: bool
        """
        return False

    @property
    def fibers(self) -> None:
        """fibers.
        No fibers for badpix
        """

    @property
    def calls_extract(self) -> bool:
        """Does this method call extract.

        :rtype: bool
        """
        return False

    # =========================================================================
    # Running the test
    # =========================================================================
    def runtest(self):
        """runtest."""

        # =====================================================================
        # Inspect outputs
        # =====================================================================
        # Checking duplicates
        comments_check1, true_dup = self.check_duplicates()

        # QC/ENDED
        comments_check4, inspect_check4 = self.check_qc(ncheck=4)
        comments_check5, inspect_check5 = self.check_ended(ncheck=5)

        dict_stop1 = self.stop_output_log(true_dup, nstop=1)

        # =====================================================================
        # Inspect calibDB
        # =====================================================================
        # Generate comment on extra master entries
        comments_check6 = ('An additional {0} {1} and {2} {3} with master = 1 '
                           'are in the master_calibDB_{4}.'.format(
                               self.master_num_entry[self.calibdb_list[0]],
                               self.calibdb_list[0],
                               self.master_num_entry[self.calibdb_list[1]],
                               self.calibdb_list[1],
                               self.instrument)
                           )


        # Get all duplicates
        calib_dup_mask = self.master_calib_df.duplicated(keep=False)
        master_mask = self.master_calib_df.MASTER
        calib_dup = self.master_calib_df[calib_dup_mask & ~master_mask]
        calib_dup = calib_dup.set_index('FILE', append=True)
        ind_names = calib_dup.index.names  # get file and key to count
        calib_dup['COUNT'] = calib_dup.groupby(ind_names).size()
        calib_dup = calib_dup.reset_index('FILE')
        calib_dup = calib_dup.drop_duplicates()


        dict_stop2 = self.stop_calibdb(calib_dup, nstop=2)

        html_dict = {
                # Summary header info
                'name': self.name,
                'setup': self.setup,
                'instrument': self.instrument,
                'date': self.date,
                'reduced_path': self.reduced_path,
                'recipe': self.recipe,
                'output_list': self.output_list,
                'calibdb_list': self.calibdb_list,
                'calibdb_path': self.calibdb_path,

                # check 1 for logs
                'recipe_num_logfits': self.log_recipe.num,
                'comments_check1': comments_check1,

                # check 2 for outputs
                'output_num_total': self.output_num_total,

                # check 3
                'output_num_unique': self.output_num_unique,

                # stop 1: output==log
                'dict_stop1': dict_stop1,

                # check 4: QC failed
                'log_qc_failed': self.log_all.qc_failed,
                'comments_check4': comments_check4,
                'inspect_check4': inspect_check4,

                # check 5: not ended
                'log_ended_false': self.log_all.ended_false,
                'comments_check5': comments_check5,
                'inspect_check5': inspect_check5,

                # Check 6: calibdb entries
                'output_num_entry': self.output_num_entry,
                'comments_check6': comments_check6,

                # Check 7: calibdb outputs
                'output_num_calibdb': self.output_num_calibdb,
                'output_dict': self.calib_output_dict,

                # Stop 2: calib output == entries
                'dict_stop2': dict_stop2,
                }

        self.gen_html(html_dict)

if __name__ == '__main__':
    test = BadPixTest()
    test.runtest()
