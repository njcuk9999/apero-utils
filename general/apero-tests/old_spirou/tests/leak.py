"""
Check if leak calib worked fine.

Tests preformed
check1: how many recipes were run (cal_leak_{instrument} in log.fits)?
        how many in the master directory?
check2: how many of each output do we have?
        output1: {ODOMETER_CODE}_pp_e2ds_{FIBER}.fits
        output2: {ODOMETER_CODE}_pp_e2dsff_{FIBER}.fits
        output3: {ODOMETER_CODE}_pp_s1d_w_{FIBER}.fits
        output4: {ODOMETER_CODE}_pp_s1d_v_{FIBER}.fits 
check3: how many of each unique output do we have?
        output1: {ODOMETER_CODE}_pp_e2ds_{FIBER}.fits
        output2: {ODOMETER_CODE}_pp_e2dsff_{FIBER}.fits
        output3: {ODOMETER_CODE}_pp_s1d_w_{FIBER}.fits
        output4: {ODOMETER_CODE}_pp_s1d_v_{FIBER}.fits 
stop1:  check3 == check1?
check4: using the log.fits how many files failed to finish? Which odometers?
        Which nights? Why (using the ERRORS and LOGFILE columns)?

@author: charles
"""
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .tests import CalibTest


class LeakTest(CalibTest):
    """LeakTest."""
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
        return 'Leak Correction Recipe #1'

    @property
    def test_id(self) -> str:
        """test_id.

        :rtype: str
        """
        return 'leak_test1'

    @property
    def output_list(self) -> List[str]:
        """List of output string patterns

        :return: output_list
        :rtype: list[str]
        """

        return ['*_pp_e2dsff_{FIBER}_wave_night_{FIBER}.fits',
                '*_pp_e2dsff_{FIBER}_wave_hclines_{FIBER}.fits',
                '*_pp_e2dsff_{FIBER}_wave_fplines_{FIBER}.fits',
                '*_pp_e2dsff_{FIBER}ccf_{FIBER}.fits']

    @property
    def calibdb_list(self) -> List[str]:
        """List of calibDB entries

        :return: calibdb_list
        :rtype: list[str]
        """
        return ['']

    @property
    def previous_calibs(self) -> List[str]:
        """previous_calibs.

        :rtype: List[str]
        """
        return [
            'CDBBAD', 'CDBBACK', 'CDBORDP', 'CDBLOCO', 'CDBSHAPX', 'CDBSHAPY',
            'CDBSHAPL', 'CDBFLAT', 'CDBBLAZE', 'CDBTHERM'
        ]

    @property
    def recipe(self) -> List[str]:
        """Recipe name

        :return: output_list
        :rtype: list[str]
        """
        return 'cal_wave_night_{}'.format(self.instrument.lower())

    @property
    def ismaster(self) -> bool:
        """Is the test for a master recipe.

        :rtype: bool
        """
        return False

    @property
    def fibers(self) -> List[str]:
        """fibers.

        :rtype: List[str]
        """
        return ['AB', 'A', 'B', 'C']

    @property
    def calls_extract(self) -> bool:
        """Does this method call extract.

        :rtype: bool
        """
        return False

    # =========================================================================
    # Run the full test
    # =========================================================================

    def runtest(self):
        """runtest."""

        comments_check1, dup = self.check_duplicates()

        comments_check4, inspect_check4 = self.check_ended(ncheck=4)

        dict_stop1 = self.stop_output_log(dup, nstop=1)


        html_dict = {
            # Summary header info
            'name': self.name,
            'setup': self.setup,
            'instrument': self.instrument,
            'recipe': self.recipe,
            'date': self.date,
            'reduced_path': self.reduced_path,
            'output_list': self.output_list,
            'calibdb_list': self.calibdb_list,
            'calibdb_path': self.calibdb_path,

            # check 1 for logs
            'recipe_num_logfits': self.log_recipe.tot_num,

            # check 2 for outputs
            'output_num_total': self.output_num_total,

            # check 3
            'output_num_unique': self.output_num_unique,

            # stop 1: output==log
            'dict_stop1': dict_stop1,

            # check 4: not ended
            'log_ended_false': self.log_all.ended_false,
            'comments_check4': comments_check4,
            'inspect_check4': inspect_check4,

        }

        self.gen_html(html_dict)


if __name__ == '__main__':
    test = LeakTest()
    test.runtest()
