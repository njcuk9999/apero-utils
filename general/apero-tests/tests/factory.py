"""
Factory function for test

@author: vandal
"""

from apero_tests import Test
from preprocessing import PPTest
from darkmaster import DarkMTest
from badpixel import BadPixTest
from localisation import LocTest
from shapemaster import ShapeMTest


TESTDICT = {
        # Preprocessing
        'preprocessing_test1': PPTest,

        # Calibration
        'darkmaster_test1': DarkMTest,
        'badpixel_test1': BadPixTest,
        'localisation_test1': LocTest,
        'shapemaster_test1': ShapeMTest,
        'shape_test1': None,
        'flat_test1': None,
        'thermal_test1': None,
        'masterleak_test1': None,
        'leak_test1': None,
        'masterwavelength_test1': None,
        'wavelength_test1': None,
        'extraction_test1': None,
        'extraction_test2': None,
        'extraction_test3': None,

        # Science
        'maketellu_test1': None,
        'fittellu_test1': None,
        'maketemplate_test1': None,
        'ccf_test1': None,
        }


def get_test(testid: str) -> Test:
    """get_test

    Function to generate any subclass of Test based on the test ID.

    :param testid: ID of the test to create
    :type testid: str

    :return: New test object
    :rtype: Test
    """

    test = TESTDICT.get(testid)

    return test() if test is not None else test
