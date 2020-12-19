"""
Factory function for test
"""

from apero_tests import Test


from preprocessing_test1 import PPTest
from darkmaster_test1 import DarkMTest
from badpixel_test1 import BadPixTest
from localisation_test1 import LocTest


TESTDICT = {
        # Preprocessing
        'preprocessing_test1': PPTest,

        # Calibration
        'darkmaster_test1': DarkMTest,
        'badpixel_test1': BadPixTest,
        'localisation_test1': LocTest,
        'shapemaster_test1': None,
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

    return test()
