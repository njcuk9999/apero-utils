"""
Factory function for test

@author: vandal
"""

from .tests import Test
from .preprocessing import PPTest
from .darkmaster import DarkMTest
from .badpixel import BadPixTest
from .localisation import LocTest
from .shapemaster import ShapeMTest
from .shape import ShapeTest
from .flat import FlatTest
from .thermal import ThermalTest
from .leakmaster import LeakMTest


# Dictionary of all tests with id and definition
TESTDICT = {
        # Preprocessing
        'preprocessing_test1': PPTest,

        # Calibration
        'darkmaster_test1': DarkMTest,
        'badpixel_test1': BadPixTest,
        'localisation_test1': LocTest,
        'shapemaster_test1': ShapeMTest,
        'shape_test1': ShapeTest,
        'flat_test1': FlatTest,
        'thermal_test1': ThermalTest,
        'leakmaster_test1': LeakMTest,
        'masterwavelength_test1': None,
        'wavelength_test1': None,
        'leak_test1': None,
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
