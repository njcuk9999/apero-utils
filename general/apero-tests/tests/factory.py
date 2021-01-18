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


# Dictionary of all tests with id and definition
TESTDICT = {
        # Preprocessing
        'preprocessing': PPTest,

        # Calibration
        'darkmaster': DarkMTest,
        'badpixel': BadPixTest,
        'localisation': LocTest,
        'shapemaster': ShapeMTest,
        'shape': ShapeTest,
        'flat': FlatTest,
        'thermal': None,
        'masterleak': None,
        'leak': None,
        'masterwavelength': None,
        'wavelength': None,
        'extraction_test1': None,
        'extraction_test2': None,
        'extraction_test3': None,

        # Science
        'maketellu': None,
        'fittellu': None,
        'maketemplate': None,
        'ccf': None,
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
