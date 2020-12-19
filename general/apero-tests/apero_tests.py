"""
Test abstract class and factory for APERO.

@author: vandal
"""
from abc import ABC, abstractmethod


class Test(ABC):
    """Test."""

    @property
    @abstractmethod
    def name(self):
        """Test full unique name."""

    # @property
    # @abstractmethod
    # def outlist(self):
    #     """outlist."""

    # @abstractmethod
    # def htmlgen(self):
    #     """htmlgen."""

    @abstractmethod
    def runtest(self):
        """Method called to run the tests on APERO products"""

    # @property
    # @abstractmethod
    # def outdict(self):
    #     """outdict."""
