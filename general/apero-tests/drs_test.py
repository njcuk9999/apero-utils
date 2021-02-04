"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import os
from datetime import datetime
from typing import Optional

class DrsTest:
    def __init__(
        self,
        instrument: Optional[str] = None,
        drs_recipe = None,
        setup: Optional[str] = None,
    ):

        # get instrument
        self.instrument = instrument

        # Get setup path
        if setup is None:
            self.setup = os.environ["DRS_UCONFIG"]
        else:
            self.setup = setup

        # Set date at start of test
        self.date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # test date

        # Extract info from recipe
        if drs_recipe is not None:
            # Get recipe corresponding to test
            self.recipe = drs_recipe

            self.params = self.recipe.drs_params
            self.recipe_name = self.recipe.name.rstrip('.py')
            self.ismaster = self.recipe.master

        else:
            self.recipe = 'UnknownRecipe'
            self.params = None

        # TODO: Load log, outputs and calibdb in a nice way with APERO
