"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import os
from datetime import datetime
from typing import Optional
from apero.core.instruments.spirou.recipe_definitions import recipes

RECIPE_DICT = dict(zip(list(map(lambda x: x.name, recipes)), recipes))


class DrsTest:
    def __init__(
        self,
        recipe: Optional[str] = None,
        setup: Optional[str] = None,
    ):

        # Get setup path
        if setup is None:
            self.setup = os.environ["DRS_UCONFIG"]
        else:
            self.setup = setup

        # Set date at start of test
        self.date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # test date


        # Extract info from recipe
        if recipe is not None:
            # Get recipe corresponding to test
            self.recipe = RECIPE_DICT[recipe]

            self.params = self.recipe.drs_params
            self.instrument = self.recipe.instrument
            self.recipe_name = self.recipe_name.rstrip('.py')
            self.ismaster = self.recipe.master

        else:
            self.recipe = None
            self.params = None
            self.instrument = 'UnknownInstrument'


        # TODO: Load log, outputs and calibdb in a nice way with APERO
