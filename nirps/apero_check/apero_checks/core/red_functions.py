#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 17:03

@author: cook
"""
from typing import Any, Dict

import os

from apero_checks import raw_tests
from apero_checks import red_tests
from apero_checks.core import base
from apero_checks.core import io
from apero_checks.core import misc

# =============================================================================
# Define variables
# =============================================================================
# version, date, author
__VERSION__ = base.__VERSION__
__DATE__ = base.__DATE__
__AUTHOR__ = base.__AUTHOR__


# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def update_apero_profile(params: dict):
    """
    Update the apero profile with the correct paths

    :param profile: dict, the profile to update
    :return:
    """
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = params['apero profile']
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # make sure parameters is reloaded (and not cached)
    return constants.load(cache=False)


def get_apero_proxy_recipe(params: dict):
    """
    Get a proxy recipe for apero (when we need one for checks)
    It doesn't need any arguments and just replaces the standard recipe
    call

    :param params: parameters dictionary of constants (from cmd/yaml/param file)

    :return: Apero Recipe Class object
    """
    from apero.tools.module.testing import drs_dev
    # set up recipe definitions (overwrites default one)
    RMOD = drs_dev.RecipeDefinition(instrument=params['APERO_INSTRUMENT'])
    # define a recipe for this tool
    proxy_recipe = drs_dev.TmpRecipe()
    proxy_recipe.name = 'ProxyRecipe'
    proxy_recipe.instrument = params['APERO_INSTRUMENT']
    proxy_recipe.in_block_str = 'red'
    proxy_recipe.out_block_str = 'red'
    proxy_recipe.extension = 'fits'
    proxy_recipe.description = 'Proxy recipe when we need one for checks'
    proxy_recipe.kind = 'misc'
    # add recipe to recipe definition
    RMOD.add(proxy_recipe)
    # return the proxy recipe
    return proxy_recipe


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # print 'Hello World!'
    print("Hello World!")

# =============================================================================
# End of code
# =============================================================================
