#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-07-03 at 17:03

@author: cook
"""
import os
import sys
from typing import Any

from apero_checks.core import base

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
def update_apero_profile(params: dict) -> Any:
    """
    Update the apero profile with the correct paths

    :param profile: dict, the profile to update
    :return:
    """
    # use os to add DRS_UCONFIG to the path
    os.environ['DRS_UCONFIG'] = params['apero profile']
    # allow getting apero
    sys.path.append(params['apero install'])
    # ------------------------------------------------------------------
    # cannot import until after DRS_UCONFIG set
    from apero.base import base
    from apero.core import constants
    from apero.core.constants import param_functions
    from apero.core.utils import drs_startup
    # reload DPARAMS and IPARAMS
    base.DPARAMS = base.load_database_yaml()
    base.IPARAMS = base.load_install_yaml()
    # ------------------------------------------------------------------
    # invalidate cache
    param_functions.CONFIG_CACHE = dict()
    # make sure parameters is reloaded (and not cached)
    apero_params = constants.load(cache=False)
    # set apero pid
    apero_params['PID'], apero_params['DATE_NOW'] = drs_startup.assign_pid()
    # no inputs
    apero_params['INPUTS'] = dict()
    apero_params['OBS_DIR'] = None
    # return apero parameters
    return apero_params


def get_apero_proxy_recipe(apero_params: Any):
    """
    Get a proxy recipe for apero (when we need one for checks)
    It doesn't need any arguments and just replaces the standard recipe
    call

    :param params: parameters dictionary of constants (from cmd/yaml/param file)

    :return: Apero Recipe Class object
    """
    from apero.tools.module.testing import drs_dev
    from apero.core import constants
    # get pseudo constants
    pconst = constants.pload()
    # set up recipe definitions (overwrites default one)
    RMOD = drs_dev.RecipeDefinition(instrument=apero_params['INSTRUMENT'])
    # define a recipe for this tool
    proxy_recipe = drs_dev.TmpRecipe()
    proxy_recipe.name = 'ProxyRecipe'
    proxy_recipe.instrument = apero_params['INSTRUMENT']
    proxy_recipe.in_block_str = 'red'
    proxy_recipe.out_block_str = 'red'
    proxy_recipe.extension = 'fits'
    proxy_recipe.description = 'Proxy recipe when we need one for checks'
    proxy_recipe.kind = 'misc'
    proxy_recipe.recipemod = pconst.RECIPEMOD()
    # add recipe to recipe definition
    RMOD.add(proxy_recipe)
    # return the proxy recipe
    return proxy_recipe


def add_run_ini_params(apero_params, apero_recipe, run_file):
    from apero.core.utils import drs_startup
    # use the apero functions to update the apero_params
    apero_params, _ = drs_startup.read_runfile(apero_params, apero_recipe,
                                               runfile=run_file, rkind='run',
                                               log_overwrite=True)
    # return the update apero params
    return apero_params


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
