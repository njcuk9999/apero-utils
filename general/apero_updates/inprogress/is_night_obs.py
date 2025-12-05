#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2023-01-31 at 11:24

@author: cook
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from tqdm import tqdm
import warnings


# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def function1():
    return 0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    from astropy.io import fits
    import astropy.coordinates as coord
    import astropy.units as u
    from astropy.coordinates import get_sun, AltAz
    from astropy.time import Time

    file = 'NIRPS_2023-01-20T08_49_12_510_pp_e2dsff_tcorr_A.fits'

    # get the header
    h = fits.getheader(file)
    lat = h['BC_LAT']  # latitude
    lon = h['BC_LONG']  # longitude
    mjd = h['MJD-OBS']  # Modified Julian day

    sun_time = Time(mjd, format='mjd')  # UTC time
    loc = coord.EarthLocation(lon=lon * u.deg,
                              lat=lat * u.deg)

    altaz = AltAz(obstime=sun_time, location=loc)
    sun_elevation = 90 - get_sun(sun_time).transform_to(altaz).zen.value

    # Leval definition of twilight angles
    CIV_TWIL = sun_elevation < (-6)  # suggestion for Civil twilight keyword
    NAU_TWIL = sun_elevation < (-12)  # suggestion for Nautical twilight keyword
    AST_TWIL = sun_elevation < (-18)  # suggestion for Astronomical twilight keyword

    print('Civil twilight : {}\n'
          'Nautical twilight : {}\n'
          'Astron twilight : {}'.format(CIV_TWIL, NAU_TWIL, AST_TWIL))
    print('Sun elevation : {:.1f} deg'.format(sun_elevation))



# =============================================================================
# End of code
# =============================================================================
