#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Example test of barycorrpy

Created on 2025-06-11 at 09:56

@author: cook
"""
import os
from typing import Union, Tuple

import barycorrpy
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as uu
from astropy.time import Time
from astropy.utils import iers


# =============================================================================
# Define variables
# =============================================================================

# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
class Star:
    def __init__(self, name: str):
        self.name = name
        self.ra = None
        self.dec = None
        self.epoch = 2451545.0  # J2000
        self.plx = 0.0 * uu.mas  # parallax in mas
        self.pmra = 0.0 * uu.mas / uu.yr  # proper motion in RA in mas/yr
        self.pmde = 0.0 * uu.mas / uu.yr   # proper motion in Dec in mas/yr
        self.rv = 0.0 * uu.km / uu.s  # radial velocity in km/s



def use_barycorrpy(times: np.ndarray, star: Star,
                   long: float, lat: float, alt: float,
                   leap_update: bool = False,
                   iersfile: Union[str, None] = None,
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Use the barycorrpy module to calculate BERV/BJD

    :param times: array of times in UTC [JD]
    :param ra: right ascension of the target [degrees]
    :param dec: declination of the target [degrees]
    :param lat: latitude of the observatory [degrees]
    :param long: longitude of the observatory [degrees]
    :param alt: altitude of the observatory [m]
    :param epoch: epoch of the target [default is J2000 = 2451545.0]
    :param plx: parallax of the target [mas] default = 0
    :param pmra: proper motion in RA of the target [mas/yr] default = 0
    :param pmde: proper motion in Dec of the target [mas/yr] default = 0
    :param rv: radial velocity of the target [km/s] -  default = 0
    :param leap_update: whether to update the leap seconds (default False)

    :return: two numpy arrays, the array of bervs for the times [km/s], and
             array of bjds [julien date]
    """
    # make barycorrpy directory an absolute path
    bc_dir = '.'
    # get args
    bkwargs = dict(ra=star.ra.to(uu.deg).value,
                   dec=star.dec.to(uu.deg).value,
                   epoch=star.epoch,
                   px=star.plx.to(uu.mas).value,
                   pmra=star.pmra.to(uu.mas/uu.yr).value,
                   pmdec=star.pmde.to(uu.mas/uu.yr).value,
                   longi=long, lat=lat, alt=alt,
                   leap_update=leap_update)
    # try to set iers file
    try:
        # iers.IERS_A_URL = iers_a_url
        iers_a_file = os.path.join(bc_dir, iersfile)
        iers.IERS.iers_table = iers.IERS_A.open(iers_a_file)
    except Exception as e:
        print('IERS_A_FILE Warning:' + str(e))
    # -------------------------------------------------------------------------
    print(f'Using barycorrpy on {len(times)} times for {star.name}')
    out1 = barycorrpy.get_BC_vel(JDUTC=times, zmeas=0.0, **bkwargs)
    out2 = barycorrpy.utc_tdb.JDUTC_to_BJDTDB(times, **bkwargs)
    # return the bervs and bjds
    bervs = out1[0] / 1000.0
    bjds = out2[0]
    return bervs, bjds


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # observatory
    obs_long = -70.731330408
    obs_lat =  -29.261165622
    obs_alt = 2400

    # Barnard's star (from astrometric database)
    gl699 = Star(name='GL699')
    gl699.ra = 269.4485025 * uu.deg
    gl699.dec = 4.739420051 * uu.deg
    gl699.epoch = 2457388.5
    gl699.pmra = -801.5509784 * uu.mas/uu.yr
    gl699.pmde = 10362.39421 * uu.mas/uu.yr
    gl699.plx = 546.9759397 * uu.mas
    gl699.rv = -110.11 * uu.km/uu.s

    # set up times
    times = Time(np.arange(60000, 61000, 0.1), format='mjd').jd
    # get the bervs and bjds
    bervs, bjds = use_barycorrpy(times, gl699, long=obs_long, lat=obs_lat,
                                 alt=obs_alt)
    # plot the results
    plt.plot(bjds, bervs, 'o-')

    plt.show()



# =============================================================================
# End of code
# =============================================================================
