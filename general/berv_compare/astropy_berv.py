#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-02-10

@author: cook
"""
from astropy.time import Time
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance, EarthLocation
import numpy as np
import warnings


# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def get_bervs(obstime, coordtime, ra, dec, pmra, pmde, plx,
             obs_long, obs_lat, obs_alt):
    """
    wrapper around get_berv (if required)

    :param obstime: observation time in JD (float)
    :param coordtime: time of the coordinates
    :param ra: ra in deg at coordtime
    :param dec: dec in deg at coordtime
    :param pmra: pmra in mas/yr at coordtime
    :param pmde: pmde in mas/yr at coordtime
    :param plx: plx in mas at coordtime (0 for None)
    :param obs_long: longitude of observatory
    :param obs_lat: latitude of observatory
    :param obs_alt: altitude of observatory
    :return: berv in km/s and bjd
    :return:
    """
    return get_berv(obstime, coordtime, ra, dec, pmra, pmde, plx,
                    obs_long, obs_lat, obs_alt)


def get_berv(obstime, coordtime, ra, dec, pmra, pmde, plx,
             obs_long, obs_lat, obs_alt):
    """
    :param obstime: observation time in JD (float)
    :param coordtime: time of the coordinates
    :param ra: ra in deg at coordtime
    :param dec: dec in deg at coordtime
    :param pmra: pmra in mas/yr at coordtime
    :param pmde: pmde in mas/yr at coordtime
    :param plx: plx in mas at coordtime (0 for None)
    :param obs_long: longitude of observatory
    :param obs_lat: latitude of observatory
    :param obs_alt: altitude of observatory
    :return: berv in km/s and bjd
    """
    # deal with distance
    if plx == 0:
        distance = None
    else:
        distance = Distance(parallax=plx * uu.mas)
    # need to propagate ra and dec to J2000
    coords = SkyCoord(ra=ra * uu.deg, dec=dec * uu.deg,
                      distance=distance,
                      pm_ra_cosdec=pmra * uu.mas/uu.yr,
                      pm_dec=pmde * uu.mas/uu.yr,
                      obstime=Time(coordtime, format='jd'))
    # get location
    loc = EarthLocation(lon=obs_long, lat=obs_lat, height=obs_alt)

    berv = coords.radial_velocity_correction(obstime=Time(obstime, format='jd'),
                                             location=loc)
    # bjd = np.nan

    return berv.value / 1000.0


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
