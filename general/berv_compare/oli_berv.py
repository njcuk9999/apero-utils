#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-02-10

@author: cook
"""
import numpy as np

# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define oli functions
# =============================================================================
def solQuick(JD):
    # in: julian day of obs
    # out: Longitude of the sun. Super-crude, OK within 5deg.
    JD0 = 2459659.1368  # 20 mars, spring
    au = 150.e6  # 1AU
    sun_long = (360. * (JD - JD0) / 365.25) % 360.
    return sun_long


def EarthVel(sunLong):
    # in: longitude of the sun
    # out: velocity vector of the Earth in equatorial coordinates [km/s]

    obliquityr = np.radians(23.44)  # obliquity of ecliptic
    co = np.cos(obliquityr)
    so = np.sin(obliquityr)

    au = 1.495978707e8  # AU in km
    year = 365.25 * 86400.  # year in sec
    eVel = 2. * np.pi * au / year  # in km/s

    eLon = sunLong
    eEc = eVel * np.array(
        [np.sin(np.radians(eLon)), -np.cos(np.radians(eLon)), 0. * eLon])  # velocity in ecliptic coord
    eEq = np.array([eEc[0], co * eEc[1] - so * eEc[2], so * eEc[1] + co * eEc[2]])  # converted to equatorial
    return eEq


def olibarycorr(JD, star):
    # JD julian day of obs
    # star: [RA, Dec] in degrees
    # out: barycentric correction in km/s

    sR = np.radians(star)
    sEq = np.array(
        [np.cos(sR[0]) * np.cos(sR[1]), np.sin(sR[0]) * np.cos(sR[1]), np.sin(sR[1])])  # equatorial vec to star

    sunLong = solQuick(JD)
    eEq = EarthVel(sunLong)  # e
    v = np.dot(sEq, eEq)  # scalar product between the 2 vectors

    return v


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
    :param ra: ra at 2000.0
    :param dec: dec at 2000.0
    :param pmra: not used
    :param pmde: not used
    :param plx: not used
    :param obs_long: longitude of observatory
    :param obs_lat: latitude of observatory
    :param obs_alt: altitude of observatory
    :return: berv in km/s and bjd
    """
    _ = pmra, pmde, plx
    _ = obs_long, obs_lat, obs_alt
    _ = coordtime

    berv = olibarycorr(obstime, np.array([ra,dec]))

    # bjd = None

    return berv


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # print hello world
    print('Hello World')

# =============================================================================
# End of code
# =============================================================================
