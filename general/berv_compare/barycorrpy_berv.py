#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2022-02-10

@author: cook
"""
import barycorrpy


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
    out1 = barycorrpy.get_BC_vel(JDUTC=obstime, zmeas=0.0, ra=ra, dec=dec,
                                 pmra=pmra, pmdec=pmde, px=plx,
                                 epoch=coordtime, lat=obs_lat,
                                 longi=obs_long, alt=obs_alt)
    # out2 = barycorrpy.utc_tdb.JDUTC_to_BJDTDB(obstime, ra=ra, dec=dec,
    #                                           pmra=pmra, pmdec=pmde, px=plx,
    #                                           epoch=coordtime, lat=obs_lat,
    #                                           longi=obs_long, alt=obs_alt)
    # return the bervs and bjds
    berv = out1[0] / 1000.0
    # bjd = out2[0]
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
