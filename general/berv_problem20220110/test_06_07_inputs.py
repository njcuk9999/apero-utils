#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CODE DESCRIPTION HERE

Created on 1/10/22

@author: cook
"""
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
from astropy.io import fits
from astropy.time import Time
from barycorrpy import get_BC_vel
import glob
import numpy as np
import os
from typing import Any, Dict, Tuple
import warnings

# =============================================================================
# Define variables
# =============================================================================
# define paths to files
PATH_07 = '/space/spirou/LBL/SPIROU07179/science/'
PATH_06 = '/space/spirou/LBL/SPIROU06132/science/'

OBJECT = 'GL699'

apply_motion_06 = True
apply_motion_07 = True
shift_date = 'obs_time'

# =============================================================================
# Define functions
# =============================================================================
def apply_space_motion(ra: float, dec: float, pmra: float, pmde: float,
                       plx: float, rv: float, coordtime_jd: float,
                       requested_time_jd: float) -> Tuple[float, float]:
    # deal with distance
    if plx == 0:
        distance = None
    else:
        distance = Distance(parallax=plx * uu.mas)
    # need to propagate ra and dec to J2000
    coords1 = SkyCoord(ra=ra * uu.deg, dec=dec * uu.deg,
                       distance=distance,
                       pm_ra_cosdec=pmra * uu.mas / uu.yr,
                       pm_dec=pmde * uu.mas / uu.yr,
                       obstime=Time(coordtime_jd, format='jd'))
    # subtract time and get delta time (in jd)
    delta_time = requested_time_jd - coordtime_jd
    # get the coordinates in J2000
    with warnings.catch_warnings(record=True) as _:
        coords2 = coords1.apply_space_motion(dt=delta_time * uu.day)
    # extract the ra and dec from SkyCoords
    ra2 = coords2.ra.value
    dec2 = coords2.dec.value

    return ra2, dec2


def calc_berv(ra: float, dec: float, pmra: float, pmde: float,
              plx: float, rv: float, long: float, lat: float, alt: float,
              obstime_jd: float, coordtime_jd: float) -> float:
    # calculate berv with barycorrpy
    bervout = get_BC_vel(obstime_jd, ra=ra, dec=dec, pmra=pmra, pmdec=pmde,
                         px=plx, rv=rv, longi=long, lat=lat, alt=alt,
                         epoch=coordtime_jd)
    # return berv in km/s
    return bervout[0][0] / 1000.0


def params_from_06(filename: str) -> Dict[str, Any]:
    # read header
    hdr = fits.getheader(filename)

    # set up params
    params = dict()

    params['ra'] = hdr['BC_RA']
    params['dec'] = hdr['BC_DEC']
    params['pmra'] = hdr['BC_PMRA']
    params['pmde'] = hdr['BC_PMDE']
    params['plx'] = hdr['BC_PLX']
    params['rv'] = hdr['BC_RV']
    params['coord_time'] = hdr['BC_EPOCH']
    params['ALT'] = hdr['BC_ALT']
    params['LAT'] = hdr['BC_LAT']
    params['LONG'] = hdr['BC_LONG']
    params['obs_time'] = Time(hdr['MJDMID'], format='mjd').jd
    params['apero_berv'] = hdr['BERV']

    return params


def params_from_07(filename: str) -> Dict[str, Any]:
    # read header
    hdr = fits.getheader(filename)

    # set up params
    params = dict()

    params['ra'] = hdr['PP_RA']
    params['dec'] = hdr['PP_DEC']
    params['pmra'] = hdr['PP_PMRA']
    params['pmde'] = hdr['PP_PMDE']
    params['plx'] = hdr['PP_PLX']
    params['rv'] = hdr['PP_RV']
    params['coord_time'] = hdr['PP_EPOCH']
    params['ALT'] = hdr['BC_ALT']
    params['LAT'] = hdr['BC_LAT']
    params['LONG'] = hdr['BC_LONG']
    params['obs_time'] = Time(hdr['MJDMID'], format='mjd').jd
    params['apero_berv'] = hdr['BERV']

    return params



# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # get files from 0.7 and 0.6
    files07 = glob.glob(os.path.join(PATH_07, OBJECT, '*AB.fits'))
    files06 = glob.glob(os.path.join(PATH_06, OBJECT, '*AB.fits'))

    basenames07 = np.array(list(map(lambda x: os.path.basename(x), files07)))
    basenames06 = np.array(list(map(lambda x: os.path.basename(x), files06)))

    for it in range(len(files07[:10])):

        # find 0.6 file for 0.7
        mask = basenames06 == basenames07[it]
        # skip if no match
        if np.sum(mask) == 0:
            continue
        pos = np.where(mask)[0][0]

        # get parameters from 0.6 and 0.7
        params07 = params_from_07(files07[it])
        params06 = params_from_06(files06[pos])


        # set pmra,pmde,plx the same
        # params07['pmra'] = params06['pmra']
        # params07['pmde'] = params06['pmde']
        # params07['plx'] = params06['plx']

        # set 0.7 coord_time
        # params07['coord_time'] = Time(2015.5, format='decimalyear').jd

        # push both to the same date
        if shift_date in params06 and shift_date in params07:
            requested_time06 = params06[shift_date]
            requested_time07 = params07[shift_date]
        else:
            requested_time06 = Time(shift_date, format='decimalyear').jd
            requested_time07 = Time(shift_date, format='decimalyear').jd

        ra07, dec07 = apply_space_motion(params07['ra'], params07['dec'],
                                         params07['pmra'], params07['pmde'],
                                         params07['plx'], params07['rv'],
                                         coordtime_jd=params07['coord_time'],
                                         requested_time_jd=requested_time06)
        ra06, dec06 = apply_space_motion(params06['ra'], params06['dec'],
                                         params06['pmra'], params06['pmde'],
                                         params06['plx'], params06['rv'],
                                         coordtime_jd=params06['coord_time'],
                                         requested_time_jd=requested_time07)

        if apply_motion_06:
            params06['ra'] = ra06
            params06['dec'] = dec06
            params06['coord_time'] = requested_time06

        if apply_motion_07:
            params07['ra'] = ra07
            params07['dec'] = dec07
            params07['coord_time'] = requested_time07

        berv_new06 = calc_berv(params06['ra'], params06['dec'],
                               params06['pmra'], params06['pmde'],
                               params06['plx'], params06['rv'],
                               params06['LONG'], params06['LAT'],
                               params06['ALT'], obstime_jd=params06['obs_time'],
                               coordtime_jd=params06['coord_time'])

        berv_new07 = calc_berv(params07['ra'], params07['dec'],
                               params07['pmra'], params07['pmde'],
                               params07['plx'], params07['rv'],
                               params07['LONG'], params07['LAT'],
                               params07['ALT'], obstime_jd=params07['obs_time'],
                               coordtime_jd=params07['coord_time'])

        dra = (ra07-ra06)*3600
        ddec = (dec07-dec06)*3600

        dberv06 = berv_new06 - params06['apero_berv']
        dberv07 = berv_new07 - params07['apero_berv']

        dtime06 = Time(params06["coord_time"], format='jd').decimalyear
        dtime07 = Time(params07["coord_time"], format='jd').decimalyear

        print(f'{basenames06[it]}'
              f'\n\tdra={dra} arcsec\tddec={ddec} arcsec'
              f'\n\t0.6: apply_motion={apply_motion_06}\tberv={berv_new06} km/s\tdberv={dberv06} km/s'
              f'\tcoord_time={dtime06}'    
              f'\n\t0.7: apply_motion={apply_motion_07}\tberv={berv_new07} km/s\tdberv={dberv07} km/s'
              f'\tcoord_time={dtime07}')






# =============================================================================
# End of code
# =============================================================================
