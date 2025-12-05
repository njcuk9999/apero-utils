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
from astropy.table import Table
from barycorrpy import get_BC_vel
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
from typing import Any, Dict, Tuple
import warnings

# =============================================================================
# Define variables
# =============================================================================
# define paths to files
PATH_07 = '/space/spirou/LBL/SPIROU07179/science/'
PATH_06 = '/space/spirou/LBL/SPIROU06132/science/'
workspace = '/data/nas1/spirou/cook/berv_problem20220110'
OBJECT = 'GL699'

apply_motion_06 = True
apply_motion_07 = True

reset = False
verbose = False
# shift_date = 'obs_time'
shift_date = 2015.5

table06name = f'berv_06_table_{OBJECT}.fits'
table07name = f'berv_07_table_{OBJECT}.fits'

epoch06 = '2015.0'
epoch07 = '2000.0'

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

    params['abspath'] = filename
    params['basename'] = os.path.basename(filename)
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
    params['mjdmid'] = hdr['MJDMID']
    params['apero_berv'] = hdr['BERV']

    return params


def params_from_07(filename: str) -> Dict[str, Any]:
    # read header
    hdr = fits.getheader(filename)

    # set up params
    params = dict()

    params['abspath'] = filename
    params['basename'] = os.path.basename(filename)
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
    params['mjdmid'] = hdr['MJDMID']
    params['apero_berv'] = hdr['BERV']

    return params



# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # generate table names
    abstbl06 = os.path.join(workspace, table06name)
    abstbl07 = os.path.join(workspace, table07name)

    if reset:
        if os.path.exists(abstbl06):
            os.remove(abstbl06)
        if os.path.exists(abstbl07):
            os.remove(abstbl07)


    if not os.path.exists(abstbl06) or not os.path.exists(abstbl07):
        # get files from 0.7 and 0.6
        files07 = glob.glob(os.path.join(PATH_07, OBJECT, '*AB.fits'))
        files06 = glob.glob(os.path.join(PATH_06, OBJECT, '*AB.fits'))
        basenames07 = np.array(list(map(lambda x: os.path.basename(x),
                                        files07)))
        basenames06 = np.array(list(map(lambda x: os.path.basename(x),
                                        files06)))
        # get a list of files
        storage06 = dict()
        storage07 = dict()

        print('Opening all headers...')

        for it in tqdm(range(len(files07))):

            # find 0.6 file for 0.7
            mask = basenames06 == basenames07[it]
            # skip if no match
            if np.sum(mask) == 0:
                continue
            pos = np.where(mask)[0][0]
            # get parameters from 0.6 and 0.7
            params07 = params_from_07(files07[it])
            params06 = params_from_06(files06[pos])
            # push into storage
            for key in params06:
                if key in storage06:
                    storage06[key].append(params06[key])
                else:
                    storage06[key] = [params06[key]]
            for key in params07:
                if key in storage07:
                    storage07[key].append(params07[key])
                else:
                    storage07[key] = [params07[key]]

        # convert storage to tbl and save
        tbl06 = Table()
        for col in storage06:
            tbl06[col] = storage06[col]
        tbl07 = Table()
        for col in storage07:
            tbl07[col] = storage07[col]

        # write to disk
        tbl06.write(abstbl06, format='fits')
        tbl07.write(abstbl07, format='fits')
    else:
        tbl06 = Table.read(abstbl06)
        tbl07 = Table.read(abstbl07)

    tbl06['new_berv'] = np.zeros(len(tbl06)).astype(float)
    tbl07['new_berv'] = np.zeros(len(tbl06)).astype(float)

    if apply_motion_06:
        str06 = f'recal_coord[{shift_date}]'
    else:
        str06 = f'recal_coord[{epoch06}]'
    if apply_motion_07:
        str07 = f'recal_coord[{shift_date}]'
    else:
        str07 = f'recal_coord[{epoch07}]'


    hdr06 = f'hdr_coord[{epoch06}]'
    hdr07 = f'hdr_coord[{epoch07}]'

    print('Calculating BERV / Shifts')
    for it in tqdm(range(len(tbl07))):

        basename = tbl06['basename'][it]

        # set pmra,pmde,plx the same
        # params07['pmra'] = params06['pmra']
        # params07['pmde'] = params06['pmde']
        # params07['plx'] = params06['plx']

        # set 0.7 coord_time
        # params07['coord_time'] = Time(2015.5, format='decimalyear').jd

        # push both to the same date
        if shift_date in tbl07.colnames and shift_date in tbl07.colnames:
            requested_time06 = tbl06[shift_date][it]
            requested_time07 = tbl07[shift_date][it]
        else:
            requested_time06 = Time(shift_date, format='decimalyear').jd
            requested_time07 = Time(shift_date, format='decimalyear').jd

        ra07, dec07 = apply_space_motion(tbl07['ra'][it], tbl07['dec'][it],
                                         tbl07['pmra'][it], tbl07['pmde'][it],
                                         tbl07['plx'][it], tbl07['rv'][it],
                                         coordtime_jd=tbl07['coord_time'][it],
                                         requested_time_jd=requested_time06)
        ra06, dec06 = apply_space_motion(tbl06['ra'][it], tbl06['dec'][it],
                                         tbl06['pmra'][it], tbl06['pmde'][it],
                                         tbl06['plx'][it], tbl06['rv'][it],
                                         coordtime_jd=tbl06['coord_time'][it],
                                         requested_time_jd=requested_time07)

        if apply_motion_06:
            tbl06['ra'][it] = ra06
            tbl06['dec'][it] = dec06
            tbl06['coord_time'][it] = requested_time06

        if apply_motion_07:
            tbl07['ra'][it] = ra07
            tbl07['dec'][it] = dec07
            tbl07['coord_time'][it] = requested_time07

        berv_new06 = calc_berv(tbl06['ra'][it], tbl06['dec'][it],
                               tbl06['pmra'][it], tbl06['pmde'][it],
                               tbl06['plx'][it], tbl06['rv'][it],
                               tbl06['LONG'][it], tbl06['LAT'][it],
                               tbl06['ALT'][it],
                               obstime_jd=tbl06['obs_time'][it],
                               coordtime_jd=tbl06['coord_time'][it])

        berv_new07 = calc_berv(tbl07['ra'][it], tbl07['dec'][it],
                               tbl07['pmra'][it], tbl07['pmde'][it],
                               tbl07['plx'][it], tbl07['rv'][it],
                               tbl07['LONG'][it], tbl07['LAT'][it],
                               tbl07['ALT'][it],
                               obstime_jd=tbl07['obs_time'][it],
                               coordtime_jd=tbl07['coord_time'][it])

        dra = (ra07-ra06)*3600
        ddec = (dec07-dec06)*3600

        dberv06 = berv_new06 - tbl06['apero_berv'][it]
        dberv07 = berv_new07 - tbl07['apero_berv'][it]

        dtime06 = Time(tbl06["coord_time"][it], format='jd').decimalyear
        dtime07 = Time(tbl07["coord_time"][it], format='jd').decimalyear

        if verbose:
            print(f'{basename}'
                  f'\n\tdra={dra} arcsec\tddec={ddec} arcsec'
                  f'\n\t0.6: apply_motion={apply_motion_06}'
                  f'\tberv={berv_new06} km/s\tdberv={dberv06} km/s'
                  f'\tcoord_time={dtime06}'    
                  f'\n\t0.7: apply_motion={apply_motion_07}'
                  f'\tberv={berv_new07} km/s\tdberv={dberv07} km/s'
                  f'\tcoord_time={dtime07}')
        # update table
        tbl06['new_berv'][it] = berv_new06
        tbl07['new_berv'][it] = berv_new07


    plt.close()
    fig, frames = plt.subplots(nrows=6, sharex='all')

    frames[0].plot(tbl07['mjdmid'],
                   1000*(tbl06['apero_berv'] - tbl07['apero_berv']),
                   color='purple', marker='.', ls='None',
                   label=f'0.6_{hdr06} - 0.7_{hdr07}')
    frames[1].plot(tbl06['mjdmid'], 1000*(tbl06['apero_berv']), 'g.',
                   alpha=0.4, label=f'0.6_{hdr06}')
    frames[1].plot(tbl06['mjdmid'], 1000*(tbl06['new_berv']), 'r.',
                   alpha=0.4, label=f'0.6_{str06}')
    frames[2].plot(tbl06['mjdmid'],
                   1000*(tbl06['apero_berv'] - tbl06['new_berv']), 'b.',
                   label=f'0.6_{hdr06} - 0.6_{str06}')
    frames[3].plot(tbl07['mjdmid'], 1000*(tbl07['apero_berv']), 'g.',
                   alpha=0.4, label=f'0.7_{hdr07}')
    frames[3].plot(tbl07['mjdmid'], 1000*(tbl07['new_berv']), 'r.',
                   alpha=0.4, label=f'0.7_{str07}')
    frames[4].plot(tbl07['mjdmid'],
                   1000*(tbl07['apero_berv'] - tbl07['new_berv']), 'b.',
                   label=f'0.7_{hdr07} - 0.7_{str07}')
    frames[5].plot(tbl07['mjdmid'],
                   1000*(tbl06['apero_berv'] - tbl07['new_berv']),
                   color='orange', marker='.', ls='None',
                   label=f'0.6_{hdr06} - 0.7_{str07}')

    frames[0].set(ylabel='dBERV [m/s]')
    frames[0].legend(loc=0)

    frames[1].set(ylabel='BERV [m/s]')
    frames[2].set(ylabel='dBERV [m/s]')
    frames[1].legend(loc=0)
    frames[2].legend(loc=0)

    frames[3].set(ylabel='BERV [m/s]')
    frames[4].set(ylabel='dBERV [m/s]')
    frames[3].legend(loc=0)
    frames[4].legend(loc=0)

    frames[5].set(xlabel='MJD', ylabel='dBERV [m/s]')
    frames[5].legend(loc=0)

    frames[0].ticklabel_format(style='plain')
    frames[1].ticklabel_format(style='plain')
    frames[2].ticklabel_format(style='plain')
    frames[3].ticklabel_format(style='plain')
    frames[4].ticklabel_format(style='plain')
    frames[5].ticklabel_format(style='plain')

    plt.subplots_adjust(hspace=0, left=0.075, bottom=0.05, top=0.99,
                        right=0.99)
    plt.show()


# =============================================================================
# End of code
# =============================================================================
