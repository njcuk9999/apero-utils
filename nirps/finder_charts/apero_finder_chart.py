import argparse
import warnings
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.visualization.wcsaxes import WCSAxes
from astropy.wcs import WCS
from astroquery.utils.tap.core import TapPlus
from dateutil.parser import parse
from matplotlib.transforms import Affine2D
from tqdm import tqdm

from apero import lang
from apero.core import constants
from apero.core.core import drs_database
from apero.core.core import drs_log
from apero.core.utils import drs_startup

# =============================================================================
# Define variables
# =============================================================================
# get WLOG
WLOG = drs_log.wlog
# columns
ASTROM_COLS = ['OBJNAME', 'RA_DEG', 'DEC_DEG', 'PMRA', 'PMDE', 'PLX', 'RV',
               'EPOCH']
# List of surveys to plot
SURVEYS = dict()
SURVEYS['DSS'] = ['DSS']
SURVEYS['2MASS'] = ['2MASS-J', '2MASS-H', '2MASS-K']
# define the time in years before and after "obs date" to plot trail for
DTMIN = -5
DTMAX = 5
DTVAL = 1
# FIELD OF VIEW
FIELD_OF_VIEW = 2 * uu.arcmin
# RADIUS
RADIUS = FIELD_OF_VIEW / np.sqrt(2)
# PIXEL SCALE
PIXEL_SCALE = 0.3 * uu.arcsec / uu.pixel
# FWHM
FWHM = 1.0 * uu.arcsec

# Transform the image
TRANSFORM_ROTATE = 45 * uu.deg

# 1-sigma magnitude limit G
SIGMA_LIMIT = dict()
SIGMA_LIMIT['G'] = 22
SIGMA_LIMIT['J'] = 17

# URL for Gaia
GAIA_URL = 'https://gea.esac.esa.int/tap-server/tap'
# noinspection SqlNoDataSourceInspection,SqlDialectInspection
GAIA_QUERY = """
SELECT
   source_id,ra,dec,parallax,pmra,pmdec,phot_g_mean_mag,phot_rp_mean_mag,
   phot_bp_mean_mag, phot_g_mean_flux, phot_rp_mean_flux, phot_bp_mean_flux
   FROM gaiadr3.gaia_source
   WHERE 
      1=CONTAINS(POINT('ICRS', ra, dec),
                 CIRCLE('ICRS', {ra}, {dec}, {radius}))
"""
GAIA_EPOCH = 2016.0
# URL for 2MASS
TWOMASS_URL = ''
# URL for 2MASS+Gaia crossmatch
GAIA_TWOMASS_URL = 'https://gaia.obspm.fr/tap-server/tap'


# =============================================================================
# Define functions
# =============================================================================
def get_args():
    """
    Define the command line arguments

    :return: argparse namespace object containing arguments
    """
    parser = argparse.ArgumentParser(description='Rene\'s magic trigger')
    # add obs dir
    parser.add_argument('objname', type=str, help='The object name')
    parser.add_argument('--date', type=str, default='None',
                        help='The date of observation')
    # load arguments with parser
    args = parser.parse_args()
    # deal with date
    if args.date == 'None':
        date = Time.now()
    else:
        try:
            date = Time(parse(str(args.date)))
        except Exception as _:
            raise ValueError(f'Cannot parse --date={args.date}')

    args.date = date

    # return arguments
    return args


def propagate_coords(objdata, obs_time: Time):
    # deal with distance
    if objdata['PLX'] == 0:
        distance = None
    else:
        distance = Distance(parallax=objdata['PLX'] * uu.mas)
    # need to propagate ra and dec to J2000
    coords = SkyCoord(ra=objdata['RA_DEG'] * uu.deg,
                      dec=objdata['DEC_DEG'] * uu.deg,
                      distance=distance,
                      pm_ra_cosdec=objdata['PMRA'] * uu.mas / uu.yr,
                      pm_dec=objdata['PMDE'] * uu.mas / uu.yr,
                      obstime=Time(objdata['EPOCH'], format='jd'))
    # work out the delta time between epoch and J2000.0
    jepoch = Time(objdata['EPOCH'], format='jd')
    delta_time = (obs_time.jd - jepoch.jd) * uu.day
    # get a list of times to work out the coordiantes at
    delta_times = delta_time + np.arange(DTMIN, DTMAX + 1, DTVAL) * uu.year
    # get the coordinates
    with warnings.catch_warnings(record=True) as _:
        new_coords = coords.apply_space_motion(dt=delta_times)
        new_times = jepoch + delta_times
        # center the image on the current coordinates
        curr_coords = coords.apply_space_motion(dt=delta_time)
        curr_time = obs_time
    # return some stuff
    return new_coords, new_times, curr_coords, curr_time


# noinspection PyUnresolvedReferences
def setup_wcs(image_shape, cent_coords, pixel_scale):
    # need to set up a wcs
    naxis2, naxis1 = image_shape
    pix_scale = pixel_scale.to(uu.deg / uu.pixel).value
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [naxis1 / 2, naxis2 / 2]
    wcs.wcs.cdelt = np.array([-pix_scale, pix_scale])
    wcs.wcs.crval = [cent_coords.ra.deg, cent_coords.dec.deg]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs


def get_gaia_sources(coords: SkyCoord, obstime: Time) -> Dict[str, List[float]]:
    # define the gaia Tap instance
    gaia = TapPlus(url=GAIA_URL)
    # launch the query
    job = gaia.launch_job(GAIA_QUERY.format(ra=coords.ra.deg,
                                            dec=coords.dec.deg,
                                            radius=RADIUS.to(uu.deg).value))
    # get the query results
    table = job.get_results()
    # get the gaia time
    gaia_time = Time(GAIA_EPOCH, format='decimalyear')
    # storage for sources
    sources = dict()
    sources['ra'] = []
    sources['dec'] = []
    sources['G'] = []
    sources['Rp'] = []
    sources['Bp'] = []
    sources['J'] = []
    # get parallax mask
    parallax_mask = table['parallax'].mask
    # -------------------------------------------------------------------------
    # get all stars in the region where they would be now
    for row in range(len(table)):
        # get table row
        table_row = table[row]
        # deal with distances
        if parallax_mask[row]:
            distance = None
        elif table_row['parallax'] <= 0:
            distance = None
        else:
            distance = Distance(parallax=table_row['parallax'] * uu.mas)
        # need to propagate ra and dec to J2000
        coords = SkyCoord(ra=table_row['ra'] * uu.deg,
                          dec=table_row['dec'] * uu.deg,
                          distance=distance,
                          pm_ra_cosdec=table_row['pmra'] * uu.mas / uu.yr,
                          pm_dec=table_row['pmdec'] * uu.mas / uu.yr,
                          obstime=gaia_time)
        # work out the delta time between epoch and J2000.0
        delta_time = (obstime.jd - gaia_time.jd) * uu.day
        # center the image on the current coordinates
        with warnings.catch_warnings(record=True) as _:
            curr_coords = coords.apply_space_motion(dt=delta_time)

        sources['ra'].append(curr_coords.ra.deg)
        sources['dec'].append(curr_coords.dec.deg)
        sources['G'].append(table_row['phot_g_mean_mag'])
        sources['Rp'].append(table_row['phot_rp_mean_mag'])
        sources['Bp'].append(table_row['phot_bp_mean_mag'])
        sources['J'].append(np.nan)
    # -------------------------------------------------------------------------
    # convert all to numpy arrays
    sources['ra'] = np.array(sources['ra'])
    sources['dec'] = np.array(sources['dec'])
    sources['G'] = np.array(sources['G'])
    sources['Rp'] = np.array(sources['Rp'])
    sources['Bp'] = np.array(sources['Bp'])
    sources['J'] = np.array(sources['J'])
    # -------------------------------------------------------------------------
    return sources


def out_of_bounds(frame, coord):
    xlims = frame.get_xlim()
    ylims = frame.get_ylim()

    lims = [xlims, ylims]
    values = [coord.ra.value, coord.dec.value]

    for value, lim in zip(values, lims):
        if lim[0] > lim[1]:
            if value < lim[1]:
                return True
            if value > lim[0]:
                return True
        else:
            if value < lim[0]:
                return True
            if value > lim[1]:
                return True
    return False


# =============================================================================
# Start of code
# =============================================================================
if __name__ == '__main__':
    # get the parameter dictionary of constants from apero
    params = constants.load()
    pconst = constants.pload()
    # Get the text types
    textentry = lang.textentry
    # load arguments
    args = get_args()
    # set apero pid
    params['PID'], params['DATE_NOW'] = drs_startup.assign_pid()
    # -------------------------------------------------------------------------
    # print progress
    print('Getting object from database')
    # load object database
    objdbm = drs_database.AstrometricDatabase(params)
    objdbm.load_db()
    # get the clean / alias-safe version of the object name
    objname, _ = objdbm.find_objname(pconst, args.objname)
    # get the data for this object
    objdata = objdbm.get_entries(','.join(ASTROM_COLS),
                                 condition=f'OBJNAME="{objname}"').iloc[0]
    # -------------------------------------------------------------------------
    print('Progating coords')
    pout = propagate_coords(objdata, args.date)
    all_coords, all_times, obs_coords, obs_time = pout
    # -------------------------------------------------------------------------
    print('Getting Gaia sources for this region')
    # get all Gaia data around the current coordinates
    gaia_sources = get_gaia_sources(obs_coords, obs_time)
    # -------------------------------------------------------------------------
    # TODO: Get 2mass mags for this gaia sources (or all available)
    #   Put all sources not in 2mass at mag limit
    # -------------------------------------------------------------------------
    print('Seeding image')

    # number of pixels in each direction
    npixel = int(FIELD_OF_VIEW.to(uu.deg) // (PIXEL_SCALE*uu.pixel).to(uu.deg))


    wcs = setup_wcs((npixel, npixel), obs_coords, PIXEL_SCALE)

    image = np.random.normal(size=(npixel, npixel), scale = 1.0, loc = 0)

    nsig_psf = np.array(10**((SIGMA_LIMIT['G'] - gaia_sources['G']) / 2.5))

    # convert all coords to pixel coords
    x_sources, y_sources = wcs.all_world2pix(gaia_sources['ra'],
                                             gaia_sources['dec'], 0)
    # create a grid of all pixels
    y, x = np.mgrid[0:npixel, 0:npixel]

    ew = FWHM.value / (2 * np.sqrt(2 * np.log(2)))

    for i in tqdm(range(len(x_sources))):
        # core of PSF
        image += nsig_psf[i] * np.exp(-((x - x_sources[i]) ** 2 + (y - y_sources[i]) ** 2) / (2 * ew ** 2))
        # halo of PSF
        image += nsig_psf[i] * np.exp(-((x - x_sources[i]) ** 2 + (y - y_sources[i]) ** 2) / (2 * (ew * 2) ** 2)) * 1e-3
        # spike in y
        image += nsig_psf[i] * np.exp(-np.abs((x - x_sources[i]) ** 2 + 0.1 * np.abs(y - y_sources[i]) ** 2) / (2 * ew ** 2)) * 1e-2
        # spike in x
        image += nsig_psf[i] * np.exp(-(0.1 * np.abs(x - x_sources[i]) ** 2 + np.abs(y - y_sources[i]) ** 2) / (2 * ew ** 2)) * 1e-2

    # -------------------------------------------------------------------------
    # plot figure
    # -------------------------------------------------------------------------
    # get the extent of the image
    extent = [wcs.wcs.crval[0] - wcs.wcs.cdelt[0] * image.shape[1] / 2,
              wcs.wcs.crval[0] + wcs.wcs.cdelt[0] * image.shape[1] / 2,
              wcs.wcs.crval[1] - wcs.wcs.cdelt[1] * image.shape[0] / 2,
              wcs.wcs.crval[1] + wcs.wcs.cdelt[1] * image.shape[0] / 2]
    dra = extent[1] - extent[0]
    ddec = extent[3] - extent[2]



    transform = Affine2D()
    transform.rotate(TRANSFORM_ROTATE.to(uu.rad).value)

    # Set up metadata dictionary
    coord_meta = {}
    coord_meta['name'] = 'ra', 'dec'
    coord_meta['type'] = 'longitude', 'latitude'
    coord_meta['wrap'] = 180, None
    coord_meta['unit'] = uu.deg, uu.deg
    coord_meta['format_unit'] = None, None

    plt.close()
    fig = plt.figure()
    # frame = WCSAxes(fig, extent, aspect='equal',
    #                 transform=transform, coord_meta=coord_meta)
    # fig.add_axes(frame)
    frame = plt.subplot(projection=wcs)


    frame.imshow(image, origin='lower', vmin=-3, vmax=20,
                 extent=extent, cmap='gist_heat')

    # plot the trail
    frame.plot(all_coords.ra, all_coords.dec, color='yellow', ls='None',
               marker='+', lw=2, markersize=10, alpha=0.5)
    # plot the year labels
    for it, coord in enumerate(all_coords):
        # check whether text is out of bounds
        if out_of_bounds(frame, coord):
            continue
        # text is offset by 1% of the ra range
        frame.text(coord.ra.value + 0.01 * dra, coord.dec.value,
                   f'{all_times[it].decimalyear:.1f}', color='yellow', alpha=0.5,
                   horizontalalignment='left', verticalalignment='center',)

    # plot the current position as an open blue circle
    frame.plot(obs_coords.ra, obs_coords.dec, marker='o',
               color='yellow', ms=20, mfc='none')
    # plot a grid in white
    frame.grid(color='white', ls='--', lw=1, alpha=0.5)

    plt.show()




    stop


