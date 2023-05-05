import argparse
import warnings
from typing import Dict, List

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.utils.tap.core import TapPlus
from dateutil.parser import parse
from scipy.optimize import minimize_scalar
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
# -----------------------------------------------------------------------------
# define the time in years before and after "obs date" to plot trail for
DTMIN = -6
DTMAX = 6
DTVAL = 3
# -----------------------------------------------------------------------------
# FIELD OF VIEW
FIELD_OF_VIEW = dict()
FIELD_OF_VIEW['G'] = 77.6 * uu.arcsec
FIELD_OF_VIEW['J'] = 77.6 * uu.arcsec
# SCALE FACTOR (size to be plotted to compared to field of view)
SCALE_FACTOR = 1.5
# RADIUS
RADIUS = dict()
RADIUS['G'] = FIELD_OF_VIEW['G'] / np.sqrt(2)
RADIUS['J'] = FIELD_OF_VIEW['J'] / np.sqrt(2)
# PIXEL SCALE
PIXEL_SCALE = dict()
PIXEL_SCALE['G'] = 0.1514 * uu.arcsec / uu.pixel
PIXEL_SCALE['J'] = 0.1514 * uu.arcsec / uu.pixel

#
# #FIELD OF VIEW
# FIELD_OF_VIEW = dict()
# FIELD_OF_VIEW['G'] = 2 * uu.arcmin
# FIELD_OF_VIEW['J'] = 2 * uu.arcmin
# # RADIUS
# RADIUS = dict()
# RADIUS['G'] = FIELD_OF_VIEW['G'] / np.sqrt(2)
# RADIUS['J'] = FIELD_OF_VIEW['J'] / np.sqrt(2)
# # PIXEL SCALE
# PIXEL_SCALE = dict()
# PIXEL_SCALE['G'] = 0.3 * uu.arcsec / uu.pixel
# PIXEL_SCALE['J'] = 0.3 * uu.arcsec / uu.pixel


# Transform the image
TRANSFORM_ROTATE = dict()
TRANSFORM_ROTATE['G'] = 90 * uu.deg - 37.4 * uu.deg
TRANSFORM_ROTATE['J'] = 90 * uu.deg - 37.4 * uu.deg


# FWHM
FWHM = dict()
FWHM['G'] = 2.0 * uu.arcsec
FWHM['J'] = 2.0 * uu.arcsec
# 1-sigma magnitude limit G
SIGMA_LIMIT = dict()
SIGMA_LIMIT['G'] = 22
SIGMA_LIMIT['J'] = 17
SIGMA_LIMIT['H'] = 17
SIGMA_LIMIT['K'] = 17
# mag limit above 1-sigma limit for null sources
MAG_LIMIT = -1
# max proper motion to search for
MAX_PM = 11 * uu.arcsec / uu.year
# compass size
COMPASS_SIZE = 10 * uu.arcsec
# FLIP X
FLIP_X = True
# FLIP Y
FLIP_Y = True


# -----------------------------------------------------------------------------
# Gaia specific
# -----------------------------------------------------------------------------
# URL for Gaia
GAIA_URL = 'https://gea.esac.esa.int/tap-server/tap'
# noinspection SqlNoDataSourceInspection,SqlDialectInspection
GAIA_QUERY = """
SELECT
    source_id,ra,dec,parallax,pmra,pmdec,phot_g_mean_mag,phot_rp_mean_mag,
    phot_bp_mean_mag, phot_g_mean_flux, phot_rp_mean_flux, phot_bp_mean_flux
FROM gaiadr3.gaia_source
WHERE 
    1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius}))
"""
# define epoch
GAIA_EPOCH = 2016.0
# -----------------------------------------------------------------------------
# Gaia specific
# -----------------------------------------------------------------------------
# define mean epoch
TMASS_EPOCH = 2000.0
# URL for Gaia
TMASS_URL = 'https://irsa.ipac.caltech.edu/TAP'
# noinspection SqlNoDataSourceInspection,SqlDialectInspection
TMASS_QUERY = """
SELECT
    {TMASS_COLS}
FROM fp_2mass.fp_psc
WHERE 
    1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius}))
"""
TMASS_COLS = 'ra,dec,j_m,h_m,k_m,jdate'
# CROSSMATCH RADIUS (SHOULD BE SMALL)
TMASS_RADIUS = 2 * uu.arcsec


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
    with warnings.catch_warnings(record=True) as _:
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
def setup_wcs(image_shape, cent_coords, pixel_scale, rotation,
              flip_x, flip_y):
    # need to set up a wcs
    naxis2, naxis1 = image_shape
    pix_scale = pixel_scale.to(uu.deg / uu.pixel).value
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [naxis1 / 2, naxis2 / 2]

    if flip_y and flip_x:
        wcs.wcs.cdelt = np.array([pix_scale, -pix_scale])
    elif flip_x:
        wcs.wcs.cdelt = np.array([pix_scale, pix_scale])
    elif flip_y:
        wcs.wcs.cdelt = np.array([pix_scale, -pix_scale])
    else:
        wcs.wcs.cdelt = np.array([-pix_scale, pix_scale])
    wcs.wcs.crval = [cent_coords.ra.deg, cent_coords.dec.deg]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    # add a rotation
    wcs.wcs.pc = np.array([[np.cos(rotation), -np.sin(rotation)],
                           [np.sin(rotation), np.cos(rotation)]])

    return wcs


def get_gaia_sources(coords: SkyCoord, obstime: Time, radius) -> Dict[str, List[float]]:
    # get the gaia time
    gaia_time = Time(GAIA_EPOCH, format='decimalyear')
    # work out the delta time between epoch and J2000.0
    with warnings.catch_warnings(record=True) as _:
        delta_time_now = (obstime.jd - gaia_time.jd) * uu.day
    # modify radius to catch all sources
    search = abs(delta_time_now * MAX_PM).to(uu.deg)
    radius = radius.to(uu.deg) + search
    # -------------------------------------------------------------------------
    # define the gaia Tap instance
    gaia = TapPlus(url=GAIA_URL)
    # launch the query
    job = gaia.launch_job(GAIA_QUERY.format(ra=coords.ra.deg,
                                            dec=coords.dec.deg,
                                            radius=radius.to(uu.deg).value))
    # get the query results
    table = job.get_results()
    # -------------------------------------------------------------------------
    # storage for sources
    sources = dict()
    sources['gaia_id'] = []
    sources['ra'] = []
    sources['dec'] = []
    sources['G'] = []
    sources['Rp'] = []
    sources['Bp'] = []
    sources['J'] = []
    sources['H'] = []
    sources['K'] = []
    sources['ra_gaia'] = []
    sources['dec_gaia'] = []
    sources['pmra'] = []
    sources['pmdec'] = []
    sources['parallax'] = []
    sources['separation'] = []
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
            plx = np.nan
        elif table_row['parallax'] <= 0:
            distance = None
            plx = np.nan
        else:
            plx = table_row['parallax']
            distance = Distance(parallax=plx * uu.mas)
        # need to propagate ra and dec to J2000
        gaia_coords = SkyCoord(ra=table_row['ra'] * uu.deg,
                               dec=table_row['dec'] * uu.deg,
                               distance=distance,
                               pm_ra_cosdec=table_row['pmra'] * uu.mas / uu.yr,
                               pm_dec=table_row['pmdec'] * uu.mas / uu.yr,
                               obstime=gaia_time, frame='icrs')
        # center the image on the current coordinates
        with warnings.catch_warnings(record=True) as _:
            curr_coords = gaia_coords.apply_space_motion(dt=delta_time_now)
        # work out the separation between this point and the center point
        separation = coords.separation(curr_coords)
        # add to sources
        sources['gaia_id'].append(table_row['source_id'])
        sources['ra'].append(curr_coords.ra.deg)
        sources['dec'].append(curr_coords.dec.deg)
        sources['G'].append(table_row['phot_g_mean_mag'])
        sources['Rp'].append(table_row['phot_rp_mean_mag'])
        sources['Bp'].append(table_row['phot_bp_mean_mag'])
        sources['J'].append(np.nan)
        sources['H'].append(np.nan)
        sources['K'].append(np.nan)
        sources['ra_gaia'].append(table_row['ra'])
        sources['dec_gaia'].append(table_row['dec'])
        sources['pmra'].append(table_row['pmra'])
        sources['pmdec'].append(table_row['pmdec'])
        sources['parallax'].append(plx)
        sources['separation'].append(separation.value)
    # -------------------------------------------------------------------------
    # convert all to numpy arrays
    sources['gaia_id'] = np.array(sources['gaia_id'])
    sources['ra'] = np.array(sources['ra'])
    sources['dec'] = np.array(sources['dec'])
    # warning is for the fact some columns are masked
    with warnings.catch_warnings(record=True) as _:
        sources['G'] = np.array(sources['G'])
        sources['Rp'] = np.array(sources['Rp'])
        sources['Bp'] = np.array(sources['Bp'])
        sources['J'] = np.array(sources['J'])
        sources['H'] = np.array(sources['H'])
        sources['K'] = np.array(sources['K'])
        sources['ra_gaia'] = np.array(sources['ra_gaia'])
        sources['dec_gaia'] = np.array(sources['dec_gaia'])
        sources['pmra'] = np.array(sources['pmra'])
        sources['pmdec'] = np.array(sources['pmdec'])
        sources['parallax'] = np.array(sources['parallax'])
        sources['separation'] = np.array(sources['separation'])
    # -------------------------------------------------------------------------
    return sources


def get_2mass_field(obs_coords, obs_time, radius):
    """
    We don't know the epoch of the 2MASS data, so we need to get the epoch
    for the field
    :param gaia_sources:
    :return:
    """
    # define the gaia Tap instance
    tmass = TapPlus(url=TMASS_URL)
    # get the 2MASS time
    tmass_time = Time(TMASS_EPOCH, format='decimalyear')
    # work out the delta time between 2mass and Gaia
    with warnings.catch_warnings(record=True) as _:
        delta_time_2mass = (tmass_time.jd - obs_time.jd) * uu.day
    # need to get the observation coords again (as a copy)
    obs_coords = SkyCoord(obs_coords)
    # get the 2mass coords at the 2mass epoch
    with warnings.catch_warnings(record=True) as _:
        tmass_coords = obs_coords.apply_space_motion(dt=delta_time_2mass)

    # modify radius to catch all sources
    search = abs(delta_time_2mass * MAX_PM).to(uu.deg)
    radius = radius.to(uu.deg) + search
    # launch the query
    job0 = tmass.launch_job(TMASS_QUERY.format(TMASS_COLS=TMASS_COLS,
                                               ra=tmass_coords.ra.deg,
                                               dec=tmass_coords.dec.deg,
                                               radius=radius.value))
    # get the query results
    table0 = job0.get_results()
    if len(table0) == 0:
        raise ValueError('No sources found in 2MASS catalog')
    else:
        for it, col in enumerate(TMASS_COLS.split(',')):
            table0[col] = table0[f'col_{it}']
            del table0[f'col_{it}']
    # get the mean date for field - we don't know when these ra and dec were
    # so we had to "guess" the date above
    mean_jdate = np.mean(table0['jdate'])
    # make this a time
    jdate_time = Time(mean_jdate, format='jd')
    return jdate_time, table0


def get_2mass_sources(gaia_sources, obs_coords, obs_time, radius):
    """
    for each source in gaia sources cross match with the 2mass/Gaia catalog

    :param gaia_sources:
    :return:
    """
    # get the mean jdate for the field
    jdate_time, tmass_table = get_2mass_field(obs_coords, obs_time, radius)
    # get all 2MASS coords
    tmass_coords = SkyCoord(ra=tmass_table['ra'], dec=tmass_table['dec'],
                            distance=None, pm_ra_cosdec=None, pm_dec=None,
                            obstime=jdate_time, frame='icrs')
    # get the gaia time
    gaia_time = Time(GAIA_EPOCH, format='decimalyear')

    # loop around objects
    for row in tqdm(range(len(gaia_sources['parallax']))):
        # deal with distances
        if np.isnan(gaia_sources['parallax'][row]):
            distance = None
        else:
            distance = Distance(parallax=gaia_sources['parallax'][row] * uu.mas)
        # need to get the gaia coords again
        gaia_coords = SkyCoord(ra=gaia_sources['ra_gaia'][row] * uu.deg,
                               dec=gaia_sources['dec_gaia'][row] * uu.deg,
                               distance=distance,
                               pm_ra_cosdec=gaia_sources['pmra'][row] * uu.mas / uu.yr,
                               pm_dec=gaia_sources['pmdec'][row] * uu.mas / uu.yr,
                               obstime=gaia_time, frame='icrs')
        # work out the delta time between 2mass and Gaia
        with warnings.catch_warnings(record=True) as _:
            delta_time_jdate = (jdate_time.jd - gaia_time.jd) * uu.day
        # get the 2mass coords at the 2mass epoch
        with warnings.catch_warnings(record=True) as _:
            jdate_coords = gaia_coords.apply_space_motion(dt=delta_time_jdate)

        # Need to find the closest sources in 2MASS to the jdate coords
        separation = jdate_coords.separation(tmass_coords)
        # find the closest source
        closest = np.argmin(separation)

        # if there is no match then we set the value to the sigma limit
        if separation[closest] > TMASS_RADIUS:
            jmag = SIGMA_LIMIT['J'] + MAG_LIMIT
            hmag = SIGMA_LIMIT['H'] + MAG_LIMIT
            kmag = SIGMA_LIMIT['K'] + MAG_LIMIT
        else:
            # get the magnitudes
            jmag = tmass_table['j_m'][closest]
            hmag = tmass_table['h_m'][closest]
            kmag = tmass_table['k_m'][closest]

        gaia_sources['J'][row] = jmag
        gaia_sources['H'][row] = hmag
        gaia_sources['K'][row] = kmag

    return gaia_sources


def seed_image(gaia_sources, pixel_scale, obs_coords, fwhm, field_of_view,
               sigma_limit, band, rotation, flip_x, flip_y, scale_factor):
    """
    Create the image WCS and seed all gaia soirces for band at their co
    """
    # plot out to the scale factor
    # TODO: fix this
    field_of_view = field_of_view * scale_factor

    # number of pixels in each direction
    npixel = int(field_of_view.to(uu.deg) // (pixel_scale * uu.pixel).to(uu.deg))

    fwhm_pix = fwhm.to(uu.arcsec) / (pixel_scale  * uu.pixel).to(uu.arcsec)

    wcs = setup_wcs((npixel, npixel), obs_coords, pixel_scale, rotation,
                    flip_x, flip_y)

    image = np.random.normal(size=(npixel, npixel), scale=1.0, loc=0)

    nsig_psf = np.array(10 ** ((sigma_limit - gaia_sources[band]) / 2.5))

    # convert all coords to pixel coords
    x_sources, y_sources = wcs.all_world2pix(gaia_sources['ra'],
                                             gaia_sources['dec'], 0)
    # create a grid of all pixels
    y, x = np.mgrid[0:npixel, 0:npixel]
    # get the width of the PSF from the fwhm
    ew = fwhm_pix.value / (2 * np.sqrt(2 * np.log(2)))
    # loop around all sources
    for i in tqdm(range(len(x_sources))):

        if np.isnan(nsig_psf[i]):
            continue

        xdiff0 = x - x_sources[i]
        ydiff0 = y - y_sources[i]

        xdiff = xdiff0 * np.cos(rotation) - ydiff0 * np.sin(rotation)
        ydiff = xdiff0 * np.sin(rotation) + ydiff0 * np.cos(rotation)
        # core of PSF
        exp1 = np.exp(-(xdiff ** 2 + ydiff ** 2) / (2 * ew ** 2))
        image += nsig_psf[i] * exp1
        # halo of PSF
        exp_halo = np.exp(-(xdiff ** 2 + ydiff ** 2) / (2 * (ew * 3) ** 2))
        image += nsig_psf[i] * exp_halo * 1e-3
        # # spike in y
        # exp_spike_y = np.exp(-np.abs(xdiff ** 2 + 0.1 * np.abs(ydiff) ** 2) / (2 * ew ** 2))
        # image += nsig_psf[i] * exp_spike_y * 1e-2
        # # spike in x
        # exp_spike_x = np.exp(-(0.1 * np.abs(xdiff) ** 2 + np.abs(ydiff) ** 2) / (2 * ew ** 2))
        # image += nsig_psf[i] * exp_spike_x * 1e-2

    # display with an arcsinh stretch
    image = np.arcsinh(image)

    # return the image and the wcs
    return image, wcs


def plot_map1(frame, image, wcs, obs_coords, title,
              field_of_view):
    # plot image
    frame.imshow(image, origin='lower', vmin=np.arcsinh(-3),
                 vmax=np.arcsinh(200),
                 cmap='gist_heat', interpolation='nearest')
    # pick a color for the annotations
    color = 'cyan'
    # -------------------------------------------------------------------------
    # plot trail
    # -------------------------------------------------------------------------
    # x_all_pix, y_all_pix = wcs.all_world2pix(all_coords.ra.value,
    #                                          all_coords.dec.value, 0)
    #
    # extent = [0, image.shape[0], 0, image.shape[1]]
    # dx = extent[1] - extent[0]
    # dy = extent[3] - extent[2]
    # # plot the year labels
    # for it, coord in enumerate(all_coords):
    #     # check whether text is out of bounds
    #     if out_of_bounds(extent, x_all_pix[it], y_all_pix[it]):
    #         continue
    #
    #     # text is offset by 1% of the ra range
    #     frame.text(float(x_all_pix[it] + 0.05 * dx),
    #                float(y_all_pix[it] + 0.05 * dy),
    #                f'{all_times[it].decimalyear:.1f}', color='yellow',
    #                alpha=0.5,
    #                horizontalalignment='left', verticalalignment='center')
    #
    #     frame.plot(x_all_pix[it], y_all_pix[it], color='y', ms=10, marker='+',
    #                ls='None', alpha=0.5)

    # -------------------------------------------------------------------------
    # plot current position
    # -------------------------------------------------------------------------
    x_curr_pix, y_curr_pix = wcs.all_world2pix(obs_coords.ra.value,
                                               obs_coords.dec.value, 0)

    frame.plot(x_curr_pix, y_curr_pix, marker='o',
               color=color, ms=30, mfc='none')
    frame.plot(x_curr_pix, y_curr_pix, marker='+',
               color=color, ms=30, mfc='none')
    # -------------------------------------------------------------------------
    # plot compass
    # -------------------------------------------------------------------------
    x_comp0 = 0.9 * image.shape[0]
    y_comp0 = 0.1 * image.shape[0]
    # convert compass center to ra/dec from pixel space
    ra_comp, dec_comp = wcs.all_pix2world(x_comp0, y_comp0, 0)
    # get the cos(dec) term
    cos_dec = np.cos(dec_comp * uu.deg)
    # get the compas
    ra_comp_e = ra_comp + (field_of_view.to(uu.deg).value/10) / cos_dec
    dec_comp_n = dec_comp + field_of_view.to(uu.deg).value/10
    # convert back to pixel space
    x_comp, y_comp = wcs.all_world2pix(ra_comp, dec_comp, 0)
    x_comp_e, y_comp_e = wcs.all_world2pix(ra_comp_e, dec_comp, 0)
    x_comp_n, y_comp_n = wcs.all_world2pix(ra_comp, dec_comp_n, 0)
    # plot the compass
    kwargs = dict(ha='center', va='center', color=color,
                  arrowprops=dict(arrowstyle='<-', color=color,
                                  shrinkA=0.0, shrinkB=0.0))
    frame.annotate('N', (x_comp, y_comp), (x_comp_n, y_comp_n), **kwargs)
    frame.annotate('E', (x_comp, y_comp), (x_comp_e, y_comp_e), **kwargs)

    # -------------------------------------------------------------------------
    # plot scale bar
    # -------------------------------------------------------------------------
    # get the length of the scale in pixels
    length = (COMPASS_SIZE / PIXEL_SCALE['G']).value

    x_scale_px = 0.1 * image.shape[0]
    y_scale_px = 0.1 * image.shape[1]
    x_stext_px = 0.1 * image.shape[0] + 0.5 * length
    y_stext_px = 0.1 * image.shape[1]

    patch = FancyArrowPatch((x_scale_px, y_scale_px),
                            (x_scale_px + length, y_scale_px),
                            arrowstyle='-', shrinkA=0.0, shrinkB=0.0,
                            capstyle='round', color=color, linewidth=2)
    frame.text(x_stext_px, y_stext_px,
               f'{COMPASS_SIZE.value:g} {COMPASS_SIZE.unit:unicode}',
               ha='center', va='bottom', color=color)
    frame.add_patch(patch)

    # -------------------------------------------------------------------------
    # remove ticks
    # -------------------------------------------------------------------------
    frame.set_xticks([])
    frame.set_yticks([])

    # -------------------------------------------------------------------------
    # add title
    # -------------------------------------------------------------------------
    frame.set_title(title, y=1.0, pad=-14, color=color)


def out_of_bounds(extent, x, y):
    xlims = extent[:2]
    ylims = extent[2:]

    lims = [xlims, ylims]
    values = [x, y]

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
    # get all Gaia data around the current coordinates
    print('Getting Gaia sources for this region')
    # TODO: Deal with >2000 sources
    gaia_sources = get_gaia_sources(obs_coords, obs_time,
                                    radius=RADIUS['G'])
    # -------------------------------------------------------------------------
    # Get the 2MASS source for each Gaia source
    print('Getting 2MASS sources for this region')
    # TODO: Deal with >2000 sources
    gaia_sources = get_2mass_sources(gaia_sources, obs_coords, obs_time,
                                     radius=RADIUS['J'])


    # -------------------------------------------------------------------------
    # seed Gaia images
    # -------------------------------------------------------------------------
    kwargs_gaia = dict(pixel_scale=PIXEL_SCALE['G'], obs_coords=obs_coords,
                        fwhm=FWHM['G'], field_of_view=FIELD_OF_VIEW['G'],
                        sigma_limit=SIGMA_LIMIT['G'], band='G',
                        scale_factor=SCALE_FACTOR)
    print('\nSeeding Gaia image')
    image1, wcs1 = seed_image(gaia_sources, flip_x=FLIP_X, flip_y=FLIP_Y,
                              rotation=TRANSFORM_ROTATE['G'].to(uu.deg),
                              **kwargs_gaia)

    print('\nSeeding Gaia image (no rotation)')
    image3, wcs3 = seed_image(gaia_sources, flip_x=False, flip_y=False,
                              rotation=0 * uu.deg, **kwargs_gaia)
    # -------------------------------------------------------------------------
    # seed 2MASS images
    # -------------------------------------------------------------------------
    kwargs_tmass = dict(pixel_scale=PIXEL_SCALE['J'], obs_coords=obs_coords,
                        fwhm=FWHM['J'], field_of_view=FIELD_OF_VIEW['J'],
                        sigma_limit=SIGMA_LIMIT['J'], band='J',
                        scale_factor=SCALE_FACTOR)
    print('\nSeeding 2MASS image')
    image2, wcs2 = seed_image(gaia_sources, flip_x=FLIP_X, flip_y=FLIP_Y,
                              rotation=TRANSFORM_ROTATE['J'].to(uu.deg),
                              **kwargs_tmass)

    print('\nSeeding 2MASS image (no rotation)')
    image4, wcs4 = seed_image(gaia_sources, flip_x=False, flip_y=False,
                              rotation=0*uu.deg, **kwargs_tmass)

    # -------------------------------------------------------------------------
    # plot figure
    # -------------------------------------------------------------------------

    fig, frames = plt.subplots(ncols=2, nrows=2, figsize=(16, 16))
    # plot the gaia map
    # noinspection PyTypeChecker
    frame11 = plot_map1(frames[0][0], image1, wcs1, obs_coords,
                        'Gaia', FIELD_OF_VIEW['G'])
    # plot the 2mass map
    # noinspection PyTypeChecker
    frame12 = plot_map1(frames[0][1], image2, wcs2, obs_coords,
                        '2MASS', FIELD_OF_VIEW['J'])

    # noinspection PyTypeChecker
    frame21 = plot_map1(frames[1][0], image3, wcs3, obs_coords,
                        'Gaia [No rotation]', FIELD_OF_VIEW['G'])
    # plot the 2mass map
    # noinspection PyTypeChecker
    frame22 = plot_map1(frames[1][1], image4, wcs4, obs_coords,
                        '2MASS [No rotation]', FIELD_OF_VIEW['J'])

    # -------------------------------------------------------------------------
    # add title
    # -------------------------------------------------------------------------
    # find the object closest to our source
    closest = np.argmin(gaia_sources['separation'])
    title = (f'Object: {objname}\n'
             f'Date: {obs_time.iso}   [{obs_time.jd}]\n'
             f'RA: {obs_coords.ra.to_string(uu.hourangle, sep=":")}   '
             f'Dec: {obs_coords.dec.to_string(uu.deg, sep=":")}   \n'
             f'Gmag: {gaia_sources["G"][closest]:.2f}   '
             f'Jmag: {gaia_sources["J"][closest]:.2f}   ')
    plt.suptitle(title)

    plt.subplots_adjust(top=0.9)

    plt.show()
    plt.close()
