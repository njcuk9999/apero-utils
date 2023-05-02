import argparse
import warnings
from typing import Dict, List

import ligo.skymap.plot as ligo_plot
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.wcs import WCS
from astroquery.utils.tap.core import TapPlus
from dateutil.parser import parse
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
FIELD_OF_VIEW['G'] = 2 * uu.arcmin
FIELD_OF_VIEW['J'] = 2 * uu.arcmin
# RADIUS
RADIUS = dict()
RADIUS['G'] = FIELD_OF_VIEW['G'] / np.sqrt(2)
RADIUS['J'] = FIELD_OF_VIEW['J'] / np.sqrt(2)
# PIXEL SCALE
PIXEL_SCALE = dict()
PIXEL_SCALE['G'] = 0.3 * uu.arcsec / uu.pixel
PIXEL_SCALE['J'] = 0.3 * uu.arcsec / uu.pixel
# FWHM
FWHM = dict()
FWHM['G'] = 2.0 * uu.arcsec
FWHM['J'] = 2.0 * uu.arcsec
# Transform the image
TRANSFORM_ROTATE = dict()
TRANSFORM_ROTATE['G'] = 22 * uu.deg
TRANSFORM_ROTATE['J'] = 22 * uu.deg
# 1-sigma magnitude limit G
SIGMA_LIMIT = dict()
SIGMA_LIMIT['G'] = 22
SIGMA_LIMIT['J'] = 17
SIGMA_LIMIT['H'] = 17
SIGMA_LIMIT['K'] = 17
# mag limit above 1-sigma limit for null sources
MAG_LIMIT = -1

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
                                            radius=RADIUS['G'].to(uu.deg).value))
    # get the query results
    table = job.get_results()
    # get the gaia time
    gaia_time = Time(GAIA_EPOCH, format='decimalyear')
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
        # work out the delta time between epoch and J2000.0
        delta_time_now = (obstime.jd - gaia_time.jd) * uu.day
        # center the image on the current coordinates
        with warnings.catch_warnings(record=True) as _:
            curr_coords = gaia_coords.apply_space_motion(dt=delta_time_now)
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
    # -------------------------------------------------------------------------
    # convert all to numpy arrays
    sources['gaia_id'] = np.array(sources['gaia_id'])
    sources['ra'] = np.array(sources['ra'])
    sources['dec'] = np.array(sources['dec'])
    sources['G'] = np.array(sources['G'])
    sources['Rp'] = np.array(sources['Rp'])
    sources['Bp'] = np.array(sources['Bp'])
    sources['J'] = np.array(sources['J'])
    sources['H'] = np.array(sources['H'])
    sources['K'] = np.array(sources['K'])
    # -------------------------------------------------------------------------
    return sources


def get_2mass_field(obs_coords, obs_time):
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
    delta_time_2mass = (tmass_time.jd - obs_time.jd) * uu.day
    # need to get the observation coords again (as a copy)
    obs_coords = SkyCoord(obs_coords)
    # get the 2mass coords at the 2mass epoch
    with warnings.catch_warnings(record=True) as _:
        tmass_coords = obs_coords.apply_space_motion(dt=delta_time_2mass)
    # launch the query
    job0 = tmass.launch_job(TMASS_QUERY.format(TMASS_COLS=TMASS_COLS,
                                               ra=tmass_coords.ra.deg,
                                               dec=tmass_coords.dec.deg,
                                               radius=RADIUS['J'].to(uu.deg).value))
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


def get_2mass_sources(gaia_sources, obs_coords, obs_time):
    """
    for each source in gaia sources cross match with the 2mass/Gaia catalog

    :param gaia_sources:
    :return:
    """
    # get the mean jdate for the field
    jdate_time, tmass_table = get_2mass_field(obs_coords, obs_time)
    # get all 2MASS coords
    tmass_coords = SkyCoord(ra=tmass_table['ra'],  dec=tmass_table['dec'],
                            distance=None, pm_ra_cosdec=None, pm_dec=None,
                            obstime=jdate_time, frame='icrs')
    # get the gaia time
    gaia_time = Time(GAIA_EPOCH, format='decimalyear')

    # TODO: This does not work - no sources matching to 2 arc sec
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


def seed_image(gaia_sources, pixel_sacle, obs_coords, fwhm, field_of_view,
               sigma_limit, band):
    # number of pixels in each direction
    npixel = int(field_of_view.to(uu.deg) // (pixel_sacle * uu.pixel).to(uu.deg))

    wcs = setup_wcs((npixel, npixel), obs_coords, pixel_sacle)

    image = np.random.normal(size=(npixel, npixel), scale=1.0, loc=0)

    nsig_psf = np.array(10 ** ((sigma_limit - gaia_sources[band]) / 2.5))

    # convert all coords to pixel coords
    x_sources, y_sources = wcs.all_world2pix(gaia_sources['ra'],
                                             gaia_sources['dec'], 0)
    # create a grid of all pixels
    y, x = np.mgrid[0:npixel, 0:npixel]
    # get the width of the PSF from the fwhm
    ew = fwhm.value / (2 * np.sqrt(2 * np.log(2)))
    # loop around all sources
    for i in tqdm(range(len(x_sources))):
        xdiff = x - x_sources[i]
        ydiff = y - y_sources[i]
        # core of PSF
        exp1 = np.exp(-(xdiff ** 2 + ydiff ** 2) / (2 * ew ** 2))
        image += nsig_psf[i] * exp1
        # halo of PSF
        exp_halo = np.exp(-(xdiff ** 2 + ydiff ** 2) / (2 * (ew * 2) ** 2))
        image += nsig_psf[i] * exp_halo * 1e-3
        # spike in y
        exp_spike_y = np.exp(-np.abs(xdiff ** 2 + 0.1 * np.abs(ydiff) ** 2) / (2 * ew ** 2))
        image += nsig_psf[i] * exp_spike_y * 1e-2
        # spike in x
        exp_spike_x = np.exp(-(0.1 * np.abs(xdiff) ** 2 + np.abs(ydiff) ** 2) / (2 * ew ** 2))
        image += nsig_psf[i] * exp_spike_x * 1e-2
    # return the image and the wcs
    return image, wcs


def plot_map(frame, image, wcs, obs_coords, all_coords, title,
             field_of_view, rotation):
    _ = ligo_plot
    # get the extent of the image
    extent = [wcs.wcs.crval[0] - wcs.wcs.cdelt[0] * image.shape[1] / 2,
              wcs.wcs.crval[0] + wcs.wcs.cdelt[0] * image.shape[1] / 2,
              wcs.wcs.crval[1] - wcs.wcs.cdelt[1] * image.shape[0] / 2,
              wcs.wcs.crval[1] + wcs.wcs.cdelt[1] * image.shape[0] / 2]
    dra = extent[1] - extent[0]
    ddec = extent[3] - extent[2]

    frame.imshow(image, origin='lower', vmin=-3, vmax=20,
                 cmap='gist_heat', extent=extent, interpolation='nearest',
                 transform=frame.get_transform('world'))

    # plot the trail
    frame.plot(all_coords.ra.deg, all_coords.dec.deg, color='yellow', ls='None',
               marker='+', lw=2, markersize=10, alpha=0.5,
               transform=frame.get_transform('world'))
    # plot the year labels
    for it, coord in enumerate(all_coords):
        # check whether text is out of bounds
        if out_of_bounds(extent, coord):
            continue
        # text is offset by 1% of the ra range
        frame.text(coord.ra.deg + 0.05 * dra, coord.dec.deg + 0.05 * ddec,
                   f'{all_times[it].decimalyear:.1f}', color='yellow',
                   alpha=0.5,
                   horizontalalignment='left', verticalalignment='center',
                   transform=frame.get_transform('world'))

    # plot the current position as an open blue circle
    frame.plot(obs_coords.ra.deg, obs_coords.dec.deg, marker='o',
               color='yellow', ms=20, mfc='none',
               transform=frame.get_transform('world'))

    # scale bar
    frame.scalebar((0.1, 0.1), field_of_view / 10, color='white').label(color='white')
    # custom compass
    compass(frame, 0.9, 0.1, 0.2, color='white', rotation=rotation)
    # plot a grid in white
    frame.grid(color='white', ls='--', lw=1, alpha=0.5)
    # set labels
    frame.set(xlabel='RA', ylabel='Dec', title=title)

    return frame


def compass(frame, x, y, size, color, rotation=0.0 * uu.deg):
    """Add a compass to indicate the north and east directions.

    Parameters
    ----------
    x, y : float
        Position of compass vertex in axes coordinates.
    size : float
        Size of compass in axes coordinates.

    """
    # get the rotations
    dx_n = size * np.sin(rotation.to(uu.rad).value)
    dy_n = size * np.cos(rotation.to(uu.rad).value)
    dx_e = size * np.sin(rotation.to(uu.rad).value - np.pi / 2)
    dy_e = size * np.cos(rotation.to(uu.rad).value - np.pi / 2)
    # get the end of the compass points
    xy = x, y
    xy_n = x + dx_n, y + dy_n
    xy_e = x + dx_e, y + dy_e
    # set up the kwargs for the annotations
    kwargs = dict(xycoords=frame.transAxes, textcoords=frame.transAxes,
                  ha='center', va='center', color=color,
                  arrowprops=dict(arrowstyle='<-', color=color,
                                  shrinkA=0.0, shrinkB=0.0))
    # add the annotations
    north = frame.annotate('N', xy, xy_n, **kwargs)
    east = frame.annotate('E', xy, xy_e, **kwargs)
    # return the annotations
    return north, east


def out_of_bounds(extent, coord):
    xlims = extent[:2]
    ylims = extent[2:]

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
    # get all Gaia data around the current coordinates
    print('Getting Gaia sources for this region')
    gaia_sources = get_gaia_sources(obs_coords, obs_time)
    # -------------------------------------------------------------------------
    # Get the 2MASS source for each Gaia source
    print('Getting 2MASS sources for this region')
    gaia_sources = get_2mass_sources(gaia_sources, obs_coords, obs_time)
    # -------------------------------------------------------------------------
    print('Seeding Gaia image')
    image1, wcs1 = seed_image(gaia_sources, PIXEL_SCALE['G'],
                              obs_coords, FWHM['G'], FIELD_OF_VIEW['G'],
                              SIGMA_LIMIT['G'], band='G')
    print('Seeding 2MASS image')
    image2, wcs2 = seed_image(gaia_sources, PIXEL_SCALE['J'],
                              obs_coords, FWHM['J'], FIELD_OF_VIEW['J'],
                              SIGMA_LIMIT['J'], band='J')

    # -------------------------------------------------------------------------
    # plot figure
    # -------------------------------------------------------------------------
    plt.close()
    fig = plt.figure()
    # setup the frames manually
    frame1 = plt.axes([0.05, 0.05, 0.4, 0.90],
                      projection='astro zoom',
                      center=obs_coords,
                      radius=FIELD_OF_VIEW['G'].to(uu.deg) / 2,
                      rotate=TRANSFORM_ROTATE['G'].to(uu.deg),
                      wcs=wcs1)
    frame2 = plt.axes([0.55, 0.05, 0.4, 0.90],
                      projection='astro zoom',
                      center=obs_coords,
                      radius=FIELD_OF_VIEW['J'].to(uu.deg) / 2,
                      rotate=TRANSFORM_ROTATE['J'].to(uu.deg),
                      wcs=wcs2)
    # plot the gaia map
    # noinspection PyTypeChecker
    frame1 = plot_map(frame1, image1, wcs1, obs_coords, all_coords, 'Gaia',
                      FIELD_OF_VIEW['G'], rotation=TRANSFORM_ROTATE['G'])
    # plot the 2mass map
    # noinspection PyTypeChecker
    frame2 = plot_map(frame2, image2, wcs2, obs_coords, all_coords, '2MASS',
                      FIELD_OF_VIEW['J'], rotation=TRANSFORM_ROTATE['J'])

    plt.show()
    plt.close()
