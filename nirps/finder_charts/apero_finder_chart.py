import argparse
from typing import Dict, List, Tuple
import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.visualization import ImageNormalize, LinearStretch
from astropy.wcs import WCS
from astroquery.utils.tap.core import TapPlus
from dateutil.parser import parse
from photutils.psf import IntegratedGaussianPRF

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
# Define the thumbnail image size
RADIUS = 5 * uu.arcmin
# TODO: Need these to simulate the image at specific band
PSF = 1.5
PIXEL_SCALE = 0.1 * uu.arcsec / uu.pixel
IMAGE_SIZE = (512, 512)
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


class Source:
    def __init__(self, gaia_id, skycoord: SkyCoord,
                 magnitudes: Dict[str, float],
                 fluxes: Dict[str, float]):
        self.gaia_id = gaia_id
        self.skycoord = skycoord
        self.magnitudes = magnitudes
        self.fluxes = fluxes


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


def get_gaia_sources(coords: SkyCoord, obstime: Time) -> List[Source]:
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
    sources = []
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

        # get magnitudes
        magnitudes = dict()
        magnitudes['G'] = table_row['phot_g_mean_mag']
        magnitudes['Rp'] = table_row['phot_rp_mean_mag']
        magnitudes['Bp'] = table_row['phot_bp_mean_mag']
        fluxes = dict()
        fluxes['G'] = table_row['phot_g_mean_flux']
        fluxes['Rp'] = table_row['phot_rp_mean_flux']
        fluxes['Bp'] = table_row['phot_bp_mean_flux']
        # add source
        source = Source(table_row['source_id'], curr_coords, magnitudes,
                        fluxes)
        # append source to sources
        sources.append(source)
    # -------------------------------------------------------------------------
    return sources


def seed_map(image_shape: Tuple[int, int], cent_coords: SkyCoord,
             sources: List[Source], band: str, sigma_psf: float,
             pixel_scale: uu.Quantity):
    # create a WCS object
    print('Settings up WCS')
    wcs = setup_wcs(image_shape, cent_coords, pixel_scale)
    # find the x and y pixel coordinates of the source
    x_cent, y_cent = wcs.all_world2pix(cent_coords.ra.deg,
                                       cent_coords.dec.deg, 0)

    # TODO: Can and should be replaced with the actual PSF model
    #       psfmodel = to_griddedpsfmodel(hdu, ext=0)
    psfmodel = IntegratedGaussianPRF(sigma=sigma_psf, x_0=x_cent, y_0=y_cent)

    image = np.zeros(image_shape)
    # create image
    print('\tAdding sources to image...')
    # loop around sources
    for it, source in enumerate(sources):
        print('\t\tAdding source {0} of {1}'.format(it + 1, len(sources)))
        flux = source.fluxes[band]
        # find the x and y pixel coordinates of the source
        x_source, y_source = wcs.all_world2pix(source.skycoord.ra.deg,
                                               source.skycoord.dec.deg, 0)

        y, x = np.mgrid[0:image_shape[0], 0:image_shape[1]]

        image[y, x] += psfmodel.evaluate(x=x, y=y, flux=flux,
                                         x_0=x_source, y_0=y_source,
                                         sigma=sigma_psf)
    # return the image
    return image, wcs


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
    # -------------------------------------------------------------------------

    print('Seeding image')

    # seed sources onto map
    seed_image, wcs = seed_map(IMAGE_SIZE, obs_coords, gaia_sources, 'G',
                               PSF, PIXEL_SCALE)

    # -------------------------------------------------------------------------
    # plot figure
    # -------------------------------------------------------------------------
    # get the extent of the image
    extent = [wcs.wcs.crval[0] - wcs.wcs.cdelt[0] * seed_image.shape[1] / 2,
              wcs.wcs.crval[0] + wcs.wcs.cdelt[0] * seed_image.shape[1] / 2,
              wcs.wcs.crval[1] - wcs.wcs.cdelt[1] * seed_image.shape[0] / 2,
              wcs.wcs.crval[1] + wcs.wcs.cdelt[1] * seed_image.shape[0] / 2]
    dra = extent[1] - extent[0]
    ddec = extent[3] - extent[2]

    # remove extreme values
    low_lim, high_lim = np.percentile(seed_image, [0.5, 99.5])
    seed_image[seed_image < low_lim] = low_lim
    seed_image[seed_image > high_lim] = high_lim

    plt.close()
    fig = plt.figure()
    frame = plt.subplot(111, projection=wcs)

    frame.imshow(seed_image, origin='lower', extent=extent, cmap='gist_heat')
    frame.set(xlabel='RA', ylabel='DEC', xlim=extent[:2], ylim=extent[2:],
              title=f'{objname} ({obs_time.datetime.strftime("%Y-%m-%d")})')
    # plot the trail
    frame.plot(all_coords.ra, all_coords.dec, color='yellow', ls='-',
               marker='+', lw=2, markersize=10)
    # plot the year labels
    for it, coord in enumerate(all_coords):
        # check whether text is out of bounds
        if out_of_bounds(frame, coord):
            continue
        # text is offset by 1% of the ra range
        frame.text(coord.ra.value + 0.01 * dra, coord.dec.value,
                   f'{all_times[it].decimalyear:.1f}', color='yellow')

    # plot the current position as an open blue circle
    frame.plot(obs_coords.ra, obs_coords.dec, marker='o',
               color='yellow', ms=20, mfc='none')


    plt.show()


