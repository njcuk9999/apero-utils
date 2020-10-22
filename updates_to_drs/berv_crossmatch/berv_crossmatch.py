#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CODE DESCRIPTION HERE

Created on 2020-10-2020-10-22 11:49

@author: cook
"""
from astropy.table import Table
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
import numpy as np
from typing import Union
import warnings

from apero.base import base
from apero.base import drs_text
from apero.core import constants
from apero.core.core import drs_database
from apero.core.core import drs_log
from apero import lang

# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'science.extract.crossmatch.py'
__INSTRUMENT__ = 'None'
__PACKAGE__ = base.__PACKAGE__
__version__ = base.__version__
__author__ = base.__author__
__date__ = base.__date__
__release__ = base.__release__
# Get Astropy Time and Time Delta
Time, TimeDelta = base.AstropyTime, base.AstropyTimeDelta
# get param dict
ParamDict = constants.ParamDict
# Get Logging function
WLOG = drs_log.wlog
# Get the text types
TextEntry = lang.core.drs_lang_text.TextEntry
TextDict = lang.core.drs_lang_text.TextDict



OBJECT_FILE = 'unique_objs.csv'
CMRADIUS = 3 * 60

QUERY_GAIA = 'SELECT {0} FROM {1} WHERE {2}'
QCOLS = ('ra as ra, dec as dec, source_id as gaiaid, parallax as plx, '
         'pmdec as pmde, pmra as pmra, radial_velocity as rv, '
         'phot_g_mean_mag as gmag, phot_bp_mean_mag as bpmag, '
         'phot_rp_mean_mag as rpmag')
QSOURCE = 'gaiadr2.gaia_source'

QCIRCLE = ('(1=CONTAINS(POINT(\'ICRS\', ra, dec), CIRCLE(\'ICRS\', {ra}, '
           '{dec}, {radius})))')

# get object database
ObjectDatabase = drs_database.ObjectDatabase

# =============================================================================
# Define functions
# =============================================================================



class AstroObject(object):
    objname: Union[str, None]
    gaia_id: Union[str, None]
    ra: Union[float, None]
    dec: Union[float, None]
    pmra: Union[float, None]
    pmde: Union[float, None]
    plx: Union[float, None]
    rv: Union[float, None]
    gmag: Union[float, None]
    bpmag: Union[float, None]
    rpmag: Union[float, None]
    teff: Union[float, None]
    aliases: Union[str, None]
    used: int

    def __init__(self, params, pconst, gaia_id: Union[str, None],
                 ra: Union[str, float], dec: Union[str, float],
                 database: ObjectDatabase,
                 objname: Union[str, None], pmra: Union[float, None],
                 pmde: Union[float, None], plx: Union[float, None],
                 rv: Union[float, None]):
        """
        Construct an astrophysical object instance

        :param gaia_id: str or None, input gaia id - if None require ra and dec
        :param ra: float, right ascension, only required if gaia_id is None
        :param dec: float, declination, only required if gaia_id is None
        :param database:
        """
        # properties from input
        self.input_gaiaid = gaia_id
        self.input_ra = ra
        self.input_dec = dec
        self.database = database
        if objname is None:
            self.input_objname = None
        else:
            self.input_objname = pconst.DRS_OBJ_NAME(objname)
        self.input_pmra = pmra
        self.input_pmde = pmde
        self.input_plx = plx
        self.input_rv = rv
        # other information
        self.pconst = pconst
        self.url = params['OBJ_LIST_GAIA_URL']
        self.radius = params['OBJ_LIST_CROSS_MATCH_RADIUS']
        self.maglimit = params['OBJ_LIST_GAIA_MAG_CUT']
        # properties we need to find
        self.objname = None
        self.gaia_id = None
        self.ra = None
        self.dec = None
        self.pmra = None
        self.pmde = None
        self.plx = None
        self.rv = None
        self.gmag = None
        self.bpmag = None
        self.rpmag = None
        self.epoch = None
        self.teff = None
        self.aliases = None
        self.used = 0

    def resolve_target(self, mjd=None):
        # deal with database not loaded
        self.database.load_db()
        # ---------------------------------------------------------------------
        # 1. try gaia id and objname against database
        # ---------------------------------------------------------------------
        # self._resolve_from_database()
        # ---------------------------------------------------------------------
        # 2. try gaia id against gaia query (only if gaia_id is still None)
        # ---------------------------------------------------------------------
        #if self.gaia_id is None:
        #    self._resolve_from_gaia_id()
        # ---------------------------------------------------------------------
        # 3. use ra + dec to get gaia id (only if gaia_id is still None)
        # ---------------------------------------------------------------------
        if self.gaia_id is None:
            self._resolve_from_coords(mjd)
        # ---------------------------------------------------------------------
        # 4. if we still cannot find gaia id use input + default parameters
        # ---------------------------------------------------------------------
        if self.gaia_id is None:
            self._use_inputs()

    def _resolve_from_database(self):
        """
        Use gaia id (from input) to check if this object is already in database
        - if it is then update all parameters from here
        - if it isn't (or isn't set) check against known names in the database
        - if names aren't found and gaia id not found do not update parameters
        """
        # ---------------------------------------------------------------------
        # deal with no input_gaiaid
        if self.input_gaiaid is not None:
            # condition condition
            condition = 'GAIAID=="{0}"'.format(self.input_gaiaid)
            # get the entries for this condition
            entries = self.database.get_entries('*', condition=condition)
        else:
            entries = None
        # ---------------------------------------------------------------------
        # deal with no entries (try resolving from name in known aliases)
        if entries is None or len(entries) == 0:
            entries = self._resolve_from_names()
        if entries is None or len(entries) == 0:
            return
        # ---------------------------------------------------------------------
        # fill out required information if available
        self.objname = entries['OBJNAME']
        self.gaia_id = entries['GAIAID']
        self.ra = entries['RA']
        self.dec = entries['DEC']
        # assign pmra
        if not drs_text.null_text(entries['PMRA'], ['None']):
            self.pmra = entries['PMRA']
        # assign pmde
        if not drs_text.null_text(entries['PMDE'], ['None']):
            self.pmde = entries['PMDE']
        # assign pmde
        if not drs_text.null_text(entries['PLX'], ['None']):
            self.plx = entries['PLX']
        # assign rv
        if not drs_text.null_text(entries['RV'], ['None']):
            self.rv = entries['RV']
        # assign gmag
        if not drs_text.null_text(entries['GMAG'], ['None']):
            self.gmag = entries['GMAG']
        # assign bpmag
        if not drs_text.null_text(entries['BPMAG'], ['None']):
            self.bpmag = entries['BPMAG']
        # assign rpmag
        if not drs_text.null_text(entries['RPMAG'], ['None']):
            self.rpmag = entries['RPMAG']
        # assign epoch
        if not drs_text.null_text(entries['EPOCH'], ['None']):
            self.epoch = entries['EPOCH']
        # assign teff
        if not drs_text.null_text(entries['TEFF'], ['None']):
            self.teff = entries['TEFF']
        # assign aliases
        if not drs_text.null_text(entries['ALIASES'], ['None']):
            self.teff = entries['ALIASES'].split(',')
        # set used
        self.used = 1

    def _resolve_from_names(self):
        """
        Search the database for a named column

        :return:
        """
        # deal with no object name (shouldn't be possible)
        if self.input_objname is None:
            return None
        # get aliases from database
        gaia_table = self.database.get_entries('GAIAID, OBJNAME, ALIASES')
        # extract required columns
        gaia_id = gaia_table['GAIAID']
        objnames = gaia_table['OBJNAME']
        alias_sets = gaia_table['ALIASES']

        # ---------------------------------------------------------------------
        # 1. check direct object name
        # ---------------------------------------------------------------------
        for row, objname in enumerate(objnames):
            # get cleaned alias
            cobjname = pconst.DRS_OBJ_NAME(objname)
            # compare to input_objname
            if cobjname == self.input_objname:
                # condition condition
                condition = 'GAIAID=="{0}"'.format(gaia_id[row])
                # return the entries for this gaia id
                return self.database.get_entries('*', condition=condition)

        # ---------------------------------------------------------------------
        # 2. check aliases
        # ---------------------------------------------------------------------
        # loop around each set of aliases and see if
        for row, alias_set in enumerate(alias_sets):
            # split the names by a comma
            aliases = alias_set.split(',')
            # loop around aliases - if alias found return gaia id for this
            for alias in aliases:
                # ignore None
                if drs_text.null_text(alias, ['None']):
                    continue
                # get cleaned alias
                calias = pconst.DRS_OBJ_NAME(alias)
                # compare to input_objname
                if calias == self.input_objname:
                    # condition condition
                    condition = 'GAIAID=="{0}"'.format(gaia_id[row])
                    # return the entries for this gaia id
                    return self.database.get_entries('*', condition=condition)
        # if we have reached this point we cannot match to name
        #   therefore return None
        return None

    def _resolve_from_gaia_id(self):
        """
        Use input_gaiaid to query gaia and them update all parameters based
        on this id

        :return:
        """
        # deal with no input_gaiaid
        if self.input_gaiaid is None:
            return
        # set up ID query
        condition = 'source_id = {0}'.format(self.input_gaiaid)
        # construct sql query
        query = QUERY_GAIA.format(QCOLS, QSOURCE, condition)
        # return results
        entries = query_gaia(self.url, query)
        # deal with no entries
        if entries is None:
            return
        if len(entries) == 0:
            return
        # fill out required information if available
        self.gaia_id = entries['gaiaid'][0]
        self.ra = entries['ra'][0]
        self.dec = entries['dec'][0]
        # assign pmra
        if not entries['pmra'].mask[0]:
            self.pmra = entries['pmra'][0]
        # assign pmde
        if not entries['pmde'].mask[0]:
            self.pmra = entries['pmde'][0]
        # assign plx
        if not entries['plx'].mask[0]:
            self.plx = entries['plx'][0]
        # assign rv
        if not entries['rv'].mask[0]:
            self.plx = entries['rv'][0]
        # assign gmag
        if not entries['gmag'].mask[0]:
            self.plx = entries['gmag'][0]
        # assign bpmag
        if not entries['bpmag'].mask[0]:
            self.plx = entries['bpmag'][0]
        # assign rpmag
        if not entries['rpmag'].mask[0]:
            self.plx = entries['rpmag'][0]
        # assign epoch
        self.epoch = 2015.5
        # set used
        self.used = 1

    def _resolve_from_coords(self, mjd=None):
        # deal with ra and dec (want it in degrees)
        ra, dec = parse_coords(self.input_ra, self.input_dec)
        # get radius in degrees
        radius = (self.radius * uu.arcsec).to(uu.deg).value
        # set up ra / dec crossmatch
        condition = QCIRCLE.format(ra=ra, dec=dec, radius=radius)
        # add additional criteria (to narrow search)
        condition += r' AND (phot_rp_mean_mag < {0})'.format(self.maglimit)
        # construct sql query
        query = QUERY_GAIA.format(QCOLS, QSOURCE, condition)
        # return results of query
        entries = query_gaia(self.url, query)
        # deal with no entries
        if entries is None:
            return
        if len(entries) == 0:
            return
        # get closest to ra and dec (propagated)
        position = best_gaia_entry(ra, dec, mjd, entries)
        # fill out required information if available
        self.gaia_id = entries['gaiaid']
        self.ra = entries['ra'][position]
        self.dec = entries['dec'][position]
        # assign pmra
        if not entries['pmra'].mask[position]:
            self.pmra = entries['pmra'][position]
        # assign pmde
        if not entries['pmde'].mask[position]:
            self.pmra = entries['pmde'][position]
        # assign plx
        if not entries['plx'].mask[position]:
            self.plx = entries['plx'][position]
        # assign rv
        if not entries['rv'].mask[position]:
            self.plx = entries['rv'][position]
        # assign gmag
        if not entries['gmag'].mask[position]:
            self.plx = entries['gmag'][position]
        # assign bpmag
        if not entries['bpmag'].mask[position]:
            self.plx = entries['bpmag'][position]
        # assign rpmag
        if not entries['rpmag'].mask[position]:
            self.plx = entries['rpmag'][position]
        # assign epoch
        self.epoch = 2015.5
        # set used
        self.used = 1

    def _use_inputs(self):
        pass



def query_gaia(url, query) -> Union[Table, None]:
    """
    Query Gaia via a TapPlus query

    :param url: str, the URL to the SQL database
    :param query: str, the SQL query to run

    :return: astropy.table.Table or None - the results of the gaia TAP query
    """
    # set fucntion name
    func_name = __NAME__ + '.query_gaia()'
    # check for astroquery and return a fail and warning if not installed
    try:
        from astroquery.utils.tap.core import TapPlus

    except Exception as e:
        eargs = [type(e), str(e), func_name]
        WLOG(params, 'warning', TextEntry('10-016-00009', args=eargs))
        return None
    # ------------------------------------------------------------------
    # try running gaia query
    try:
        with warnings.catch_warnings(record=True) as _:
            # construct gaia TapPlus instance
            gaia = TapPlus(url=url)
            # launch gaia job
            job = gaia.launch_job(query=query)
            # get gaia table
            table = job.get_results()
    except Exception as e:
        wargs = [url, query, type(e), e, func_name]
        WLOG(params, 'warning', TextEntry('10-016-00008', args=wargs))
        # return No row and True to fail
        return None
    # ------------------------------------------------------------------
    # if we have no entries we did not find object
    if len(table) == 0:
        # return None
        return None
    # else we return result
    return table


def parse_coords(ra: float, dec: float, ra_unit='deg', dec_unit='deg'):
    """
    Convert coordinates into degrees via SkyCoord

    :param ra: right ascension (with units "ra_units")
    :param dec: declination (with units "dec_units"
    :param ra_unit: units for right ascension
    :param dec_unit: units for declination
    :return:
    """
    # get Sky Coord instances
    coord = SkyCoord(ra, dec, unit=[ra_unit, dec_unit])

    return coord.ra.value, coord.dec.value



def best_gaia_entry(ra: float, dec: float, mjd: float, entries: Table):
    """
    Using the originally supplied ra and dec choose the closest
    entries (propagating all entries in time to match current ra and dec
    at time = 'mjd')

    :param ra: float, the right ascension in degrees
    :param dec: float, the declination in degrees
    :param mjd: float, the modified julien date for observation
    :param entries:
    :return:
    """

    print(ra, dec)

    # get gaia time and observation time
    gaia_time = Time('2015.5', format='decimalyear')
    obs_time = Time(mjd, format='mjd')

    ra_arr = np.array(entries['ra']) * uu.deg
    dec_arr = np.array(entries['dec']) * uu.deg
    pmra_arr = np.array(entries['pmra']) * uu.mas/uu.yr
    pmde_arr = np.array(entries['pmde']) * uu.mas/uu.yr
    plx_arr = np.array(entries['plx']) * uu.mas


    # propagate all entries ra and dec to mjd
    coords0 = SkyCoord(ra_arr, dec_arr,
                       pm_ra_cosdec=pmra_arr, pm_dec=pmde_arr,
                       distance=Distance(parallax=plx_arr),
                       obstime=gaia_time)
    # apply space motion
    coords1 = coords0.apply_space_motion(obs_time)

    # crossmatch with ra and dec and keep closest

    # TODO: Got to here



# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # load table
    table = Table.read(OBJECT_FILE, format='csv')

    # load params
    params = constants.load('SPIROU')
    pconst = constants.pload('SPIROU')
    # for tests lets set the cross match radius
    params.set('OBJ_LIST_CROSS_MATCH_RADIUS', CMRADIUS)
    # load database
    objdbm = ObjectDatabase(params)
    objdbm.load_db()

    # loop around objects (later)
    row = 0
    # get properties
    gaia_id = '4472832130942575872'
    ra = table['RA_DEG'][row]
    dec = table['DEC_DEG'][row]
    objname = table['OBJECT'][row]
    # mjdmid
    exptime_days = (table['EXPTIME'][row] / 2.0) / (3600 * 24.0)
    mjdmid = table['MJDEND'][row] - exptime_days

    # set up an object instance for this target
    astro_obj = AstroObject(params, pconst, gaia_id, ra, dec, objdbm, objname,
                            0.0, 0.0, 0.0, 0.0)

    # resolve target (from gaia id or ra/dec)
    astro_obj.resolve_target(mjdmid)






# =============================================================================
# End of code
# =============================================================================
