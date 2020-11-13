#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CODE DESCRIPTION HERE

Created on 2020-10-2020-10-22 11:49

@author: cook
"""
from astropy.time import Time
from astropy.table import Table
from astropy import units as uu
from astropy.coordinates import SkyCoord, Distance
import numpy as np
import requests
from typing import List, Union
import warnings


# =============================================================================
# Define variables
# =============================================================================
__NAME__ = 'berv_crossmatch.py'

# Fake Logging function
def WLOG(params, level, message, *args, **kwargs):
    print(message)
    if level == 'error':
        raise Exception(message)


# Fake TextEntry
def TextEntry(code, args, **kwargs):
    return code


# Fake drs_text
class drstext:
    def null_text(self, value, nulls):
        if value is None:
            return True
        if value in nulls:
            return True
        return False


drs_text = drstext()


class Pconst:
    def DRS_OBJ_NAME(self, objname: str) -> str:
        """
        Clean and standardize an object name

        Default action: make upper case and remove white spaces

        :param objname: str, input object name
        :return:
        """
        # set function name
        _ = __NAME__ + '.DRS_OBJ_NAME()'
        # clean object name
        rawobjname = str(objname)
        objectname = rawobjname.strip()
        objectname = objectname.replace(' ', '_')
        objectname = objectname.upper()
        # return object name
        return objectname

pconst = Pconst()

# Fake params
params = dict()
params['OBJ_LIST_GAIA_URL'] = 'https://gea.esac.esa.int/tap-server/tap'
params['OBJ_LIST_CROSS_MATCH_RADIUS'] = 180
params['OBJ_LIST_GAIA_MAG_CUT'] = 15.0
params['OBJ_LIST_GAIA_PLX_LIM'] = 0.5
params['OBJ_LIST_GOOGLE_SHEET_URL'] = ('1jwlux8AJjBMMVrbg6LszJIpFJ'
                                       'rk6alhbT5HA7BiAHD8')
params['OBJ_LIST_GOOGLE_SHEET_WNUM'] = 0

USE_DATABASE = False

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

# cache for google sheet
GOOGLE_TABLES = dict()
# define standard google base url
GOOGLE_BASE_URL = ('https://docs.google.com/spreadsheets/d/{}/gviz/'
                   'tq?tqx=out:csv&sheet={}')


# =============================================================================
# Define functions
# =============================================================================
class AstroObject(object):
    aliases: List[str]
    used: int

    def __init__(self, params, pconst, gaia_id: Union[str, None],
                 ra: Union[str, float], dec: Union[str, float],
                 database, objname: Union[str, None], pmra: float = np.nan,
                 pmde: float = np.nan, plx: float = np.nan,
                 rv: float = np.nan, teff: float = np.nan):
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
        self.input_teff = teff
        # other information
        self.pconst = pconst
        self.gaia_url = params['OBJ_LIST_GAIA_URL']
        self.radius = params['OBJ_LIST_CROSS_MATCH_RADIUS']
        self.maglimit = params['OBJ_LIST_GAIA_MAG_CUT']
        self.plxlimit = params['OBJ_LIST_GAIA_PLX_LIM']
        self.gsheet_url = params['OBJ_LIST_GOOGLE_SHEET_URL']
        self.gsheet_wnum = params['OBJ_LIST_GOOGLE_SHEET_WNUM']
        # properties we need to find
        self.objname = None
        self.gaia_id = None
        self.ra = np.nan
        self.dec = np.nan
        self.pmra = np.nan
        self.pmde = np.nan
        self.plx = np.nan
        self.rv = np.nan
        self.gmag = np.nan
        self.bpmag = np.nan
        self.rpmag = np.nan
        self.epoch = np.nan
        self.teff = np.nan
        self.aliases = []
        self.source = None
        self.used = 0

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if self.input_objname is not None:
            return 'AstroObject[{0}]'.format(self.input_objname)
        elif self.input_gaiaid is not None:
            return 'AstroObject[Gaia DR2 {0}]'.format(self.input_gaiaid)
        else:
            return 'AstroObject[ra={0},dec={1}]'.format(self.input_ra,
                                                        self.input_dec)

    def resolve_target(self, mjd=None):
        # deal with database not loaded
        if USE_DATABASE:
            self.database.load_db()
        # ---------------------------------------------------------------------
        # 1. try gaia id and objname against database
        # ---------------------------------------------------------------------
        self._resolve_from_database()
        # ---------------------------------------------------------------------
        # 2. try gaia id against gaia query (only if gaia_id is still None)
        # ---------------------------------------------------------------------
        if self.gaia_id is None:
           self._resolve_from_gaia_id()
        # ---------------------------------------------------------------------
        # 3. try to get gaia id from google sheet of gaia id (with object name)
        # ---------------------------------------------------------------------
        if self.gaia_id is None:
            self._resolve_from_glist()
        # ---------------------------------------------------------------------
        # 4. use ra + dec to get gaia id (only if gaia_id is still None)
        # ---------------------------------------------------------------------
        # if self.gaia_id is None:
        #     self._resolve_from_coords(mjd)
        # ---------------------------------------------------------------------
        # 5. if we still cannot find gaia id use input + default parameters
        # ---------------------------------------------------------------------
        if self.gaia_id is None:
            self._use_inputs()

    def get_simbad_aliases(self):
        # deal with aliases already set
        if len(self.aliases) > 0:
            return
        # storage aliases
        self.aliases = []
        # set aliases to objname and input objname if different
        if self.objname is not None:
            self.aliases += [self.objname]
        if self.input_objname is not None:
            if self.objname != self.input_objname:
                self.aliases += [self.input_objname]
        # only do this is we have a gaia-id
        if self.gaia_id is None:
            return
        obj_id = 'Gaia DR2 {0}'.format(self.gaia_id)
        # get entries for this gaia id
        entries = query_simbad_id(obj_id)
        # deal with no entries
        if entries is None:
            return
        if len(entries) == 0:
            return
        # get the aliases
        raw_aliases = entries['IDS'][0].decode('utf-8')
        # slit via the pipe (|)
        self.aliases += raw_aliases.split('|')

    def _resolve_from_database(self):
        """
        Use gaia id (from input) to check if this object is already in database
        - if it is then update all parameters from here
        - if it isn't (or isn't set) check against known names in the database
        - if names aren't found and gaia id not found do not update parameters
        """
        if not USE_DATABASE:
            return
        # ---------------------------------------------------------------------
        # deal with no input_gaiaid
        if self.input_gaiaid is not None:
            # condition condition
            condition = 'GAIAID=="{0}"'.format(self.input_gaiaid)
            # get the entries for this condition
            entries = self.database.get_entries('*', condition=condition)
            # set source
            self.source = 'database-gaia-id'
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
            self.aliases = entries['ALIASES'].split('|')
        # set used
        self.used = 1

    def _resolve_from_names(self):
        """
        Search the database for a named column

        :return:
        """
        if not USE_DATABASE:
            return
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
                # set source
                self.source = 'database-objname'
                # return the entries for this gaia id
                return self.database.get_entries('*', condition=condition)

        # ---------------------------------------------------------------------
        # 2. check aliases
        # ---------------------------------------------------------------------
        # loop around each set of aliases and see if
        for row, alias_set in enumerate(alias_sets):
            # split the names by a comma
            aliases = alias_set.split('|')
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
                    # set source
                    self.source = 'database-aliases'
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
        entries = query_gaia(self.gaia_url, query)
        # deal with no entries
        if entries is None:
            return
        if len(entries) == 0:
            return
        # fill out required information if available
        self.objname = str(self.input_objname)
        self.gaia_id = str(entries['gaiaid'][0])
        self.ra = float(entries['ra'][0])
        self.dec = float(entries['dec'][0])
        # assign pmra
        if not entries['pmra'].mask[0]:
            self.pmra = float(entries['pmra'][0])
        # assign pmde
        if not entries['pmde'].mask[0]:
            self.pmde = float(entries['pmde'][0])
        # assign plx
        if not entries['plx'].mask[0]:
            self.plx = float(entries['plx'][0])
        # assign rv
        if not entries['rv'].mask[0]:
            self.rv = float(entries['rv'][0])
        # assign gmag
        if not entries['gmag'].mask[0]:
            self.gmag = float(entries['gmag'][0])
        # assign bpmag
        if not entries['bpmag'].mask[0]:
            self.bpmag = float(entries['bpmag'][0])
        # assign rpmag
        if not entries['rpmag'].mask[0]:
            self.rpmag = float(entries['rpmag'][0])
        # assign epoch
        self.epoch = 2015.5
        # set used
        self.used = 1
        # set source
        self.source = 'gaia-query-id-input'

    def _resolve_from_glist(self):
        """
        Resolve gaia id from google sheets (using object name) and then
        retry the gaia id against the gaia database

        :return:
        """
        # try to get gaia id from google sheets
        self.input_gaiaid = query_glist(self.input_objname, self.gsheet_url,
                                        self.gsheet_wnum)
        #  try (again) gaia id against gaia query (only if gaia_id is still
        #  None)
        if self.input_gaiaid is not None:
           self._resolve_from_gaia_id()
        # set source
        self.source = 'gaia-query-id-gsheet'

    def _resolve_from_coords(self, mjd=None):
        """
        Resolve from Gaia using coordinates (and the current date)

        :param mjd: observation modified julien date

        :return:
        """
        # deal with ra and dec (want it in degrees)
        ra, dec = parse_coords(self.input_ra, self.input_dec)
        # get radius in degrees
        radius = (self.radius * uu.arcsec).to(uu.deg).value
        # set up ra / dec crossmatch
        condition = QCIRCLE.format(ra=ra, dec=dec, radius=radius)
        # add additional criteria (to narrow search)
        condition += r' AND (phot_rp_mean_mag < {0})'.format(self.maglimit)
        # add a parallax condition
        condition += r' AND (parallax > {0})'.format(self.plxlimit)
        # construct sql query
        query = QUERY_GAIA.format(QCOLS, QSOURCE, condition)
        # return results of query
        entries = query_gaia(self.gaia_url, query)
        # deal with no entries
        if entries is None:
            return
        if len(entries) == 0:
            return
        # get closest to ra and dec (propagated)
        position = best_gaia_entry(ra, dec, mjd, entries)
        # fill out required information if available
        self.objname = str(self.input_objname)
        self.gaia_id = str(entries['gaiaid'][position])
        self.ra = float(entries['ra'][position])
        self.dec = float(entries['dec'][position])
        # assign pmra
        if not entries['pmra'].mask[position]:
            self.pmra = float(entries['pmra'][position])
        # assign pmde
        if not entries['pmde'].mask[position]:
            self.pmde = float(entries['pmde'][position])
        # assign plx
        if not entries['plx'].mask[position]:
            self.plx = float(entries['plx'][position])
        # assign rv
        if not entries['rv'].mask[position]:
            self.rv = float(entries['rv'][position])
        # assign gmag
        if not entries['gmag'].mask[position]:
            self.gmag = float(entries['gmag'][position])
        # assign bpmag
        if not entries['bpmag'].mask[position]:
            self.bpmag = float(entries['bpmag'][position])
        # assign rpmag
        if not entries['rpmag'].mask[position]:
            self.rpmag = float(entries['rpmag'][position])
        # assign epoch
        self.epoch = 2015.5
        # set used
        self.used = 1
        # set source
        self.source = 'gaia-query-coords'

    def _use_inputs(self):
        """
        If all else fails use the input values

        :return:
        """
        # fill out required information if available
        self.gaia_id = None
        self.objname = self.input_objname
        self.ra = self.input_ra
        self.dec = self.input_dec
        self.pmra = self.input_pmra
        self.pmde = self.input_pmde
        self.plx = self.input_plx
        self.rv = self.input_rv
        self.gmag = np.nan
        self.bpmag = np.nan
        self.rpmag = np.nan
        self.epoch = np.nan
        self.teff = np.nan
        self.aliases = []
        self.source = None
        self.used = 1
        # set source
        self.source = 'input'

    def write_obj(self, database, commit: bool = True):
        # write to database
        database.add_entry(objname=self.objname, gaia_id=self.gaia_id,
                           ra=self.ra, dec=self.dec,
                           pmra=self.pmra, pmde=self.pmde,
                           plx=self.plx, rv=self.rv,
                           gmag=self.gmag, bpmag=self.bpmag,
                           rpmag=self.rpmag, epoch=self.epoch,
                           teff=self.input_teff, aliases=self.aliases,
                           commit=commit)

    def write_table(self, outdict):
        """
        Proxy write function (used to write to dictionary --> Table)
        :param outdict:
        :return:
        """

        columns = ['OBJNAME', 'GAIAID', 'RA', 'DEC', 'PMRA', 'PMDE', 'PLX',
                   'RV', 'GMAG', 'BPMAG', 'RPMAG', 'EPOCH', 'TEFF', 'ALIASES',
                   'USED', 'SOURCE']
        values = [self.objname, self.gaia_id, self.ra, self.dec, self.pmra,
                  self.pmde, self.plx, self.rv, self.gmag,
                  self.bpmag, self.rpmag, self.epoch, self.teff]
        # deal with aliases
        if isinstance(self.aliases, str):
            values.append(self.aliases)
        elif isinstance(self.aliases, list):
            values.append('|'.join(self.aliases))
        else:
            values.append('None')
        # add used
        values.append(self.used)
        # add source
        values.append(self.source)
        # loop around and add to outdict
        for row in range(len(values)):
            if columns[row] in outdict:
                outdict[columns[row]].append(values[row])
            else:
                outdict[columns[row]] = [values[row]]
        # return outdict
        return outdict


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


def query_simbad_id(obj_id: str) -> Union[Table, None]:
    # set fucntion name
    _ = __NAME__ + '.query_simbad()'
    # check for astroquery and return a fail and warning if not installed
    try:
        # import astroquery
        from astroquery.simbad import Simbad
        # get results
        with warnings.catch_warnings(record=True) as _:
            # add ids column
            Simbad.add_votable_fields('ids')
            # query simbad
            return Simbad.query_object(obj_id)
    # deal with all exceptions here
    except Exception as e:
        # log that there was an error with astroquery
        wargs = [obj_id, type(e), str(e)]
        WLOG(params, 'warning', TextEntry('10-016-00020', args=wargs))
        # return unset ra/dec
        return None


def query_glist(objname: str, sheet_id: str, worksheet: int = 0):

    # get the google sheet
    gtable = get_google_sheet(sheet_id, worksheet)

    # deal with empty table
    if gtable is None:
        return None
    if len(gtable) == 0:
        return None
    # set initial position to None
    position = None
    row = np.nan
    # loop around rows and look for aliases
    for row in range(len(gtable)):
        # set aliases as the objname
        aliases = [gtable['OBJNAME'][row]]
        # get the aliases for this row
        aliases += gtable['ALIASES'][row].split('|')
        # search for object name
        position = crossmatch_name(objname, aliases)
        # break if we have found a match
        if position is not None:
            break
    # if position is still None return None
    if position is None:
        return None
    # else we have our Gaia id so return it
    return gtable['GAIAID'][row]


def get_google_sheet(sheet_id: str, worksheet: int = 0,
                     cached: bool = True) -> Table:
    """
    Load a google sheet from url using a sheet id (if cached = True and
    previous loaded - just loads from memory)

    :param sheet_id: str, the google sheet id
    :param worksheet: int, the worksheet id (defaults to 0)
    :param cached: bool, if True and previous loaded, loads from memory

    :return: Table, astropy table representation of google sheet
    """
    # set google cache table as global
    global GOOGLE_TABLES
    # construct url for worksheet
    url = GOOGLE_BASE_URL.format(sheet_id, worksheet)
    # deal with table existing
    if url in GOOGLE_TABLES and cached:
        return GOOGLE_TABLES[url]
    # get data using a request
    rawdata = requests.get(url)
    # convert rawdata input table
    table = Table.read(rawdata.text, format='ascii')
    # add to cached storage
    GOOGLE_TABLES[url] = table
    # return table
    return table


def crossmatch_name(name: str, namelist: List[str]) -> Union[int, None]:
    """
    Crossmatch a name with a list of names (returning position in the list of
    names)

    :param name: str, the name to search for
    :param namelist: list of strings, the list of names to search
    :return: int, the position of the name in the namelist (if found) else
             None
    """

    # clean name
    name = name.strip().upper()
    # strip namelist as char array
    namelist = np.char.array(namelist).strip().upper()
    # search for name in list
    if name in namelist:
        position = np.where(name == namelist)[0][0]
        # return position
        return position
    # if not found return None
    return None


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
    # get the original coords in SkyCoord
    ocoord = SkyCoord(ra, dec, unit='deg')
    # get gaia time and observation time
    gaia_time = Time('2015.5', format='decimalyear')
    obs_time = Time(mjd, format='mjd')
    # get entries as numpy arrays (with units)
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
    separation = coords1.separation(ocoord)
    # find the position of the minimum separated value
    position = np.argmin(separation.value)
    # return the position
    return position


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # load table
    table = Table.read(OBJECT_FILE, format='csv')

    if USE_DATABASE:

        from apero.core import constants
        from apero.core.core import drs_database

        # get object database
        ObjectDatabase = drs_database.ObjectDatabase

        params = constants.load('SPIROU')
        pconst = constants.pload('SPIROU')

        params.set('OBJ_LIST_CROSS_MATCH_RADIUS', value=CMRADIUS)
        params.set('OBJ_LIST_GAIA_PLX_LIM', value=0.5)
        params.set('OBJ_LIST_GOOGLE_SHEET_URL',
                   value='1jwlux8AJjBMMVrbg6LszJIpFJrk6alhbT5HA7BiAHD8')
        params.set('OBJ_LIST_GOOGLE_SHEET_WNUM', value=0)

        # load database
        objdbm = ObjectDatabase(params)
        objdbm.load_db()

        # for now delete database
        columns, ctypes = pconst.OBJECT_DB_COLUMNS()
        objdbm.database.delete_table('MAIN')
        objdbm.database.add_table('MAIN', columns, ctypes)
        objdbm.load_db()
        outdict = dict()
    else:
        objdbm = None
        outdict = dict()

    # loop around objects (later)
    for row in range(len(table)):
        print('='*50)
        print('Accessing row {0} / {1}'.format(row + 1, len(table)))
        print('='*50)
        # get properties
        gaia_id = None
        ra = table['RA_DEG'][row]
        dec = table['DEC_DEG'][row]
        objname = table['OBJECT'][row]
        # mjdmid
        exptime_days = (table['EXPTIME'][row] / 2.0) / (3600 * 24.0)
        mjdmid = table['MJDEND'][row] - exptime_days
        # set up an object instance for this target
        astro_obj = AstroObject(params, pconst, gaia_id, ra, dec, objdbm,
                                objname, 0.0, 0.0, 0.0, 0.0, None)
        # resolve target (from gaia id or ra/dec)
        astro_obj.resolve_target(mjdmid)

        # get simbad aliases for this object
        astro_obj.get_simbad_aliases()

        # write to database
        if USE_DATABASE:
            astro_obj.write_obj(objdbm)
        else:
            astro_obj.write_table(outdict)

        print(astro_obj)

    # need to write table
    if not USE_DATABASE:
        # make table
        outtable = Table()
        # add columns
        for key in outdict:
            outtable[key] = outdict[key]
        # sort by parallax
        with warnings.catch_warnings(record=True) as _:
            sortmask = np.argsort(outtable['PLX'])
        outtable = outtable[sortmask]
        # write to file
        outtable.write('test_obj_database.csv', format='csv', overwrite=True)


# =============================================================================
# End of code
# =============================================================================
