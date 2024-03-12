# Pierrot, 2024-03-12

from typing import Any, Dict, List, Tuple
import numpy as np
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, Distance
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
import astropy.units as uu
import glob
import requests
import pandas as pd
from io import StringIO
import warnings

def ra_dec_to_alt_az(ra, dec, observer_latitude, observer_longitude, observation_time):
    # Create a SkyCoord object with RA and DEC
    coordinates = SkyCoord(ra=ra*uu.deg, dec=dec*uu.deg, frame='icrs')

    # Set observer's location
    observer_location = EarthLocation(observer_longitude*uu.deg, observer_latitude*uu.deg)

    # Set observation time
    observation_time = Time(observation_time)

    # Transform coordinates to Altitude and Azimuth
    altaz_coordinates = coordinates.transform_to(AltAz(obstime=observation_time, location=observer_location))

    # Return Altitude and Azimuth
    return altaz_coordinates.alt.deg*uu.deg, altaz_coordinates.az.deg*uu.deg

def coords_check(fits_file, RA, DEC, tol=1/60, verbose=False): 
    # Read in the fits file
    obj_data = fits.open(fits_file, hdu=1)
    header = obj_data[0].header

    # Fetch observation infos
    ALT, AZ = header['HIERARCH ESO TEL ALT']*uu.deg, header['HIERARCH ESO TEL AZ']*uu.deg
    date_obs = header['DATE-OBS']
    observer_latitude, observer_longitude = header['HIERARCH ESO TEL GEOLAT'], header['HIERARCH ESO TEL GEOLON']
    
    # Predict Altitude and Azimuth
    pred_ALT, pred_AZ = ra_dec_to_alt_az(RA, DEC, observer_latitude, observer_longitude, date_obs)
    pred_AZ = pred_AZ - 180*uu.deg # Southern hemisphere
    
    # Compute the error in position
    sin_error_ALT = np.abs(np.sin(ALT) - np.sin(pred_ALT)) # sin(x) = x for small x
    sin_error_AZ = np.abs(np.sin(AZ) - np.sin(pred_AZ))
    
    if verbose: print('Position error: delta ALT = {} deg, delta AZ = {} deg'.format(sin_error_ALT, sin_error_AZ))
    
    if sin_error_ALT > tol or sin_error_AZ > tol:
        if verbose: print('Position error is too large!')
        return False
    else: 
        if verbose: print('Position error is within tolerance.')
        return True
    

# ## Fetching the target's info 

obj_name = 'LHS1140'

MAIN_URL = ('https://docs.google.com/spreadsheets/d/'
            '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
            'export?format=csv&gid=0')
PENDING_URL = ('https://docs.google.com/spreadsheets/d/'
               '1dOogfEwC7wAagjVFdouB1Y1JdF9Eva4uDW6CTZ8x2FM/'
               'export?format=csv&gid=623506317')

# get main table
main_request = requests.get(MAIN_URL)
# open main table
main_dataframe = pd.read_csv(StringIO(main_request.text))
# get pending table
pending_request = requests.get(PENDING_URL)
# open pending table
pending_dataframe = pd.read_csv(StringIO(pending_request.text))
# merge these keeping all rows from main table and adding non-repeating
#  rows from pending table
astrom_dataframe = pd.concat([main_dataframe,
                                pending_dataframe]).drop_duplicates()

# if you prefer an astropy table
astrom_table = Table.from_pandas(astrom_dataframe)

itarget = np.where(astrom_table['OBJNAME'] == obj_name)[0][0]
target_table = astrom_table[itarget]


# Propogate the star w.r.t its proper motion

def propagate_coords(objdata: Dict[str, Any], obs_time: Time) -> Tuple[SkyCoord, Time]:
    # deal with distance
    if objdata['PLX'] <= 0:
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
    # get the coordinates
    with warnings.catch_warnings(record=True) as _:
        # center the image on the current coordinates
        curr_coords = coords.apply_space_motion(dt=delta_time)
        curr_time = obs_time
    # return some stuff
    return curr_coords, curr_time

# ## Checking the files of your target
# Fetch .fits files (TODO: Enter your file path here)
fits_files = glob.glob('../run_lbl/data_dir_HARPS/science/LHS1140/*.fits')

tolerance = 1/60 # 1 arcmin

success_array = []
for file in fits_files:
    time = Time(fits.getheader(file)['DATE-OBS'])
    coords, time = propagate_coords(target_table, time)
    success_array.append(coords_check(file, coords.ra.deg, coords.dec.deg, tol = np.arcsin(tolerance), verbose=False))
    
print('Success rate: {:.2f}%'.format(np.sum(success_array)/len(success_array)*100))
print('Problematic files:')
for i, file in enumerate(fits_files):
    if not success_array[i]: print(file)




