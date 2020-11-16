import os
import requests
import numpy as np
from astropy.table import Table

URL_BASE = ('https://docs.google.com/spreadsheets/d/'
            '{}/gviz/tq?tqx=out:csv&sheet={}')
SHEET_ID = '1gvMp1nHmEcKCUpxsTxkx-5m115mLuQIGHhxJCyVoZCM'
WORKSHEET = 0


def odo_from_path(patharr):
    """
    Vectorized numpy string slicing in case lists get too long
    """
    sep = os.path.sep
    a = np.array(np.char.split(patharr, sep).tolist())
    a = a[:, -1]
    a = np.char.ljust(a, 7)
    return a.astype(int)


def check_blacklist(file_list, localfile=None, url=None, odo_kwd='ODOMETER'):
    """
    Pass a file list and remove files that have their odometer
    included in the bad odometer list (google sheet or local document).
    Args:
        file_list (array): List of odometers to check
        localfile (str): Local file to use as bad odometer list. Overrides
                         any URL when not None. Default: None
        url (str): URL to google sheet document with bad odometer list.
                   Overrides default SLS url.
        odo_kwd (str): Odometer keyword in bad odo list file.
    Returns:
        good_list (list): list of odometers that passed the check.
    """

    file_list = np.array(file_list)
    odo_list = odo_from_path(file_list)  # Keep only odometer

    # Load list
    if localfile is not None:
        bad_odo = Table.read(localfile)
    else:
        if url is None:
            url = URL_BASE.format(SHEET_ID, WORKSHEET)
        data = requests.get(url)
        bad_odo = Table.read(data.text, format='ascii')

    bad_odo = np.array(bad_odo[odo_kwd].tolist())

    # find which odometers are in bad list
    # we then use this mask to filter file names
    bad_mask = np.isin(odo_list, bad_odo)

    # remove bad files from list
    if np.sum(bad_mask):
        print('We remove the following files : ')
        bad_files = file_list[bad_mask]
        for bad_file in bad_files:
            print('\t{0}'.format(bad_file))
        print()
    else:
        print('No bad odometers listed.')

    return file_list[~bad_mask]
