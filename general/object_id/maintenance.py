"""
Functions to facilitate maintenance of apero sheet.
Primarly designed for automated updates with update_

@author: vandalt
"""
import os
import glob

import tqdm
import pandas as pd
from astropy.io import fits

import newutils as ut

RAW_PATTERN = '/spirou/cfht_nights/common/raw/20*/*o.fits'


def get_object_info(fpattern, verbose=False):
    """
    Fetch object information in headers of raw data from CFHT
    """
    keys = [
            'OBJECT',
            'OBJRV',
            'OBJTEMP',
            'MJDEND',
            'RA_DEG',
            'DEC_DEG',
            ]

    # Initialize dict
    outdict = dict()
    for k in keys:
        outdict[k] = []

    # Handle files separately
    fnames = []

    # Loop through files
    # Slower than dfits | fitsort but avoids system calls
    if verbose:
        print("Fetching data from headers. This might take a few minutes.")
    for filepath in tqdm.tqdm(glob.glob(fpattern)):
        hdr = fits.getheader(filepath)
        for k in keys:
            ut.hdr_to_dict(outdict, k, hdr)
        fnames.append(os.path.basename(filepath))

    # Create dataframe from dict
    df = pd.DataFrame([])
    df['FILE'] = fnames
    for k in keys:
        df[k] = outdict[k]

    return df
