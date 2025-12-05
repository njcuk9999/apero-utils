from astropy.io import fits
from astroquery import simbad
from astropy.table import Table
import numpy as np
from astropy.coordinates import SkyCoord
import os


# Retrieve Simbad names and header names for a list of FITS files
def retrieve_simbad_name(files):
    simbad_names = []  # List to store Simbad names
    hdr_names = []     # List to store header OBJECT names

    for file in files:
        h = fits.getheader(file)  # Read FITS header
        ra = h['RA']              # Get RA from header
        dec = h['DEC']            # Get DEC from header
        hdr_name = h['OBJECT']    # Get OBJECT name from header

        c = SkyCoord(ra, dec, unit="deg")  # Create SkyCoord object

        # Query Simbad for H and V magnitudes in a small region around the coordinates
        simbad.Simbad.add_votable_fields('flux(H)', 'flux(V)')
        result = simbad.Simbad.query_region(c, radius='0d1m0s')

        # Rename columns for easier access
        if 'FLUX_H' in result.colnames:
            result.rename_column('FLUX_H', 'H')
        if 'FLUX_V' in result.colnames:
            result.rename_column('FLUX_V', 'V')
        # Capitalize all column names
        for col in result.colnames:
            if col != col.upper():
                result.rename_column(col, col.upper())

        # Set masked values to 99 for H and V bands
        if 'H' in result.colnames:
            result['H'].fill_value = 99
        if 'V' in result.colnames:
            result['V'].fill_value = 99

        imin = np.argmin(result['H'])  # Find index of minimum H mag
        name = result['MAIN_ID'][imin] # Get Simbad name for that entry

        simbad_names.append(name)
        hdr_names.append(hdr_name)

    return simbad_names, hdr_names



# Parse file names from user input, keeping only existing files
def parse_file_names():
    entries = []
    entry = ' '
    print("Enter a file name (or press enter to finish): ")
    while entry != '':
        entry = input()
        if entry != '':
            entries.append(entry)

    entries = np.array(entries)

    # Split entries by spaces to handle pasted lists
    # the input is in a random format. We just join everything with spaces and split
    # at spaces. We assume that filenames will be surrounded by ' ' and the rest is rubbish
    v = ' '.join(entries)
    v = v.split(' ')
    v = np.array(v)

    # Keep only files that exist
    keep = np.zeros_like(v, dtype=bool)
    for i in range(len(v)):
        if os.path.exists(v[i]):
            keep[i] = True
    v = v[keep]

    return v

files = parse_file_names()

simbad_names, hdr_names = retrieve_simbad_name(files)

# Print unique Simbad names and their corresponding header names
unames = np.unique(simbad_names)
for name in unames:
    print(name,'\t', hdr_names[simbad_names.index(name)])