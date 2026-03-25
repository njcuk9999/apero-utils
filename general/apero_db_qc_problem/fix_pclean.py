#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Code to add QCC_PASS = 1 to all pclean files

Created on 2026-03-25 14:10:00

@author: cook
"""
import multiprocessing as mp
import os

from astropy.io import fits
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================
FILE_PATHS = ['/cosmos99/spirou/apero-data/spirou_offline/red/',
              '/cosmos99/spirou/apero-data/spirou_offline/tellu/']

FILE_SUFFIX = 'pclean_AB.fits'

KW_QCC_PASS = 'QCC_ALL'

PARAM_TABLE_QCC_KEY = 'PASSED_ALL_QC'

PRECLEAN_QCC_KEY_PREFIX = 'TQCCP'

N_WORKERS = int(os.getenv('PCLEAN_NWORKERS', '20'))
CHUNKSIZE = int(os.getenv('PCLEAN_CHUNKSIZE', '25'))


# =============================================================================
# Define functions
# =============================================================================
def find_files(file_path):

    matched_files = []
    progress = tqdm(os.walk(file_path), desc='Scanning', unit='dir')
    for root, dirs, dir_files in progress:
        for file in dir_files:
            if file.endswith(FILE_SUFFIX):
                matched_files.append(os.path.join(root, file))
        progress.set_postfix(found=len(matched_files))

    print('Found {0} total files'.format(len(matched_files)))

    return matched_files


def get_qcc_pass_from_header(header) -> int:
    """Read all preclean QCC PASS keys from a FITS header."""
    pass_values = []
    for key in header:
        if key.startswith(PRECLEAN_QCC_KEY_PREFIX):
            pass_values.append(int(header[key]))

    if len(pass_values) == 0:
        return 0

    return min(pass_values)


def process_file(filename: str):
    """Process one FITS file and return a status tuple."""
    try:
        with fits.open(filename, mode='update') as hdul:
            header = hdul[0].header

            qcc_value = get_qcc_pass_from_header(header)
            header[KW_QCC_PASS] = qcc_value

            table = hdul['PARAM_TABLE'].data
            for row in range(len(table)):
                # this fixes a mistake I introduced
                if str(table['NAME'][row]) == '0':
                    table['NAME'][row] = PARAM_TABLE_QCC_KEY
                    table['VALUE'][row] = qcc_value == 1
                # this is the original fix for those files I didn't get to yet
                elif PARAM_TABLE_QCC_KEY in table['NAME'][row]:
                    table['VALUE'][row] = qcc_value == 1

            hdul.flush()

        return 'updated', filename, ''
    except Exception as exc:
        return 'error', filename, str(exc)


def process_files_parallel(files):
    """Process files in parallel and report aggregate status counts."""
    stats = {'updated': 0, 'skipped': 0, 'error': 0}
    errors = []

    if len(files) == 0:
        return stats, errors

    with mp.Pool(processes=N_WORKERS) as pool:
        iterator = pool.imap_unordered(process_file, files, chunksize=CHUNKSIZE)
        with tqdm(total=len(files), desc='Applying fix', unit='file') as pbar:
            for status, filename, message in iterator:
                stats[status] += 1
                if status == 'error':
                    errors.append((filename, message))

                pbar.update(1)
                pbar.set_postfix(updated=stats['updated'],
                                 skipped=stats['skipped'],
                                 error=stats['error'])

    return stats, errors


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # loop around filepaths
    for filepath in FILE_PATHS:
        print('=' * 70)
        print(f'Analysing {filepath}')
        print('=' * 70)

        files = find_files(filepath)
        stats, errors = process_files_parallel(files)

        print('Done for {0}'.format(filepath))
        print('  updated = {0}'.format(stats['updated']))
        print('  skipped = {0}'.format(stats['skipped']))
        print('  error   = {0}'.format(stats['error']))

        if len(errors) > 0:
            print('First 10 errors:')
            for filename, message in errors[:10]:
                print(f'  {filename}: {message}')

# =============================================================================
# End of code
# =============================================================================

