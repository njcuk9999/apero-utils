"""
Utility functions for the getrv.py script
"""
import os
import re

import numpy as np
from astropy.io import fits


def process_list_clarg(argument):
    """
    Process CLI argument that has a list of values which should be used as
    set of possible strings
    Args:
        argument (list): argument outputed by the parser.
    Returns
        newarg (str):
    """
    if argument == ['all']:
        newarg = '.*'
    elif len(argument) == 1:
        newarg = argument[0]
    else:
        newarg = '(?:'+'|'.join(argument)+')'

    return newarg


def get_ccf_files(path, obj, mask, sanit):
    """
    Get list of CCF file names matching all possible combinations of obj, mask,
    and sanit that correspond to existing files.
    Args:
        path (str): Path to all CCF subdirs
        obj (list of str): List of object names, can be ['all']
        mask (list of str): List of mask names, can be ['all']
        sanit (list of str): List of sanitization methods, can be ['all']
    Returns:
        fnames (list of str): List of file paths matching the input
            information.
    """
    # Pattern: OBJECT/*_SANIT_*_MASK_AB.fits
    # Bash: {obj1,obj2,...}/*_{sanit1,sanit2}_*_{mask1,mask2,...}_AB.fits
    # Bash: could have * instead of lists
    # RegEx: '(obj1|obj2|obj3)/.*_(sanit1|sanit2)_*_(mask1|mask2)_AB.fits'
    # RegEx: could have .* instead of lists, maybe [^/]*

    # Make sure we have proper list objects
    if obj == ['all']:
        obj = os.listdir(path)

    # Convert sanit and mask to regex strings
    sanit = process_list_clarg(sanit)
    mask = process_list_clarg(mask)

    # Get regex matchgin filenames
    regex = re.compile('.*_{}_.*_{}_AB.fits'.format(sanit, mask))

    # Only loop objects we are interested in
    # This will speed things up compared to os.walk and including obj in regex
    ccf_files = []
    for o in obj:

        # Get all files for one object
        flist = os.listdir(os.path.join(path, o))

        # regex will search single string, so put one line per file
        match_list = regex.findall('\n'.join(flist))

        # Add new files
        for match in match_list:
            ccf_files.append(os.path.join(path, o, match))

    return sorted(ccf_files)


def get_ccf_dict(ccf_files):
    """
    Get dictionnary mapping timeseries ID to list of files.
    Args:
        ccf_files (list of str): list of all ccf files.
    Return:
        ccf_dict (dict): dictionnary with IDs as keys and lists as elements.
    """
    file_ids = []
    file_paths = []
    for i, ccf_file in enumerate(ccf_files):
        hdr = fits.getheader(ccf_file, 1)

        # Get info from header
        obj = hdr['OBJECT']
        mask = hdr['CCFMASK'].split('.')[0]
        drs_version = hdr['VERSION']

        sani_sub = 'SANI'  # substring in some keys only if sanitized
        is_sanit = np.any([sani_sub in k for k in hdr.keys()])
        if is_sanit:
            sanit = 'sani'
        else:
            sanit = 'tcorr'
        reduction_id = '{0}__{1}__{2}__{3}'.format(obj,
                                                   mask,
                                                   sanit,
                                                   drs_version)
        file_ids.append(reduction_id)
        file_paths.append(ccf_file)

    file_ids = np.array(file_ids)
    file_paths = np.array(file_paths)
    uids = np.unique(file_ids)

    ccf_dict = {}
    for uid in uids:
        id_mask = file_ids == uid
        ccf_dict[uid] = file_paths[id_mask]

    return ccf_dict
