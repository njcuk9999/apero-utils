"""
Utility functions for the getrv.py script
"""
import os
import re
import yaml
import argparse

import numpy as np
from astropy.io import fits


DEFAULT_CFG = {
        'pattern': None,
        'outdir': None,
        'path': '.',
        'object': ['all'],
        'mask': ['all'],
        'sanitize': ['all'],
        'method': 'all',
        'exclude_orders': [-1],
        'bandpass': ['Y', 'J', 'H', 'K'],
        'snr_min': 0.0,
        'velocity_window': 10.0,
        'dvmax_per_order': 1.0,
        'force': True,
        'bin': False,
        'median': False,
        'do_plots': True,
        'show_plots': False,
        'verbose': False,
        'save_plots': False,
        'save_weights': False,
        'save_table': False,
        'save_cube': False,
        'save_rv': True,
        'do_bad_odo': False,
        }


def parse_clargs():
    parser = argparse.ArgumentParser(
            description='Analyze series of SPIRou CCF data products to obtain '
                        'radial velocity'
            )

    # Arguments used to find files
    parser.add_argument(
            '-c',
            '--config',
            help='Path to YAML configuration file describing the run.'
            )
    parser.add_argument(
            '-i',
            '--pattern',
            default=argparse.SUPPRESS,
            help='Input CCF data pattern relative to path. '
                 'This overrides object, mask and sanitize arguments. '
                 'Default is None'
            )
    parser.add_argument(
            '-t',
            '--save-to',
            dest='outdir',
            default=argparse.SUPPRESS,
            help='Directory the data is saved to. By default, '
                 'the data is saved to the parent directory of CCF files'
            )
    parser.add_argument(
            '-p',
            '--path',
            default=argparse.SUPPRESS,
            help='Path to directory containing data file. '
                 'This is CWD by default.'
            )
    parser.add_argument(
            '-o',
            '--object',
            nargs='+',
            default=argparse.SUPPRESS,
            help='Name of objects to analyze. If all, will consider all '
                 'directories in path and use their names as object names.'
                 'By default, all objects are used.'
            )
    parser.add_argument(
            '-m',
            '--mask',
            nargs='+',
            default=argparse.SUPPRESS,
            help='Masks used to calculate CCF. If all, all available masks '
                 'for each object are used. If one of the masks is not '
                 'available for a given object, it is ignored.'
                 'By default, all masks are used.'
            )
    parser.add_argument(
            '-s',
            '--sanitize',
            nargs='+',
            choices=['sani', 'tcorr', 'all'],
            default=argparse.SUPPRESS,
            help='Use sanitized (sani) or non-sanitized (tcorr) data, '
                 'or all available options (all, this is the default).'
            )
    # parser.add_argument(
    #       '-d',
    #       '--drs-version',
    #       dest='drs_version',
    #       default='all',
    #       help='DRS version (X.X.XX) used to reduced data. If all,
    #            'any version '
    #            'is accepted. Otherwise only files reduced with specified '
    #            'version are processed.'
    #       )

    # Arguments passed to get_object_rv
    parser.add_argument(
            '-r',
            '--rv-method',
            dest='method',
            default=argparse.SUPPRESS,
            help='Method used to calculate RVs. Must be either '
                 '"all" (default), "template", "gaussian" or "bisector_M_N" '
                 'where M and N are percentiles.'
            )  # choices are checked in get_object_rv
    parser.add_argument(
            '-e',
            '--exclude-orders',
            dest='exclude_orders',
            nargs='+',
            default=argparse.SUPPRESS,
            type=int,
            help='Orders between 0 and 49 to exclude from the analysis. If no '
                 'orders are supplied, none are excluded (default).'
            )
    parser.add_argument(
            '-b',
            '--bandpass',
            nargs='+',
            choices=['Y', 'J', 'H', 'K'],
            default=argparse.SUPPRESS,
            help='Bandpass to include in analyis. '
                 'Any combination of Y, J, H, K works.'
                 'All bands are used by default.'
            )
    parser.add_argument(
            '-n',
            '--snr-min',
            dest='snr_min',
            type=float,
            default=argparse.SUPPRESS,
            help='Minimum SNR (Default: 0.0).',
            )
    parser.add_argument(
            '-w',
            '--velocity-window',
            dest='velocity_window',
            type=float,
            default=argparse.SUPPRESS,
            help='Velocity window to scan (Default: 10.0)',
            )
    parser.add_argument(
            '-x',
            '--dvmax-per-order',
            dest='dvmax_per_order',
            type=float,
            default=argparse.SUPPRESS,
            help='Max dv per order (Default: 1.0)',
            )
    parser.add_argument(
            '-f',
            '--no-force',
            dest='force',
            action='store_false',
            default=argparse.SUPPRESS,
            help='Do not force new calculation.'
            )
    parser.add_argument(
            '--bin',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Bin RV output'
            )
    parser.add_argument(
            '--median',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Use median and error from MAD for binned RV timeseries. '
                 'By default, the weighted averge is used.'
            )
    parser.add_argument(
            '--show-plots',
            dest='show_plots',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Show plots'
            )
    parser.add_argument(
            '--save-plots',
            dest='save_plots',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Save plots'
            )
    parser.add_argument(
            '--save-weights',
            dest='save_weights',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Save weight table'
            )
    parser.add_argument(
            '--save-table',
            dest='save_table',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Save full result table'
            )
    parser.add_argument(
            '--save-cube',
            dest='save_cube',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Save CCF cube.'
            )
    parser.add_argument(
            '--no-rv',
            dest='save_rv',
            action='store_false',
            default=argparse.SUPPRESS,
            help='Save RV timeseries.'
            )
    parser.add_argument(
            '--no-plots',
            dest='do_plots',
            action='store_false',
            default=argparse.SUPPRESS,
            help='Completely skip plotting routines.'
            )
    parser.add_argument(
            '--bad-odo',
            dest='do_bad_odo',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Check bad odometer list.'
            )
    parser.add_argument(
            '-v',
            '--verbose',
            action='store_true',
            default=argparse.SUPPRESS,
            help='Print debug info.'
            )

    return parser.parse_args()


def get_config(clargs):
    """
    Generate config directory from CLI arguments.
    Args:
        clargs (argparse.Namespace): Namespace with CLI args.
    Returns:
        cfg (dict): Dictionnary with run information.
    """
    # Load config file
    if clargs.config is not None:
        with open(clargs.config, 'r') as ymlfile:
            cfg = yaml.load(ymlfile, Loader=yaml.FullLoader)

        # Make sure we have all keys
        for k in DEFAULT_CFG.keys():
            if k not in cfg:
                cfg[k] = DEFAULT_CFG[k]

    else:
        cfg = DEFAULT_CFG.copy()

    # Overwrite config file if key has as CLI arg
    clargs_dict = vars(clargs)
    for k in cfg.keys():
        if k in clargs:
            cfg[k] = clargs_dict[k]

    return cfg


def process_list_clarg(argument, to_re=True):
    """
    Process CLI argument that has a list of values
    Args:
        argument (list): argument outputed by the parser.
        to_re (bool): Convert list to regex matching elements
    Returns
        newarg (str):
    """
    # Make sure list
    if not isinstance(argument, list):
        argument = [argument]

    if to_re:
        if argument == ['all']:
            newarg = '.*'
        elif len(argument) == 1:
            newarg = argument[0]
        else:
            newarg = '(?:'+'|'.join(argument)+')'
    else:
        newarg = argument

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
    obj = process_list_clarg(obj, to_re=False)
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
