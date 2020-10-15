"""
Created on Oct. 14 2020

Author: Thomas Vandal <thomas.vandal@umontreal.ca>

The spirou_ccf_analysis script by Eder Martioli was used as a starting point
for this script (https://github.com/edermartioli/spirou-ccf-analysis).

The tools used under the hood to calculate RV timeseries from CCFs were written
by Étienne Artigau and edited by Eder Martioli and Thomas Vandal.
"""
from argparse import ArgumentParser


# CLI arguments
parser = ArgumentParser(
        description='Analyze series of SPIRou CCF data products to obtain '
                    'radial velocity'
        )

# Arguments used to find files
parser.add_argument(
        '-p',
        '--path',
        default='.',
        help='Path to directory containing data file. This is CWD by default.'
        )
parser.add_argument(
        '-o',
        '--object',
        nargs='+',
        default='all',
        help='Name of objects to analyze. If all, will consider all '
             'directories in path and use their names as object names.'
        )
parser.add_argument(
        '-m',
        '--mask',
        nargs='+',
        default='all',
        help='Masks used to calculate CCF. If all, all available masks for '
             'each object is used. If one of the masks is not available for a '
             'given object, it is skipped.'
        )
parser.add_argument(
        '-s',
        '--sanitize',
        choices=['sanit', 'tcorr', 'all'],
        default='all',
        help='Use sanitized (sanit) or non-sanitized (tcorr) data, '
             'or all available options (all).'
        )
parser.add_argument(
        '-i',
        '--pattern',
        help='Input CCF data pattern relative to path. This overrides object, '
             'mask and sanitize arguments.'
        )
parser.add_argument(
        '-d',
        '--drs-version',
        dest='drs_version',
        default='all',
        help='DRS version (X.X.XX) used to reduced data. If all, any version '
        'is accepted. Otherwise only files reduced with specified '
        'version are processed.'
        )

# Arguments passed to get_object_rv
parser.add_argument(
        '-r',
        '--rv-method',
        dest='rv_method',
        default='all',
        help='Method used to calculate RVs. Must be either "all", "template" '
             '"gaussian" or "bisector_M_N" where M and N are percentiles. '
        )  # choices are checked in get_object_rv
parser.add_argument(
        '-e',
        '--exclude-orders',
        dest='exclude_orders',
        nargs='+',
        type=int,
        default=[-1],
        help='Orders between 0 and 49 to exclude from the analysis. If no '
             'orders are supplied, none are excluded.'
        )
parser.add_argument(
        '-b',
        '--bandpass',
        nargs='+',
        default=['Y', 'J', 'H', 'K'],
        choices=['Y', 'J', 'H', 'K'],
        help='Bandpass to include in analyis. Any combination of YJHK works.'
        )
parser.add_argument(
        '-n',
        '--min-snr',
        dest='min_snr',
        type=float,
        default=0.0,
        help='Minimum SNR',
        )
parser.add_argument(
        '-v',
        '--verbose',
        help='Print debug info.'
        )

clargs = parser.parse_args()

# Get list of files
