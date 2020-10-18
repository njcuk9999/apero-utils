"""
Created on Oct. 14 2020

Author: Thomas Vandal <thomas.vandal@umontreal.ca>

The spirou_ccf_analysis script by Eder Martioli was used as a starting point
for this script (https://github.com/edermartioli/spirou-ccf-analysis).

The tools used under the hood to calculate RV timeseries from CCFs were written
by Étienne Artigau and edited by Eder Martioli and Thomas Vandal.
"""
import os
import glob
from argparse import ArgumentParser

import numpy as np

from ccf2rv import get_object_rv
import getrv_utils as gu


# CLI arguments
parser = ArgumentParser(
        description='Analyze series of SPIRou CCF data products to obtain '
                    'radial velocity'
        )

# Arguments used to find files
parser.add_argument(
        '-i',
        '--pattern',
        help='Input CCF data pattern relative to path. This overrides object, '
             'mask and sanitize arguments.'
        )
parser.add_argument(
        '-t',
        '--save-to',
        dest='outdir',
        default=None,
        help='Directory the data is saved to. The default is None, in which '
             'case the data is saved to the parent directory of CCF files'
        )
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
        default=['all'],
        help='Name of objects to analyze. If all, will consider all '
             'directories in path and use their names as object names.'
        )
parser.add_argument(
        '-m',
        '--mask',
        nargs='+',
        default=['all'],
        help='Masks used to calculate CCF. If all, all available masks for '
             'each object is used. If one of the masks is not available for a '
             'given object, it is skipped.'
        )
parser.add_argument(
        '-s',
        '--sanitize',
        nargs='+',
        choices=['sani', 'tcorr', 'all'],
        default=['all'],
        help='Use sanitized (sani) or non-sanitized (tcorr) data, '
             'or all available options (all).'
        )
# parser.add_argument(
#       '-d',
#       '--drs-version',
#       dest='drs_version',
#       default='all',
#       help='DRS version (X.X.XX) used to reduced data. If all, any version '
#       'is accepted. Otherwise only files reduced with specified '
#       'version are processed.'
#       )

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
        '--snr-min',
        dest='snr_min',
        type=float,
        default=0.0,
        help='Minimum SNR',
        )
parser.add_argument(
        '-w',
        '--velocity-window',
        dest='velocity_window',
        type=float,
        default=10.0,
        help='Velocity window to scan',
        )
parser.add_argument(
        '-x',
        '--dvmax-per-order',
        dest='dvmax_per_order',
        type=float,
        default=1.0,
        help='Max dv per order',
        )
parser.add_argument(
        '-f',
        '--no-force',
        dest='force',
        action='store_false',
        default=True,
        help='Do not force new calculation.'
        )
parser.add_argument(
        '--bin',
        action='store_true',
        default=False,
        help='Bin RV output'
        )
parser.add_argument(
        '--show-plots',
        dest='show_plots',
        action='store_true',
        default=False,
        help='Show plots'
        )
parser.add_argument(
        '--save-plots',
        dest='save_plots',
        action='store_true',
        default=False,
        help='Save plots'
        )
parser.add_argument(
        '--save-weights',
        dest='save_weights',
        action='store_true',
        default=False,
        help='Save weight table'
        )
parser.add_argument(
        '--save-table',
        dest='save_table',
        action='store_true',
        default=False,
        help='Save full result table'
        )
parser.add_argument(
        '--save-cube',
        dest='save_cube',
        action='store_true',
        default=False,
        help='Save CCF cube.'
        )
parser.add_argument(
        '--no-rv',
        dest='save_rv',
        action='store_false',
        default=True,
        help='Save RV timeseries.'
        )
parser.add_argument(
        '--no-plots',
        dest='do_plots',
        action='store_false',
        default=True,
        help='Completely skip plotting routines'
        )
parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=False,
        help='Print debug info.'
        )

clargs = parser.parse_args()

# If no pattern, build it from other information
# In both cases, we load a list of CCF files
if clargs.pattern is None:
    ccf_files = gu.get_ccf_files(
            clargs.path,
            clargs.object,
            clargs.mask,
            clargs.sanitize
            )

else:
    if clargs.verbose:
        print('WARNING: pattern overrides mask, object, and sanitize')
    pattern = os.path.join(clargs.path, clargs.pattern)
    ccf_files = sorted(glob.glob(pattern))

print(ccf_files)
exit()
# Get dictionnary where each key is the "ID" of a timeseries with given
# object, mask, sanitize, DRS version. The keys have format
# object__ mask__sani__version.
ccf_dict = gu.get_ccf_dict(ccf_files)

bandpass = ''.join(np.unique(clargs.bandpass))

for run_id in ccf_dict.keys():

    list_of_files = ccf_dict[run_id]

    if clargs.verbose:
        print('Processing batch {0} containing {1} files'
              .format(run_id, len(list_of_files))
              )

    tbl = get_object_rv(
            ccf_files=list_of_files,
            outdir=clargs.outdir,
            run_id=run_id,
            method=clargs.rv_method,
            exclude_orders=clargs.exclude_orders,
            force=clargs.force,
            snr_min=clargs.snr_min,
            bandpass=bandpass,
            velocity_window=clargs.velocity_window,
            dvmax_per_order=clargs.dvmax_per_order,
            save_ccf_cube=clargs.save_cube,
            save_result_table=clargs.save_table,
            save_rv_timeseries=clargs.save_rv,
            save_weight_table=clargs.save_weights,
            bin_rv_timeseries=clargs.bin,
            doplot=clargs.do_plots,
            showplots=clargs.show_plots,
            saveplots=clargs.save_plots,
            )
