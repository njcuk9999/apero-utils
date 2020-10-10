import os
from shutil import copy
import sys
import warnings
import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from astropy.table import Table
from astropy.time import Time
from astropy.io import fits
from tqdm import tqdm

# Find which methods are from these
from bisector import bisector
from fits2wave import fits2wave
from check_blacklist import check_blacklist


def dispatch_object(
        obj,
        in_dir='all_ccfs',
        out_parent=None,
        verbose=False
        ):
    """ Dispatch CCFs per object
    Find all CCF Files and copy those matching 'obj' into a folder with the
    object's name.
    Args:
        obj (str): String with object name. If this is 'all', all objects in
                   input dir will be dispatched.
        in_dir (str): Path to directory where CCF files are stored.
        out_parent (str): Path to parent of output directories. A subdirectory
                          for each object will be created here. If none,
                          parent directory of in_dir is used.
        verbose (bool): print messages telling what the script does
    """
    # Get list of CCF files
    all_ccf_files = glob.glob(os.path.join(in_dir, '*ccf*AB.fits'))

    # Handle output dir
    if out_parent is None:
        os.path.abspath(os.path.join(in_dir, os.pardir))

    if verbose:
        print(
            'We are dispatching CCF files into the proper directory for '
            'object: {0}'.format(obj)
            )
    ngood_AB = 0
    ngood_C = 0
    for ccf_file in tqdm(all_ccf_files):
        # read ccf file
        im, h = fits.getdata(ccf_file, header=True)

        # Output directory file
        outdir = os.path.join(out_parent, h['DRSOBJN'])
        outname = os.path.join(outdir, os.path.basename(ccf_file))

        if (h['DRSOBJN'] == obj) or (obj == 'all'):
            # If the directory with the obj name does not exist, we create it
            if not os.path.isdir(outdir):
                os.makedirs(outdir)  # Create output dir and parents
            if not os.path.isfile(outname):
                copy(ccf_file, outname)
            if '_AB.fits' in ccf_file:
                ngood_AB += 1
            if '_C.fits' in ccf_file:
                ngood_C += 1

    print('We found {0} files for that obj. '
          'It includes {1} AB files and {2} C drift files'.format(
                ngood_AB + ngood_C, ngood_AB, ngood_C
                )
          )


def get_object_rv(obj,
                  ccf_parent=None,
                  mask='sept18_andres_trans50',
                  method='template',
                  exclude_orders=[-1],
                  weight_table='',
                  force=True,
                  snr_min=0.0,
                  weight_type='',
                  sanitize=False,
                  bandpass='YJHK',
                  velocity_window=10,
                  dvmax_per_order=1.0,
                  doplot=True,
                  do_blacklist=False,
                  detailed_output=False
                  ):
    """ Get RV Timeseries for a given object

    Args:
        obj (str): name of object to be analyzed. Should be the exact same name
            as the directory where files for this object are contained
        ccf_parent (str): Parent directory with per-object CCF directories. If
            None, this is taken to be the current working directory.
            Default: None
        mask (str): Name of the mask used for the CCF.
            Default: sept18_andres_trans50
        method (str): Method used to measure velocity.
            Default: template
            Supported methods:
                - template: default method, TODO
                - bisector_{N}_{M}: returns velocity between N and Mth
                  percentiles of line depth.
                - gaussfit: fits a Gaussian to the mean CCF
                - all: use all methods and store RV for each of them
        exclude_orders (list): list of orders to exclude systematically. When
            set to [-1], no orders are removed
            Default: [-1]
        weight_table (str): path to file with weight_table for orders. If no
            table is provided. The table must have 49 rows, one column with
            name "WEIGHT", it must be formatted as a proper astropy table.
            Default: '' (Empty string)
        force (bool): If set to False and CCF table exists, it is simply read
            and returned.
            Default: False
        snr_min (float): Set a threshold below which CCFs are rejected.
            The extracted SNR for order 35 is used as a threshold.
            Default: 0.0
        weight_types: TODO (Still experimental)
        sanitize (bool): If true, look for CCFs of santiized file with _sani_
            in the name. Indicates that sanitizing ocde as been run on the
            data. Otherwise, file have _tcorr_ in place of _sani_.
            Default: False
        bandpass (str): Which photometric bandpass to look at (any of 'YJHK').
            Default: 'YJHK'
        velocity_window (float): Window in km/s used around the CCF minimum
            used to measure velocity. We go from -window to +window. When
            adjusting median CCF to CCF of one observation, (in template
            method), this is the width over which the dot product between the
            derivative and the residuals are computed. When determining the Q
            value of the CCF per band, this is also the width of the integral.
            Default: 10.0
        dvmax_per_order (float): Reject any order which CCF minimum beyond this
            dvmax (in km/s). For very low SNR CCF in some orders, the median
            CCF is not a (near) Gaussian line and the minimum can be way off.
            This is used to reject those orders. User gets a warning when this
            happens.
            Default: 1.0
        doplot (True): Whether the code should produced plots.
            Default: True
        do_blacklist (True): Check if the input files are blacklisted.
            Default: False (because this adds some overheads)
    """

    if method == 'all':
        method = 'gaussian_template_bisector_20_80'

    if sanitize:
        sp_type = 'sani'
    else:
        sp_type = 'tcorr'

    dict_ccf = dict()

    if ccf_parent is None:
        ccf_parent = os.getcwd()

    fpattern = '*{0}*{1}_AB.fits'.format(sp_type, mask)
    ccf_files = np.array(
            glob.glob(os.path.join(ccf_parent, obj, fpattern))
            )

    if do_blacklist:
        ccf_files = check_blacklist(ccf_files)

    # get typical wavelength solution from first file and central wavelength
    # grid per order
    wave_middle = np.nanmean(
            fits2wave(fits.getheader(ccf_files[0], ext=1)), axis=1
            )

    keep_orders = np.zeros(49)
    if 'Y' in bandpass:
        keep_orders[(wave_middle > 938)*(wave_middle < 1113)] = True
    if 'J' in bandpass:
        keep_orders[(wave_middle > 1153)*(wave_middle < 1354)] = True
    if 'H' in bandpass:
        keep_orders[(wave_middle > 1462)*(wave_middle < 1808)] = True
    if 'K' in bandpass:
        keep_orders[(wave_middle > 1957)*(wave_middle < 2400)] = True

    for i in range(49):
        if i in exclude_orders:
            keep_orders[i] = False

    # update orders rejected because of domain as well as input requirements.
    exclude_orders = np.where(np.logical_not(keep_orders))[0]

    for i in range(49):
        if i in exclude_orders:
            keep_orders[i] = False

    # form a unique batch name with mask, obj and method
    batch_name = '{0}_mask_{1}_{2}'.format(obj, mask, method)
    batch_name = os.path.join(ccf_parent, batch_name)

    if not force:
        if os.path.isfile('{0}.csv'.format(batch_name)):
            return Table.read('{0}.csv'.format(batch_name))

    tbl = Table()  # output table to be saved as CSV file with RV measurements
    tbl['FILES'] = ccf_files

    # to keep track of the unique designation of the file
    tbl['ODOMETER'] = np.zeros_like(tbl, dtype='U7')
    for i in range(len(tbl)):
        tbl['ODOMETER'][i] = tbl['FILES'][i].split('/')[-1].split('o')[0]

    tbl['RV'] = np.zeros_like(ccf_files, dtype=float)  # measured RV
    tbl['ERROR_RV'] = np.zeros_like(
            ccf_files, dtype=float  # measured RV error
            )
    tbl['CCF_RESIDUAL_RMS'] = np.zeros_like(
            ccf_files, dtype=float  # RMS of CCF - median(CCF)
            )

    # keywords from file headers to be added to the CSV table.
    keywords = [
            'MJDATE',
            'BERV',
            'RV_DRIFT',
            'EXTSN035',
            'AIRMASS',
            'TLPEH2O',
            'TLPEOTR',
            'RV_WAVFP',
            'RV_SIMFP',
            'DATE',
            'MJDMID',
            'DATE-OBS',
            'EXPTIME'
            ]

    # add method-specific keywords
    if 'bisector' in method:
        # Vt and Vb from Perryman
        tbl['RV_BIS'] = np.zeros_like(ccf_files, dtype=float)  # mid point
        tbl['BIS_SLOPE'] = np.zeros_like(ccf_files, dtype=float)  # slope
        tbl['BIS_WIDTH'] = np.zeros_like(ccf_files, dtype=float)  # width
        tbl['Vt'] = np.zeros_like(ccf_files, dtype=float)  # velocity 'top'
        tbl['Vb'] = np.zeros_like(ccf_files, dtype=float)  # velocity 'bottom'
        tbl['BIS'] = np.zeros_like(ccf_files, dtype=float)  # velocity width

    # add method-specific keywords
    if 'gaussian' in method:
        tbl['RV_GAUSS'] = np.zeros_like(ccf_files, dtype=float)  # mean velo.
        tbl['GAUSS_WIDTH'] = np.zeros_like(ccf_files, dtype=float)  # width
        tbl['GAUSS_AMP'] = np.zeros_like(ccf_files, dtype=float)  # depth
        tbl['GAUSS_ZP'] = np.zeros_like(ccf_files, dtype=float)  # zp

    # Feel free to add columns here

    npy_file = '{0}_{1}_{2}_ccf_cube.npy'.format(obj, sp_type, mask)
    npy_file = os.path.join(obj, npy_file)  # Add dir to path

    if not os.path.isfile(npy_file):
        print('we load all CCFs into one big cube')
        ccf_RV_previous = None
        for i in (range(len(ccf_files))):
            # we loop through all files
            ccf_tbl = fits.getdata(ccf_files[i])

            # ccf velocity offset, not to be confused with measured RV
            ccf_RV = ccf_tbl['RV']

            # We must absolutely have always the same RV grid for the CCF.
            # We must check that consecutive ccf_RV are identical
            if i != 0:
                if np.sum(ccf_RV != ccf_RV_previous):
                    print('We have a big problem!'
                          'The RV vector of CCF files are not all the same')
                    print(
                        'Files {0} and {1} have different RV vectors.'.format(
                            ccf_files[i-1], ccf_files[i]
                            )
                        )
                    sys.exit()
            ccf_RV_previous = np.array(ccf_RV)

            print('V[min/max] {0:.1f} / {1:.1f} km/s, file {2}'.format(
                np.min(ccf_RV), np.max(ccf_RV), ccf_files[i]
                ))
            # if this is the first file, we create a cube that contains all
            # CCFs for all orders for all files
            if i == 0:
                ccf_cube = np.zeros([49, len(ccf_tbl), len(ccf_files)])+np.nan

            # we input the CCFs in the CCF cube
            for j in range(49):
                tmp = ccf_tbl['ORDER'+str(j).zfill(2)]

                if False not in np.isfinite(tmp):
                    # we normalize to a continuum of 1
                    tmp /= np.polyval(np.polyfit(ccf_RV, tmp, 1), ccf_RV)
                    ccf_cube[j, :, i] = tmp

        print('We save {0}, this will speed things up '
              'next time you run this code'.format(npy_file))
        np.save(npy_file, ccf_cube)

    else:
        print('We load {0}, this is speedier'.format(npy_file))
        ccf_cube = np.load(npy_file)

        # we need to load the first file just to get the velocity grid
        ccf_tbl = fits.getdata(ccf_files[0])
        ccf_RV = ccf_tbl['RV']

    for j in range(49):
        # if we need to exlude orders, we do it here.
        if j in exclude_orders:
            ccf_cube[j, :, :] = np.nan

    print('We load values from headers, slower on first run, faster later')
    for i in (range(len(ccf_files))):
        hdr = fits.getheader(ccf_files[i], ext=1)
        if i == 0:
            # Now that we have a first header,
            # we add the relevant columns to the CSV table
            for key in keywords:
                if key in hdr:

                    key_type = type(hdr[key])
                    # if we have a string, we set the table to accept long
                    # values (up to 99 characters)
                    if key_type == str:
                        key_type = '<U99'
                else:
                    # keyword not in header, we need to assume something.
                    # The safest is string
                    key_type = str

                # add the column to the CSV file
                tbl[key] = np.zeros_like(ccf_files, dtype=key_type)

        for key in keywords:
            if key in hdr:
                tbl[key][i] = hdr[key]

    # we apply the SNR threshold
    keep = tbl['EXTSN035'] > snr_min
    tbl = tbl[keep]
    ccf_cube = ccf_cube[:, :, keep]
    ccf_files = ccf_files[keep]

    with warnings.catch_warnings(record=True) as _:
        # some slices in the sum are NaNs, that's OK
        # this is the median CCF  all epochs for the 49 orders
        med_ccf = np.nanmedian(ccf_cube, axis=2)

    for iord in range(49):
        if iord not in exclude_orders:
            if not np.isfinite(np.mean(med_ccf[iord, :])):
                print('Order {0} has a CCF full of NaNs. '
                      'Added to the rejected orders '.format(iord))
                exclude_orders = np.append(exclude_orders, iord)

    # Find minimum for CCF. This is used to fit a gaussian to each order
    # and force velocity to zero
    id_min = np.nanargmin(np.nanmedian(med_ccf, axis=0))

    # find if a given CCF is off by more than the pre-defined threshold
    dv_CCF_min = (ccf_RV[np.argmin(med_ccf, axis=1)] - ccf_RV[id_min])
    bad_orders = dvmax_per_order < np.abs(dv_CCF_min)

    for iord in range(49):
        if iord not in exclude_orders:
            if bad_orders[iord]:
                print('The CCF of order {0} has its minima {1:.2f} km/s '
                      'from median CCF, above threshold of +-{2:.2f} km/s'
                      .format(iord, dv_CCF_min[iord], dvmax_per_order))
                exclude_orders = np.append(exclude_orders, iord)

    ccf_cube[exclude_orders, :, :] = np.nan

    # Find valid pixels to measure CCF properties
    g = np.abs(ccf_RV - ccf_RV[id_min]) < velocity_window

    with warnings.catch_warnings(record=True) as _:
        # Some slices in the sum are NaNs, that's OK
        ccf_Q = np.nansum(np.gradient(med_ccf[:, g], axis=1)**2, axis=1)

    ccf_Q[ccf_Q == 0] = np.nan
    ccf_depth = 1-med_ccf[:, id_min]

    if weight_type == 'DVRMS_CC':
        raise ValueError('DVRMS_CC not available yet')
        # Commented out. We can implement this on separate branch
        # weights = 1/np.nanmedian(DVRMS_CC, axis=1)**2
        # weights[np.isfinite(weights) == False] = 0

        # for iord in range(49):
        #     if iord in exclude_orders:
        #         weights[iord] = 0

        # weights = weights/np.sum(weights)
    else:
        if weight_table == "":
            # now we find the RMS of the Nth spectrum relative to the median
            rms = np.zeros([len(ccf_files), 49])
            for i in range(len(ccf_files)):
                with warnings.catch_warnings(record=True) as _:
                    # Some slices in the median are NaNs, that's OK
                    rms[i, :] = np.nanmedian(
                            np.abs(ccf_cube[:, :, i] - med_ccf),
                            axis=1
                            )
                rms[i, :] /= np.nanmedian(rms[i, :])

            rms[:, exclude_orders] = np.nan

            if doplot:
                vmin = np.nanpercentile(rms, 3)
                vmax = np.nanpercentile(rms, 97)
                plt.imshow(rms, aspect='auto', vmin=vmin, vmax=vmax)
                plt.xlabel('Nth order')
                plt.ylabel('Nth frame')
                plt.title('RMS of CCF relative to median')
                plt.show()

            with warnings.catch_warnings(record=True) as _:
                # some slices in the sum are NaNs, that's OK
                # this is the typical noise from the ccf dispersion
                ccf_rms = np.nanmedian(rms, axis=0)

            # set to NaN values that are invalid
            ccf_rms[ccf_rms == 0] = np.nan

            # Assuming that the CCF has the same depth everywhere,
            # this is the correct weighting of orders
            weights = ccf_Q/ccf_rms**2
            weights[weights == 0] = np.nan
            weights[exclude_orders] = np.nan
            # we normalize the sum of the weights to one
            weights /= np.nansum(weights)

            if doplot:
                fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True)
                ax[0].plot(weights, 'go')
                ax[0].set(title='{0}, mask {1}'.format(obj, mask),
                          xlabel='Nth order', ylabel='Relative order weight')

                ax[1].plot(ccf_Q, 'go')
                ax[1].set(xlabel='Nth order', ylabel='ccf Q')

                ax[2].plot(1/ccf_rms**2, 'go')
                ax[2].set(xlabel='Nth order', ylabel='1/$\\sigma_{CCF}^2$')
                plt.tight_layout()
                plt.savefig('{0}_weights.pdf'.format(batch_name))
                plt.show()

            tbl_weights = Table()
            tbl_weights['order'] = np.arange(49)
            tbl_weights['weights'] = weights
            tbl_weights['ccf_depth'] = ccf_depth
            tbl_weights['ccf_Q'] = ccf_Q
            tbl_weights.write('{0}_weights.csv'.format(batch_name),
                              overwrite=True)

        else:
            print('You provided a weight file, we load it and'
                  'apply weights accordingly')
            tbl_weights = Table.read(weight_table)
            weights = np.array(tbl_weights['weights'], dtype=float)
            weights /= np.nansum(weights)

    if doplot:
        plt.imshow(med_ccf, aspect='auto', vmin=0.8, vmax=1.05,
                   extent=[np.min(ccf_RV), np.max(ccf_RV), 49, 0])
        plt.xlabel('Velocity bin [km/s] ')
        plt.ylabel('Nth order')
        plt.title('Median CCF')
        plt.savefig('{0}_medianccf.pdf'.format(batch_name))
        plt.show()

        plt.imshow(ccf_cube[:, :, 0]-med_ccf, aspect='auto', vmin=-0.1,
                   vmax=0.1, extent=[np.min(ccf_RV), np.max(ccf_RV), 49, 0])
        plt.xlabel('Velocity bin [km/s]')
        plt.ylabel('Nth order')
        plt.title('Sample residual CCF map')
        plt.savefig('{0}_residualccf.pdf'.format(batch_name))
        plt.show()

        plt.plot(tbl['MJDATE'], tbl['EXTSN035'], 'g.')
        plt.xlabel('MJDATE')
        plt.ylabel('SNR for order 35\n(around 1.6 $\\mu$m)')
        plt.title('Signal-to-noise ratio')
        plt.savefig('{0}_snr.pdf'.format(batch_name))
        plt.show()

    ccf_cube_norm = np.zeros_like(ccf_cube)
    for i in range(49):
        if np.isfinite(weights[i]):
            ccf_cube_norm[i, :, :] = (ccf_cube[i, :, :] * weights[i])

    # Get a per-file weighted mean CCF
    mean_ccf = np.nansum(ccf_cube_norm, axis=0)

    if doplot:
        fig, ax = plt.subplots(nrows=1, ncols=1)
        for i in range(len(ccf_files)):
            color = [i/len(ccf_files), 1-i/len(ccf_files), 1-i/len(ccf_files)]
            ax.plot(ccf_RV, mean_ccf[:, i], color=color, alpha=0.2)

        ax.set(xlabel='Velocity [km/s]', ylabel='CCF depth', title='Mean CCFs')
        plt.tight_layout()
        plt.savefig('{0}_CCFs.pdf'.format(batch_name))
        plt.show()

    g = np.abs(ccf_RV - ccf_RV[id_min]) < velocity_window

    rv_prev = np.array(tbl['RV'])

    # We use the bisector method, the name should be something like
    # method = 'bisector_30_70' to get the mean bisector between the
    # 30th and 70th percentile
    if 'bisector' in method:
        # we find the min of CCF and will only compute bisector of +-50 km/s
        # to avoid problems at ccf edges
        imin = np.argmin(np.nanmedian(mean_ccf, axis=1))

        # just get the parameters after bisector
        params_bis = method.split('bisector_')[1]
        bis_min, bis_max = np.array(
                params_bis.split('_')[0:2], dtype=float
                ) / 100.

        for i in range(len(ccf_files)):

            # Try:
            depth, bis, width = bisector(ccf_RV, mean_ccf[:, i],
                                         low_high_cut=0.2)
            fit = np.polyfit(
                    depth[(depth > bis_min) & (depth < bis_max)]
                    - (bis_min + bis_max) / 2,
                    bis[(depth > bis_min) & (depth < bis_max)],
                    1)

            tbl['RV'][i] = fit[1]
            # just in case you want to have both bisector and
            # template, we keep a RV that is specific to this method
            tbl['RV_BIS'][i] = fit[1]
            tbl['BIS_SLOPE'][i] = fit[0]
            tbl['BIS_WIDTH'][i] = np.mean(
                    width[(depth > bis_min) & (depth < bis_max)]
                    )

            # mean 'top' CCF between 55 and 80% of depth
            tbl['Vt'][i] = np.mean(bis[(depth > 0.55)*(depth < 0.80)])

            # mean 'bottom' CCF between 20-40%
            tbl['Vb'][i] = np.mean(bis[(depth > 0.20)*(depth < 0.40)])
            tbl['BIS'][i] = tbl['Vt'][i] - tbl['Vb'][i]

            if False:
                print('We had an error with file {0} computing the bisector'
                      .format(ccf_files[i]))
                print('Values will be reported as NaN')
                tbl['RV'][i] = np.nan
                tbl['BIS_SLOPE'][i] = np.nan
                tbl['BIS_WIDTH'][i] = np.nan

        fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
        ax[0].plot(tbl['MJDATE'], tbl['RV'], 'g.')
        ax[0].set(title='Velocity', xlabel='MJDATE', ylabel='RV [km/s]')
        ax[1].plot(tbl['MJDATE'], tbl['BIS_SLOPE'], 'g.')
        ax[1].set(title='Bisector slope', xlabel='MJDATE',
                  ylabel='slope [km/s/fract. depth]')
        plt.tight_layout()
        plt.savefig('{0}_RV.pdf'.format(batch_name))
        plt.show()

    if 'gaussian' in method:
        imin = np.argmin(np.nanmedian(mean_ccf, axis=1))

        for i in range(len(ccf_files)):
            p0 = [ccf_RV[imin], 1, 1, -0.1]
            fit, pcov = curve_fit(gauss, ccf_RV, mean_ccf[:, i], p0=p0)

            tbl['RV'][i] = fit[0]

            # just in case you want to have gauss/bisector and
            # template, we keep a RV that is specific to this method
            tbl['RV_GAUSS'][i] = fit[0]

            tbl['GAUSS_WIDTH'][i] = fit[1]
            tbl['GAUSS_AMP'][i] = fit[3]
            tbl['GAUSS_ZP'][i] = fit[2]

        fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
        ax[0].plot(tbl['MJDATE'], tbl['RV'], 'g.')
        ax[0].set(title='Velocity', xlabel='MJDATE', ylabel='RV [km/s]')
        ax[1].plot(tbl['MJDATE'], tbl['GAUSS_WIDTH']*2.354, 'g.')
        ax[1].set(title='Gaussian width', xlabel='MJDATE',
                  ylabel='Gaussian FWHM [km/s]')
        plt.tight_layout()
        plt.savefig('{0}_RV.pdf'.format(batch_name))
        plt.show()

    # fitting a 'Template' ... this is always done.
    nite_max = 20
    ite = 0
    rms_rv_ite = np.inf
    # we iterate until we have an rms from iteration to iteration of <10 cm/s
    # or we reached a max of 20 iterations
    print('\n')

    corr_ccf = np.array(mean_ccf)

    fig, ax = plt.subplots(nrows=1, ncols=2)

    # funky scaling of imshow
    vmin = np.nanpercentile(corr_ccf, 3)
    vmax = np.nanpercentile(corr_ccf, 97)
    ax[0].imshow(corr_ccf,  aspect='auto', vmin=vmin, vmax=vmax,
                 extent=[0, len(ccf_files), np.min(ccf_RV), np.max(ccf_RV)])
    ax[0].set(xlabel='Nth observation', ylabel='Velocity [km/s]',
              title='Before CCF register')

    per_ccf_rms = np.ones(len(ccf_files))
    while (rms_rv_ite > 1e-4) and (ite < nite_max):
        if ite == 0:
            tbl['RV'] = 0

        w = 1/per_ccf_rms**2
        w /= np.sum(w)
        med_corr_ccf = np.zeros(len(ccf_RV))
        for i in range(len(w)):
            med_corr_ccf += (corr_ccf[:, i] * w[i])

        # normalize continuum to 1
        continuum = np.abs(ccf_RV-ccf_RV[id_min]) > velocity_window
        med_corr_ccf /= np.nanmedian(med_corr_ccf[continuum])

        fit = np.polyfit(ccf_RV[continuum], med_corr_ccf[continuum], 2)
        corr = np.polyval(fit, ccf_RV)
        corr -= np.mean(corr)
        med_corr_ccf -= corr

        for i in range(len(ccf_files)):
            spline = ius(ccf_RV, mean_ccf[:, i], ext=3, k=5)
            corr_ccf[:, i] = spline(ccf_RV+tbl['RV'][i])

            # correcting median of CCF
            med = np.nanmedian(corr_ccf[:, i] - med_corr_ccf)
            mean_ccf[:, i] -= med

            # correcting depth of CCF
            amp = (np.nansum((corr_ccf[:, i] - np.mean(corr_ccf[:, i]))
                             * (med_corr_ccf - np.mean(med_corr_ccf)))
                   / np.nansum((med_corr_ccf - np.mean(med_corr_ccf))**2)
                   )
            mean_ccf[:, i] = (
                    (mean_ccf[:, i] - np.mean(mean_ccf[:, i]))/np.sqrt(amp)
                    + np.mean(mean_ccf[:, i])
                    )

            # correcting 2rd order polynomial structures in continuum
            fit = np.polyfit(ccf_RV, med_corr_ccf-corr_ccf[:, i], 2)

            corr = np.polyval(fit, ccf_RV)
            mean_ccf[:, i] += corr/2

        deriv = np.gradient(med_corr_ccf) / np.gradient(ccf_RV)
        deriv = deriv[g]
        deriv = deriv / np.nansum(deriv ** 2)

        for i in range(len(ccf_files)):
            residu = corr_ccf[:, i] - med_corr_ccf
            per_ccf_rms[i] = np.nanstd(residu)
            tbl['RV'][i] -= np.nansum(residu[g]*deriv)

        tbl['RV'] -= np.nanmean(tbl['RV'])
        # plt.plot( tbl['RV'],'.')
        rms_rv_ite = np.nanstd(rv_prev - tbl['RV'])
        print('Template CCF iteration number {0:3}, '
              'rms RV change {1:3.4f} km/s for this step'
              .format(ite+1, rms_rv_ite))
        rv_prev = np.array(tbl['RV'])
        ite += 1

    tbl['RV_TEMPLATE'] = np.array(tbl['RV'])

    vmin = np.nanpercentile(corr_ccf, 3)
    vmax = np.nanpercentile(corr_ccf, 97)
    ax[1].imshow(corr_ccf, aspect='auto', vmin=vmin, vmax=vmax,
                 extent=[0, len(ccf_files), np.min(ccf_RV), np.max(ccf_RV)])
    ax[1].set(xlabel='Nth observation',
              ylabel='Velocity [km/s]',
              title='After CCF register')
    plt.show()

    # we get the systemic velocity from the BISECTOR between 0.3 and 0.7 depth
    depth, bis, width = bisector(ccf_RV, med_corr_ccf, low_high_cut=0.3,
                                 doplot=True,
                                 figure_title='mean CCF\ndebug plot',
                                 ccf_plot_file='ccf_{0}.pdf'.format(obj))
    tbl['RV'] += np.nanmedian(bis)

    if doplot:
        fig, ax = plt.subplots(nrows=2, ncols=1)
        for i in range(len(ccf_files)):
            color = [i/len(ccf_files), 1-i/len(ccf_files), 1-i/len(ccf_files)]
            ax[0].plot(ccf_RV, corr_ccf[:, i], color=color, alpha=0.2)
            ax[1].plot(ccf_RV, corr_ccf[:, i], color=color, alpha=0.2)

        ax[0].set(xlabel='Velocity [km/s]', ylabel='CCF depth',
                  title='Mean CCFs')
        ax[1].set(xlabel='Velocity [km/s]',
                  ylabel='CCF depth',
                  title='Mean CCFs',
                  xlim=[ccf_RV[id_min]-10, ccf_RV[id_min]+10])
        plt.tight_layout()
        plt.savefig('{0}_template.pdf'.format(batch_name))
        plt.show()

    # we add a measurement of the STDDEV of each mean CCF relative to the
    # median CCF after correcting for the measured velocity.
    # If you are going to add 'methods', add them before this line
    corr_ccf = np.array(mean_ccf)
    for i in range(len(ccf_files)):
        spline = ius(ccf_RV, mean_ccf[:, i], ext=3)
        corr_ccf[:, i] = spline(ccf_RV+tbl['RV'][i]-np.mean(tbl['RV']))

    med_corr_ccf = np.nanmedian(corr_ccf, axis=1)
    g = np.abs(ccf_RV - ccf_RV[id_min]) < 10

    if doplot:
        plt.plot(ccf_RV, med_corr_ccf, color='black', alpha=0.4,
                 label='median CCF', linewidth=2)

    # We compute the projection of the CCF residuals onto the second
    # and third derivatives of the CCF
    d2 = np.gradient(np.gradient(med_corr_ccf) / np.gradient(ccf_RV))
    d3 = np.gradient(np.gradient(
        np.gradient(med_corr_ccf) / np.gradient(ccf_RV)
        ))
    # second derivatives
    tbl['D2_RESIDUAL_CCF'] = np.zeros_like(tbl, dtype=float)
    # third derivatives
    tbl['D3_RESIDUAL_CCF'] = np.zeros_like(tbl, dtype=float)

    # pix scale expressed in CCF pixels
    # SPIRou pixels are about 2.3 km/s
    pix_scale = 2.3/np.nanmedian(np.gradient(ccf_RV))
    for i in range(len(ccf_files)):
        residual = corr_ccf[:, i] - med_corr_ccf

        tbl['D2_RESIDUAL_CCF'][i] = np.nansum(residual*d2)/np.nansum(d2)
        tbl['D3_RESIDUAL_CCF'][i] = np.nansum(residual*d3)/np.nansum(d3)

        tbl['CCF_RESIDUAL_RMS'][i] = np.std(residual[g])

        # 1/dvrms -avoids division by zero
        inv_dvrms = ((np.gradient(med_corr_ccf) / np.gradient(ccf_RV))
                     / ((np.nanstd(residual) * np.sqrt(pix_scale)))
                     )
        tbl['ERROR_RV'][i] = 1 / np.sqrt(np.nansum(inv_dvrms ** 2))

        color = [i/len(ccf_files), 1-i/len(ccf_files), 1-i/len(ccf_files)]
        if doplot:
            plt.plot(ccf_RV, residual+1, color=color, alpha=0.2)

    if doplot:
        plt.title('Residual CCFs')
        plt.xlabel('velocity [km/s]')
        plt.ylabel('CCF depth')
        plt.legend()
        plt.savefig('{0}_residual_CCF.pdf'.format(batch_name))

        plt.show()

    if doplot:
        t3 = Time(tbl['MJDATE'], format='mjd')
        plt.plot_date(t3.plot_date, tbl['D2_RESIDUAL_CCF'], 'go')
        plt.title('Second derivative \n activity indicator')
        plt.xlabel('Date')
        plt.ylabel('CCF residual projection on\nCCF 2nd derivative')
        plt.savefig('{0}_d2_activity.pdf'.format(batch_name))

        plt.show()

    # output to csv file
    tbl.write('{0}.csv'.format(batch_name), overwrite=True)

    if not detailed_output:
        return tbl
    else:

        dict_ccf['TABLE_CCF'] = tbl
        dict_ccf['MEAN_CCF'] = mean_ccf

        return dict_ccf


def gauss(v, v0, ew, zp, amp):
    # gaussian with a constant offset.
    # As we know that the ccfs are negative structures, amp will be negative
    return zp+amp*np.exp(-0.5*(v-v0)**2/ew**2)
