import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm
import etienne_tools as et
from astropy.table import Table
import os
from scipy.optimize import curve_fit
import warnings

def fit_gauss(x, y, p0):
    fit, pcov = curve_fit(gauss, x, y, p0=p0)
    return fit


def gauss(x, cen, ew, amp, zp):
    return np.exp(-0.5 * (x - cen) ** 2 / ew ** 2) * amp + zp

def remove_bgnd_e2ds_fp(e2ds):
    with warnings.catch_warnings(record=True) as _:
        for iord in range(e2ds.shape[0]):
            sp = np.array(e2ds[iord])
            low = np.zeros_like(sp)
            high = np.zeros_like(sp)

            for reg in range(e2ds.shape[1]):
                i1 = reg-50
                i2 = reg+50
                if i1 <0:
                    i1 = 0
                if i2>e2ds.shape[1]:
                    i2 = e2ds.shape[1]

                pp = np.nanpercentile(sp[i1:i2],[5,95])
                low[reg] = pp[0]
                high[reg] = pp[1]

            low = et.lowpassfilter(low, 50)
            high = et.lowpassfilter(high, 50)

            e2ds[iord] = (sp-low)/(high-low)

    return e2ds

"""
Things to do to use this code :

Create the "NIRPS_drift_folder" directory and change the variable below.

Put all your e2ds A and B files for a sequence in a given sub-folder (e.g., data_010203)
This will be your data set name afterward

1) We expect the e2ds files to be names "something_A.fits" "something_B.fits". If not,
   change the glob.glob line below so that files from each fiber are properly found by
   the code

1a) You may have to change the MJDATE keyword in the code depending on how dates are saved in the 
    e2ds headers.

2) A ref_table is created for each fiber, it will keep track of where FP peaks are located

3) Once we have the reference table, we take each e2ds and fit a gaussian to each line

4) Once we have all line positions, we dump them in a table for each e2ds file

5) The position tables are read and we get a difference to the median position for each line (~16000 lines)

6) We do some stats to find the mean motion of lines. The propagation of errors is also done here.

"""

# parameters to update
NIRPS_drift_folder = '/Volumes/courlan/nirps_drift/'
#datasets = ['data1','data2','data3','data4']
datasets = ['data1_geneva','data2_geneva','data3_geneva','data4_geneva']


for dataset in datasets:
    if not os.path.isdir(NIRPS_drift_folder + dataset):
        raise ValueError(
            et.color('\n\nHey, the folder {} redoes not exist!\nThis is bad bad bad !!!'.format(NIRPS_drift_folder + dataset),
                     'red'))

    for fiber in ['A', 'B']:
        # line to be changed is your e2ds files are not name 'something_A.fits' and 'something_B.fits'
        files = np.array(glob.glob('{0}{1}/*{2}.fits'.format(NIRPS_drift_folder, dataset, fiber)))

        if 'geneva' not in dataset:
            # creation of the file that gives the position of each FP peak to 1 pixel for the gaussian fit
            ref_table = 'tbl_FP_{0}_peaks.csv'.format(fiber)
        else:
            # creation of the file that gives the position of each FP peak to 1 pixel for the gaussian fit
            ref_table = 'tbl_FP_{0}_geneva_peaks.csv'.format(fiber)

        # only create if necessary
        if not os.path.isfile(ref_table):
            # median of all e2ds files for the fiber
            sp = fits.getdata(files[0])
            cube = np.zeros([sp.shape[0], sp.shape[1], len(files)])
            for i in tqdm(range(len(files))):
                cube[:, :, i] = fits.getdata(files[i])
                cube[:, :, i] /= np.nanpercentile(cube[:, :, i], 95)

            # median FP for that fiber
            ref = np.nanmedian(cube, axis=2)

            # find lines as domains between two consecutive minima
            all_x = []
            all_iord = []
            all_step = []
            for ord in tqdm(range(ref.shape[0]),leave = False):
                sp = ref[ord]
                g = np.where((np.gradient(sp)[0:-1] > 0) * (np.gradient(sp)[1:] < 0))[0]

                fit, _ = et.robust_polyfit(g, np.gradient(g), 3, 3)

                iord = np.zeros_like(g, dtype=int) + ord

                all_step = np.append(all_step, np.polyval(fit, g))
                all_x = np.append(all_x, g)
                all_iord = np.append(all_iord, iord)

            # write the reference table
            tbl = Table()
            tbl['step'] = all_step
            tbl['xpos'] = np.array(all_x, dtype=int)
            tbl['order'] = np.array(all_iord, dtype=int)
            tbl['xpos_measured'] = np.zeros_like(all_step, dtype=float)

            tbl.write(ref_table)

        # for each file from each fiber, measure the position of lines with a gaussian fit
        xpix = np.arange(4088)  # pixel position along order
        for i in tqdm(range(len(files)),leave = False):
            # read the table
            tbl_file = files[i].split('.fits')[0] + '.csv'
            tbl = Table.read(ref_table)
            # set position along order to NaN in case the fit fails
            tbl['xpos_measured'] = np.nan
            if not os.path.isfile(tbl_file):
                tbl = et.td_convert(Table.read(ref_table))

                sp = remove_bgnd_e2ds_fp(fits.getdata(files[i]))
                ord_prev = -1

                # loop through FP peaks and fit gaussians
                for j in tqdm(range(len(tbl['order'])),leave = False):
                    if tbl['order'][j] != ord_prev:
                        tmp = sp[int(tbl['order'][j])]
                        ord_prev = np.array(int(tbl['order'][j]))

                    # start and end of the FP peak domain
                    i1 = int(tbl['xpos'][j] - tbl['step'][j] / 2)
                    i2 = int(tbl['xpos'][j] + tbl['step'][j] / 2) + 1

                    if np.isfinite(np.sum(tmp[i1:i2])):
                        # fit but allow for errors
                        try:
                            p0 = (i1 + i2) / 2., 2, np.max(tmp[i1:i2]) - np.min(tmp[i1:i2]), 0
                            with warnings.catch_warnings(record=True) as _:
                                fit = fit_gauss(xpix[i1:i2], tmp[i1:i2], p0)
                            tbl['xpos_measured'][j] = fit[0]
                        except:
                            _ = 1

                # save table to a file
                tbl = et.td_convert(tbl)
                tbl.write(tbl_file)


    # create the compilation plot
    fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True,figsize = [16,8])
    for fiber in ['A', 'B']:
        # loop through fibers and define colors
        if fiber == 'A':
            fmt = 'g.'

            tbl_drifts = Table()
            tbl_drifts['files'] = files
            tbl_drifts['date'] = np.zeros_like(files)
            tbl_drifts['rv_A'] = np.zeros_like(files,dtype = float)
            tbl_drifts['rv_B'] = np.zeros_like(files,dtype = float)

            tbl_drifts['rv_left_A'] = np.zeros_like(files,dtype = float)
            tbl_drifts['rv_left_B'] = np.zeros_like(files,dtype = float)

            tbl_drifts['rv_right_A'] = np.zeros_like(files,dtype = float)
            tbl_drifts['rv_right_B'] = np.zeros_like(files,dtype = float)

            tbl_drifts['sig_A'] = np.zeros_like(files,dtype = float)
            tbl_drifts['sig_B'] = np.zeros_like(files,dtype = float)

            tbl_drifts['sig_right_A'] = np.zeros_like(files,dtype = float)
            tbl_drifts['sig_right_B'] = np.zeros_like(files,dtype = float)

            tbl_drifts['sig_left_A'] = np.zeros_like(files,dtype = float)
            tbl_drifts['sig_left_B'] = np.zeros_like(files,dtype = float)

            tbl_drifts['HIERARCH ESO PRO REC1 CAL9 NAME'] = np.zeros_like(files,dtype = '<U99')

            tbl_drifts['mjdate'] = np.zeros_like(files,dtype = float)

        if fiber == 'B':
            fmt = 'r.'



        # label for legend
        label = 'Fiber {}'.format(fiber)

        files = np.array(glob.glob('{0}/*{1}.fits'.format(dataset, fiber)))

        if 'geneva' not in dataset:
            # creation of the file that gives the position of each FP peak to 1 pixel for the gaussian fit
            ref_table = 'tbl_FP_{0}_peaks.csv'.format(fiber)
        else:
            # creation of the file that gives the position of each FP peak to 1 pixel for the gaussian fit
            ref_table = 'tbl_FP_{0}_geneva_peaks.csv'.format(fiber)


        # read the reference file
        tbl = Table.read(ref_table)
        # TODO: be removed, just to match the number of orders in Montreal and Geneva
        tbl = et.td_convert(tbl[tbl['order'] < 41])


        # cube that contains all lines for all files
        cube = np.zeros([len(tbl['order']), len(files)])
        mjdate = np.zeros_like(files, dtype=float)
        for i in tqdm(range(len(files)),leave = False):
            # get the MJDATE from header
            hdr = fits.getheader(files[i])

            if 'MJDATE' in hdr:
                # maybe change the MJDATE keyword if needed
                mjdate[i] = hdr['MJDATE']
            if 'MJD-OBS' in hdr:
                # maybe change the MJDATE keyword if needed
                mjdate[i] = hdr['MJD-OBS']

            tbl_drifts['date'][i] = hdr['DATE-OBS'].split('T')[0]

            if 'HIERARCH ESO PRO REC1 CAL9 NAME' in hdr:
                tbl_drifts['HIERARCH ESO PRO REC1 CAL9 NAME'][i] = hdr['HIERARCH ESO PRO REC1 CAL9 NAME']

            # read the table that contains line positions
            tbl_file = files[i].split('.fits')[0] + '.csv'



            tbl = Table.read(tbl_file)
            # TODO: be removed, just to match the number of orders in Montreal and Geneva
            tbl = tbl[tbl['order']<41]


            # add to the cube with all lines
            tmp = np.array(tbl['xpos_measured'])
            tmp[tmp == 0] = np.nan
            cube[:, i] = tmp

        # find the stddev of each line position through the entire sequence
        # we approximate the stddev as 1/2 the difference between the 84th and 16th percentile
        # of positions.
        v1, v2 = np.nanpercentile(cube, [16, 84], axis=1)

        sigma = (v2 - v1) / 2

        # median line position
        med = np.nanmedian(cube, axis=1)

        # place holder for the average difference to median position
        moy = np.zeros(len(files), dtype=float)
        moy_left = np.zeros(len(files), dtype=float)
        moy_right = np.zeros(len(files), dtype=float)

        # same for errors
        err = np.zeros(len(files), dtype=float)
        err_left = np.zeros(len(files), dtype=float)
        err_right = np.zeros(len(files), dtype=float)

        left = tbl['xpos']<2044
        right = tbl['xpos']>=2044

        # loop through files and get the mean line offset for each file
        for i in tqdm(range(len(files)),leave = False):
            cube[:, i] -= med

        #plt.close()
        #plt.imshow(cube, aspect='auto', vmin=-0.01, vmax=0.01)
        #plt.show()

        for i in tqdm(range(len(files)), leave=False):
            moy0, err0 = et.odd_ratio_mean(cube[:,i], sigma, 1e-4)
            moy[i], err[i] = moy0, err0

        cube_left = cube[left,:]
        med_left = np.nanmedian(cube_left, axis=1)
        sigma_left = sigma[left]
        # loop through files and get the mean line offset for each file
        for i in range(len(files)):
            moy0, err0 = et.odd_ratio_mean(cube_left[:, i] - med_left, sigma_left, 1e-4)
            moy_left[i], err_left[i] = moy0, err0


        cube_right = cube[right,:]
        med_right = np.nanmedian(cube_right, axis=1)
        sigma_right = sigma[right]
        # loop through files and get the mean line offset for each file
        for i in tqdm(range(len(files)),leave = False):
            moy0, err0 = et.odd_ratio_mean(cube_right[:, i] - med_right, sigma_right, 1e-4)
            moy_right[i], err_right[i] = moy0, err0

        if fiber == 'A':
            tbl_drifts['rv_A'] = moy
            tbl_drifts['sig_A'] = err

            tbl_drifts['rv_right_A'] = moy_right
            tbl_drifts['sig_right_A'] = err_right

            tbl_drifts['rv_left_A'] = moy_left
            tbl_drifts['sig_left_A'] = err_left

            tbl_drifts['mjdate'] = mjdate

        if fiber == 'B':
            tbl_drifts['rv_B'] = moy
            tbl_drifts['sig_B'] = err

            tbl_drifts['rv_right_B'] = moy_right
            tbl_drifts['sig_right_B'] = err_right

            tbl_drifts['rv_left_B'] = moy_left
            tbl_drifts['sig_left_B'] = err_left


        ax[2].errorbar(mjdate,(moy_right-moy_left)*1000, yerr = np.sqrt(err_right**2+err_left**2)*1000,
                       alpha=0.5,fmt=fmt, label='left-right, fib='+fiber)

        # plot the timeserie
        ax[0].errorbar(mjdate, moy*1000, yerr=err*1000, fmt=fmt, alpha=.5, label=label)

        # keep track of A and B fibers
        if fiber == 'A':
            rv_A = np.array(moy)
            err_A = np.array(err)

        if fiber == 'B':
            rv_B = np.array(moy)
            err_B = np.array(err)

    # plot difference between A and B
    fmt = 'k.'
    ax[1].errorbar(mjdate, (rv_A - rv_B)*1000, yerr=np.sqrt(err_A ** 2 + err_B ** 2)*1000, fmt=fmt, alpha=.5, label='A-B')

    ax[0].set(xlabel='MJD-OBS', ylabel='mean fp position [m/s]',title='Absolute drift')
    ax[1].set(xlabel='MJD-OBS', ylabel='mean fp position [m/s]',title='Science to calib drift')
    ax[2].set(xlabel='MJD-OBS', ylabel='mean fp position [m/s]',title='left-to-right drift')
    for i in range(3):
        ax[i].grid(linestyle='--', color='grey', alpha=0.3)
        ax[i].legend()

    plt.tight_layout()

    print('slope in m/s/day {}, fiber A'.format(np.polyfit(mjdate,rv_A,1)[0]))
    print('slope in m/s/day {}, fiber B'.format(np.polyfit(mjdate,rv_B,1)[0]))

    print('Fractional focal plane change : {}'.format((np.polyfit(mjdate,moy_right-moy_left,1)[0]/2044)))
    print('RMS in pix of left-to-right : {}'.format(np.nanstd(moy_right-moy_left)))

    plt.savefig('{}_fp_drift.pdf'.format(dataset))
    plt.close()
    #plt.show()

    tbl_drifts.write('{}_fp_drift.csv'.format(dataset),overwrite = True)

