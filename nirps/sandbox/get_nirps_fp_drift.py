import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm
import etienne_tools as et
from astropy.table import Table
import os
from scipy.optimize import curve_fit


def fit_gauss(x, y, p0):
    fit, pcov = curve_fit(gauss, x, y, p0=p0)
    return fit


def gauss(x, cen, ew, amp, zp):
    return np.exp(-0.5 * (x - cen) ** 2 / ew ** 2) * amp + zp


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
dataset = 'datax'

if not os.path.isdir(NIRPS_drift_folder + dataset):
    raise ValueError(
        et.color('\n\nHey, the folder {} does not exist!\nThis is bad bad bad !!!'.format(NIRPS_drift_folder + dataset),
                 'red'))

for fiber in ['A', 'B']:
    # line to be changed is your e2ds files are not name 'something_A.fits' and 'something_B.fits'
    files = np.array(glob.glob('{0}{1}/*{2}.fits'.format(NIRPS_drift_folder, dataset, fiber)))

    # creation of the file that gives the position of each FP peak to 1 pixel for the gaussian fit
    ref_table = 'tbl_FP_{0}_peaks.csv'.format(fiber)

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
        for ord in tqdm(range(ref.shape[0])):
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
    for i in tqdm(range(len(files))):
        # read the table
        tbl_file = files[i].split('.')[0] + '.csv'
        tbl = Table.read(ref_table)
        # set position along order to NaN in case the fit fails
        tbl['xpos_measured'] = np.nan
        if not os.path.isfile(tbl_file):
            tbl = et.td_convert(Table.read(ref_table))

            sp = fits.getdata(files[i])
            ord_prev = -1

            # loop through FP peaks and fit gaussians
            for j in tqdm(range(len(tbl['order']))):
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
                        fit = fit_gauss(xpix[i1:i2], tmp[i1:i2], p0)
                        tbl['xpos_measured'][j] = fit[0]
                    except:
                        _ = 1

            # save table to a file
            tbl = et.td_convert(tbl)
            tbl.write(tbl_file)

# create the compilation plot
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
for fiber in ['A', 'B']:
    # loop through fibers and define colors
    if fiber == 'A':
        fmt = 'g.'
    if fiber == 'B':
        fmt = 'r.'

    # label for legend
    label = 'Fiber {}'.format(fiber)

    files = np.array(glob.glob('{0}/*{1}.fits'.format(dataset, fiber)))
    ref_table = 'tbl_FP_{0}_peaks.csv'.format(fiber)

    # read the reference file
    tbl = et.td_convert(Table.read(ref_table))

    # cube that contains all lines for all files
    cube = np.zeros([len(tbl['order']), len(files)])
    mjdate = np.zeros_like(files, dtype=float)
    for i in tqdm(range(len(files))):
        # get the MJDATE from header
        hdr = fits.getheader(files[i])

        # maybe change the MJDATE keyword if needed
        mjdate[i] = hdr['MJDATE']

        # read the table that contains line positions
        tbl_file = files[i].split('.')[0] + '.csv'
        tbl = Table.read(tbl_file)
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
    # same for errors
    err = np.zeros(len(files), dtype=float)

    # loop through files and get the mean line offset for each file
    for i in range(len(files)):
        moy0, err0 = et.odd_ratio_mean(cube[:, i] - med, sigma, 1e-4)
        moy[i], err[i] = moy0, err0

    # plot the timeserie
    ax[0].errorbar(mjdate, moy, yerr=err, fmt=fmt, alpha=.5, label=label)

    # keep track of A and B fibers
    if fiber == 'A':
        rv_A = np.array(moy)
        err_A = np.array(err)

    if fiber == 'B':
        rv_B = np.array(moy)
        err_B = np.array(err)

# plot difference between A and B
fmt = 'k.'
ax[1].errorbar(mjdate, rv_A - rv_B, yerr=np.sqrt(err_A ** 2 + err_B ** 2), fmt=fmt, alpha=.5, label='A-B')

ax[0].set(xlabel='MJD-OBS', ylabel='mean fp position [pix]')
ax[1].set(xlabel='MJD-OBS', ylabel='mean fp position [pix]')
ax[0].grid(linestyle='--', color='grey', alpha=0.3)
ax[1].grid(linestyle='--', color='grey', alpha=0.3)
plt.tight_layout()
ax[0].legend()
ax[1].legend()
plt.savefig('{}_fp_drift.pdf'.format(dataset))
plt.show()


