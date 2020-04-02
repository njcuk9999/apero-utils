import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.signal import medfilt

# Needed on Etienne's laptop to be consistent with his matplotlib version.
# Will not hurt on other machines with other versions. This is the
# default.
plt.ioff()


# Function to get a measurement of the peak of the blaze in NIRPS images.
# We need a FLAT_FLAT and a DARK_DARK. The DARK_DARK is used to mask
# pixels that are too hot (>10 sigma), and we use the FLAT_FLAT to determine
# a Nth (default N=85..97%) percentile profile in the direction of cross-dispersion.
# These profiles give the overall shape of the blaze. A bisector of  the blaze
# is then derived and we return the median mid-point of the blaze cutouts
# as well as the std of the values. This is used to guide the adjustment of the
# grating of NIRPS.
#
# Do *not* forget that the mid-point of the blaze has to be centered on the mid-point
# of the order, *not* the center of the array. The two would be identical if we had a
# constant dispersion (in velocity) along the orders, but as there is a ~10% variation
# in the dispersion along orders, the peak of the blaze should be slightly on the bluer
# side of orders. The optimal position (from optical design) is at X ~ 1925.
#
# INPUTS --
#
# find_peak_blaze(file_flat, file_dark, percentiles = [85,97])
#
# file_flat --> file with a FLAT_FLAT image
# file_dark --> file with a DARK_DARK image to be used for masking
# percentiles -> cutouts percentiles. Leave 85-97 unless you have a
#                good reason to do so.
# axis -> Leave to 1 for simulated images. You *may* have to use axis = 0
#         if images are rotated by 90 deg. With axis = 1, we expect the cross-dispersion
#         to be along the Y axis.

def find_peak_blaze(file_flat, file_dark, percentiles = [85,97], axis = 1):

    # some output for the overeager
    print('Reading {0} (flat_flat) and {1} (dark_dark)'.format(file_flat, file_dark))

    # reading files
    image = fits.getdata(file_flat)
    dark = fits.getdata(file_dark)

    # find pixels that are more than +-10 median absolute deviations in the dark frame
    rms = np.nanmedian(abs(dark-np.nanmedian(dark)))
    # replace these pixels by NaNs
    dark[np.isfinite(dark) == False] = 1e9*rms
    image[np.abs(dark) > 10*rms] = np.nan

    # loop through percentile values to find midpoint.
    midpoints = []
    for percentile in range(percentiles[0],percentiles[1]+1):
        # some outputs as this is not instantaneous
        print('Computing blaze as {0}th percentile along columns'.format(percentile))
        # find the Nth percentile profile and filter it to avoid very hot/cold pixels
        # to bias the result
        blaze = medfilt(np.nanpercentile(image,percentile,axis = axis),31)
        # normalize by peak. We take the 99th percentile of the blaze to be its peak, again
        # to avoid spurious points to bias
        blaze /= np.nanpercentile(blaze,99)

        # finding the midpoint
        for cut in np.arange(.6,.9,.01):
            g = np.where(blaze>cut)[0]
            mid = (np.min(g)+np.max(g))/2

            # keep track of all midpoints
            midpoints.append(mid)

            plt.plot(mid,cut,'g.')
        plt.plot(blaze,label = str(percentile)+'th')

    if axis == 1:
        plt.xlabel('X pixel')
    if axis == 0:
        plt.xlabel('Y pixel')
    plt.ylabel('Blaze value & midpoints')
    plt.legend(loc = 0,fontsize = 'x-small')
    plt.show()

    # make it a numpy array for the median below
    midpoints = np.array(midpoints)

    print('********************************************************************************')
    print('')
    print('Midpoint of blaze = {0}+-{1}'.format(np.nanmedian(midpoints),np.nanstd(midpoints)))
    print('')
    print('********************************************************************************')

    return np.nanmedian(midpoints)