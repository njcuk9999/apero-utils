import numpy as np
from astropy.io import fits

# =============================================================================
# Define variables
# =============================================================================
# number of amplifiers in total
NAMPS = 32
# number of x pixels
NBXPIX = 4096
# the fits file to correct
FILENAME = '/cosmos99/nirps/apero-data/nirps_he_07276_online/raw/2022-12-06'


# =============================================================================
# Define functions
# =============================================================================
def estimate_sigma(tmp: np.ndarray) -> float:
    """
    Return a robust estimate of N sigma away from the mean

    :param tmp: np.array (1D) - the data to estimate N sigma of
    :param sigma: int, number of sigma away from mean (default is 1)

    :return: the sigma value
    """
    # get formal definition of N sigma
    sig1 = 0.6826894921370859
    # get the 1 sigma as a percentile
    p1 = (1 - (1 - sig1) / 2) * 100
    # work out the lower and upper percentiles for 1 sigma
    upper = np.nanpercentile(tmp, p1)
    lower = np.nanpercentile(tmp, 100 - p1)
    # return the mean of these two bounds
    return (upper - lower) / 2.0


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # load the fits image
    image = fits.getdata(FILENAME)
    # work out the width of the amplifiers
    width_amp = NBXPIX // NAMPS
    # find pixels in first read column and neighbouring valid column
    amp_pix = np.arange(NAMPS // 2)
    # calculate start and ending positions
    width_amp2 = width_amp * 2
    end = (width_amp * 2) - 1
    # get these as an array (for each amplifier)
    in_col = np.append(amp_pix * width_amp2,
                       amp_pix * width_amp2 + end)
    out_col = np.append(amp_pix * width_amp2 + 1,
                        amp_pix * width_amp2 + end - 1)
    # pixel-wise difference
    diff1 = image[:, in_col] - image[:, out_col]
    # find median bad column pattern
    bad_col_pattern = np.nanmedian(diff1, axis=1)
    # subtract the bad column pattern off the pixel-wise difference
    dd = np.array(diff1)
    for col in range(len(in_col)):
        dd[:, col] = dd[:, col] - bad_col_pattern
    # mask outliers not to be affected by residual flux
    bad = np.abs(dd) > 3 * estimate_sigma(dd)
    # remove these bad outliers from pixel-wise difference
    diff1[bad] = np.nan
    # median of diff given the bad column values
    bad_col_pattern = np.nanmedian(diff1, axis=1)
    # loop around columns in image and remove bad_col_pattern
    for col in in_col:
        image[:, col] = image[:, col] - bad_col_pattern
    # -------------------------------------------------------------------------
    # save the corrected image over the slope
    with fits.open(FILENAME) as hdul:
        # change the extension=1 image only
        hdul[1].data = image
        # write to new file
        hdul.writeto(FILENAME.replace('.fits', '_corrected.fits'),
                     overwrite=True)
