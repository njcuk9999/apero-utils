from astropy.io import fits
import glob
import numpy as np
import os
import sys
from scipy.signal import convolve2d
from astropy.table import Table

def fix_shift(image):
    """
    Read an image and fix a row of NaNS if required

    Status:
     status = 0 -> all good, no shit
     status = 1 -> normal shift in expected direction when we have a
                   controler issue, we apply correction
     status = -1 -> opposite direction of known shift pattern, we
                    apply correction
     status = -2 -> unknown shift, we need to look at the image

    :param image:
    :return:
    """

    # getting image size
    number_amps = 32
    width_image = image.shape[0]
    # get width of amplifier ribbon
    width_amp = int(width_image / number_amps)
    # look for a discontinuity between consecutive columns, this traces
    # reference vs science pixels
    gap = np.zeros(15)

    # we skip the first pixel as it is sometimes offset from the rest of
    # columns. To get imax = 0, we would need a shift of 3 pixels, which would
    # be very strange.
    for i in range(1,15):
        gap[i] = np.nanmedian(image[:, i] - image[:, i + 1])
    # If all is fine, we have 4 ref pixels and the 3rd difference shows a glitch
    imax = np.argmax(np.abs(gap))
    # assume no fix required by default
    status = 0
    # When images are bad, we have a value of 2 and know we need to
    # offset odd/even amplifier
    if imax == 2:
        print('\twe have a shift')

        for i in range(number_amps):
            offset = 1 - (i % 2) * 2

            start = i * width_amp
            end = i * width_amp + width_amp

            image[:, start:end] = np.roll(image[:, start:end], offset, axis=1)

        status = 1

    if imax == 3:
        print('\tall good, no shift')

    if imax == 4:
        print('\twe have a shift, but in the direction opposite to what '
              'we saw in 20200909 dataset')

        for i in range(number_amps):
            offset = (i % 2) * 2 - 1

            start = i * width_amp
            end = i * width_amp + width_amp

            image[:, start:end] = np.roll(image[:, start:end], offset, axis=1)

        status = -1

    if (imax != 2) * (imax != 3) * (imax != 4):
        print('\treally bad, rips the fabric of the Universe!')

        status = -2

    return image, status

def mk_isolated_nans(image):
    # input image is known to have its pixels in the right position
    # without shifts

    nans = np.isfinite(image) == 0
    ker = np.ones([3,3],dtype = float)
    neighbours = convolve2d(nans,ker,mode='same')

    isolated = nans*(neighbours==1)

    y,x = np.where(isolated)

    tbl = Table()
    tbl['x'] = x
    tbl['y'] = y
    tbl.write('isolated_nans.csv',overwrite = True)


def fix_shift2(image,threshold_valid = 0.9):
    """
    Read an image and fix a row of NaNS if required

    Status:
     status = 0 -> all good, no shit
     status = 1 -> normal shift in expected direction when we have a
                   controler issue, we apply correction
     status = -1 -> opposite direction of known shift pattern, we
                    apply correction
     status = -2 -> unknown shift, we need to look at the image

    threshold_valid
        is the fraction of NaNs that are in the expected position
        should be close to but not equal to 1. 0.9 seems OK, but this may require a
        confirmation with a large number of files

    :param image:

    :return:
    """

    tbl = Table.read('isolated_nans.csv')
    ref_nans = tbl['y'],tbl['x']

    #ref_nans = np.array(fits.getdata('isolated_nans.fits'), dtype = bool)

    # getting image size
    number_amps = 32
    width_image = image.shape[0]
    # get width of amplifier ribbon
    width_amp = int(width_image / number_amps)
    # look for a discontinuity between consecutive columns, this traces
    # reference vs science pixels

    frac_nans = np.nanmean(~np.isfinite(image[ref_nans]))

    print('no shift, frac of NaNs in right place: ',frac_nans)
    if frac_nans>threshold_valid:
        return image,0


    # try +1 and -1 roll
    image2 = np.array(image)

    # the usual correction
    for i in range(number_amps):
        offset = 1 - (i % 2) * 2

        start = i * width_amp
        end = i * width_amp + width_amp

        image2[:, start:end] = np.roll(image2[:, start:end], offset, axis=1)

    frac_nans = np.nanmean(~np.isfinite(image2[ref_nans]))

    print('usual shift, frac of NaNs in right place: ',frac_nans)
    if frac_nans>threshold_valid
        return image2,1

    image2 = np.array(image)
    for i in range(number_amps):
        offset = (i % 2) * 2 - 1

        start = i * width_amp
        end = i * width_amp + width_amp

        image[:, start:end] = np.roll(image[:, start:end], offset, axis=1)

    frac_nans = np.nanmean(~np.isfinite(image2[ref_nans]))
    print('opposite of usual shift, frac of NaNs in right place: ',frac_nans)
    if frac_nans>threshold_valid:
        return image2,-1

    return image, -2



if __name__ == '__main__':
    # get nightname from arguments
    args = sys.argv
    # if we have two arguments assume second argument is a night name
    if len(args) == 2:
        # get all files
        files = glob.glob(os.path.join('.', args[1], '*.fits'))
        # loop around files and see if we need to fix them
        for it, filename in enumerate(files):
            # print progress
            pargs = [os.path.basename(filename), it + 1, len(files)]
            print('Processing file {0} ({1}/{2})'.format(*pargs))
            # load data HDU
            hdu = fits.open(filename)
            # load the data into image 1
            image1 = np.array(hdu[1].data)
            # see / fix data
            image2, status = fix_shift(image1)
            # only save if status is not 0
            if status != 0:
                # push data back into hdu
                hdu[0].data = image2
                # write hdu to file
                hdu.writeto(filename, overwrite=True)
            # close hdu
            hdu.close()

