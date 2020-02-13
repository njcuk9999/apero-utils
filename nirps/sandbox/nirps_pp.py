import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
from scipy.ndimage import zoom

def med32(im):
    # input image MUST be 4096x4096
    # and have 32 amplifiers with a mirror
    # odd/even symmetry. We fold the image
    # into a 32x4096x128 cube where odd amplifiers
    # are flipped left/right. The median
    # amplifier structure is then re-padded
    # in an image that maps amp x-talk.
    # We expect the orders to be masked with NaNs
    # and the image to be high-passed so that only
    # high-frequency structures are left here.

    # cube to contain ordres in an easily
    # managed form
    cube = np.zeros([32, 4096, 128])
    for i in range(32):
        if (i % 2) == 0: # for left/right flipping
            i1 = i * 128
            i2 = i * 128 + 128
            sig = 1
        else:
            i1 = i * 128 + 127
            i2 = i * 128 -1
            sig = -1

        cube[i, :, :] = im[:, i1:i2:sig]
    # derive median amplifier structure
    med = np.nanmedian(cube, axis=0)


    # pad back onto the output image
    im2 = np.zeros_like(im)
    for i in range(32): # unflip
        if (i % 2) == 0: # for left/right flipping
            i1 = i * 128
            i2 = i * 128 + 128
            sig = 1
        else:
            i1 = i * 128 + 127
            i2 = i * 128 -1
            sig = -1

        im2[:, i1:i2:sig] = med
    return im2

def medbin(im,bx,by):
    # median-bin an image to a given size through
    # some funny np.reshape. To be used for low-pass
    # filterning of an image.
    sz = np.shape(im)
    out = np.nanmedian(np.nanmedian(im.reshape([bx,sz[0]//bx,by, sz[1]//bx]), axis=1), axis=2)
    return out

if not os.path.isfile('mask.fits'):
    # creation of the mask image
    im_flat = 'SIMU_NIRPS_HA_3f(flat_flat).fits'
    im = fits.getdata(im_flat)
    # find pixel that are more than 10 absolute deviations
    # from the image median
    im -= np.nanmedian(im)

    sig = np.nanmedian(np.abs(im))

    # generate a first estimate of the mask
    mask = im>10*sig
    # now apply a proper filtering of the image
    im2 = np.array(im)
    im2[mask] = np.nan

    ################################################
    # Same code as for an individual science frame #
    ################################################

    # median-bin and expand back to original size
    binsize = 32
    lowf = zoom(medbin(im2, binsize, binsize), 4096 // binsize)

    # subtract low-frequency from masked image
    im2 -= lowf
    # find the amplifier x-talk map
    xtalk = med32(im2)

    # subtract both low-frequency and x-talk from input image
    im -= (lowf + xtalk)

    # generate a better estimate of the mask
    mask = im>10*sig

    fits.writeto('mask.fits',np.array(mask,dtype = int), overwrite = True)


file = '20200207143909_HA_FP_FP.fits'
im = fits.getdata(file)

mask = np.array(fits.getdata('mask.fits'),dtype = bool)

im2 = np.array(im)
im2[mask] = np.nan

fits.writeto('test1.fits',im2,overwrite = True)

# we find the low level frequencies
# we bin in regions of 32x32 pixels. This CANNOT be
# smaller than the order footprint on the array
# as it would lead to a set of NaNs in the downsized
# image and chaos afterward
binsize = 32 # pixels

# median-bin and expand back to original size
lowf = zoom(medbin(im2,binsize,binsize),4096//binsize)

# subtract low-frequency from masked image
im2 -= lowf

# find the amplifier x-talk map
xtalk = med32(im2)

# subtract both low-frequency and x-talk from input image
im -= (lowf+xtalk)

fits.writeto('test2.fits',im,overwrite = True)