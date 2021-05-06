import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
from scipy.ndimage import zoom
import glob
from scipy.ndimage.morphology import binary_dilation
from scipy.ndimage import median_filter

def rot8(im,nrot):
    """
    Rotation of a 2d image with the 8 possible geometries. Rotation 0-3
    do not flip the image, 4-7 perform a flip

    nrot = 0 -> same as input
    nrot = 1 -> 90deg counter-clock-wise
    nrot = 2 -> 180deg
    nrot = 3 -> 90deg clock-wise
    nrot = 4 -> flip top-bottom
    nrot = 5 -> flip top-bottom and rotate 90 deg counter-clock-wise
    nrot = 6 -> flip top-bottom and rotate 180 deg
    nrot = 7 -> flip top-bottom and rotate 90 deg clock-wise
    nrot >=8 -> performs a modulo 8 anyway

    :param im: input image
    :param nrot: integer between 0 and 7
    :return: rotated and/or flipped image
    """
    nrot = int(nrot % 8)
    return np.rot90( im[::1-2*(nrot//4)],nrot % 4 )

def med32(im,doplot = False):
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
        cube[i, :, :] -= np.nanmedian(cube[i, :, :])
    # derive median amplifier structure

    med = np.nanmedian(cube, axis=0)
    cube2 = np.array(cube)

    diff_cuts = np.zeros(32)
    for i in range(32):
        cube2[i,:,:] -= med
        tmp = cube2[i,:,:]
        cut1 = np.nanpercentile(tmp,5)
        cut2 = np.nanpercentile(tmp,95)
        diff_cuts[i] = cut2-cut1
        tmp[tmp<cut1] = np.nan
        tmp[tmp>cut2] = np.nan
        cube2[i,:,:] = tmp
    bad_amps = diff_cuts>1.5*np.nanmedian(diff_cuts)

    cube[bad_amps,:,:] = np.nan
    cube2[bad_amps,:,:] = np.nan

    if doplot:
        plt.plot(np.arange(32),diff_cuts,'go', label = 'RMS of amplifier post-xtalk subtract')
        plt.plot(np.arange(32)[bad_amps],diff_cuts[bad_amps],'ro',label = 'removed amps')
        plt.ylabel('5th to 95th difference')
        plt.xlabel('nth amplifier')
        plt.legend()
        plt.show()
        cube[np.isfinite(cube2) ==  False] = np.nan

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

def get_mask(mask_file):
    # we need to pass a flat image to get the corresponding mask. The code creates a binary _mask file that is read
    # if it is present on disk

    # name of output mask
    outname = mask_file.split('.fits')[0]+'_mask.fits'

    if os.path.isfile(outname):
        return fits.getdata(outname)

    # creation of the mask image
    im = fits.getdata(mask_file)

    # find pixel that are more than 5 absolute deviations
    # from the image median
    im -= np.nanmedian(im)
    sig = np.nanmedian(np.abs(im))

    # generate a first estimate of the mask
    mask = im>5*sig

    # now apply a proper filtering of the image
    im2 = np.array(im)
    im2[mask] = np.nan

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
    mask = im>5*sig
    mask[np.isfinite(mask) == False] = np.nan
    mask = binary_dilation(mask, iterations = 4)

    # write binary mask file
    fits.writeto(outname,np.array(mask,dtype = int), overwrite = True)

    return np.array(mask,dtype = int)

def top_bottom(im, doplot = False):
    # remove the ramp in pixel values between the top and bottom of the array
    # corrects for amplified DC level offsets

    # map of pixel values within an amplified
    y, x = np.indices([4096, 128])
    y = y / 4096 # fraction of position on the amplifier

    med1 = np.zeros(32)
    med2 = np.zeros(32)
    for i in range(32):
        # median of bottom ref pixels
        med1[i] = np.nanmedian(im[0:4, i * 128:(i + 1) * 128])
        # median of top ref pixels
        med2[i] = np.nanmedian(im[-4:, i * 128:(i + 1) * 128])
        # subtraction of slope between top and bottom
        im[:, i * 128:(i + 1) * 128] -= (med1[i] + y * (med2[i] - med1[i]))

    if doplot:
        plt.plot(med1,'go',label = 'bottom ref pix')
        plt.plot(med2,'ro',label = 'top ref pix')
        plt.xlabel('Nth amplifier')
        plt.ylabel('median flux')
        plt.title('top/bottom ref pixel correction')
        plt.legend()
        plt.show()
    return im

def left_right(im, width = 15, doplot = False):
    # we correct for the left-right reference pixels
    # we take the median of the 8 left-right pixels to derive
    # a noise pattern that is 4096 pixels long. This noise pattern
    # is median-filtered with a kernel with a "width". Typical
    # values for this parameter are 10-20 pixels
    #
    reference_pixels_sides = im[:,[0,1,2,3,4092,4093,4094,4095]]
    for i in range(8):
        # remove DC level between pixels so that the median+axis=1 really
        # filters noise. Otherwise, the median would just select the pixel
        # that has the median DC level
        reference_pixels_sides[:,i] -= np.nanmedian(reference_pixels_sides[:,i])

    # median profile of the 8 pixels
    medprofile = np.nanmedian(reference_pixels_sides,axis=1)

    medprofile_filtered = median_filter(medprofile,width)
    # correlated noise replicated onto the output image format
    correlated_noise = np.repeat(medprofile_filtered, 4096).reshape(4096, 4096)

    if doplot:
        plt.plot(medprofile, label = 'median profile')
        plt.plot(medprofile_filtered,alpha = 0.5,label = 'after filtering by running median: {0} pixels'.format(width))
        plt.legend()
        plt.title('left-right pixel correction\ncorrelated noise correction')
        plt.show()

    return im - correlated_noise

def nirps_pp(files,mask_file = '',indir='',outdir='',doplot = False, force = False):
    # pre-processing of NIRPS images with only left/right and top/bottom pixels
    # if we pass the name of a flat_flat file as 'mask_file', then we also
    # correct for correlated noise between amps using dark pixels. This should
    # not be done on an LED image. You can set plot = True or force = True
    # if you want to see nice plots or force the overwrite of existing pre-process
    # files.

    # we provide or not a file to be used as a mask. Note that you need to point toward the FLAT+FLAT file,
    # not the mask created afterward. In other words, don't pass the _mask file here, but the input file
    # generated from it. If you pass '', you only have the top/bottom pixel correction that is applier
    if mask_file == '':
        xtalk_filter = False
    else:
        xtalk_filter = True

    # if more than one file has been passed
    if type(files) == str:
        files = glob.glob(files)

    for file in files:
        # if already a _pp file, skip
        if '_pp.' in file:
            print(file+' is a _pp file')
            continue

        # outname of files will be same as input but with a _pp suffix
        outname = '_pp.fits'.join(file.split('.fits'))
        outname = outdir+outname

        # if _pp file exists, skip unless you set force == True
        if os.path.isfile(outname) and (force == False):
            print('File : '+outname +' exists')
            continue
        else:
            print('We pre-process '+file)

            hdr = fits.getheader(indir+file)
            im = fits.getdata(indir+file)
            # before we even get started, we remove top/bottom ref pixels
            # to reduce DC level differences between ampliers
            for ite in range(3):
                im = top_bottom(im,doplot = doplot)
                im = left_right(im,width = 15,doplot = doplot)

            # we find the low level frequencies
            # we bin in regions of 32x32 pixels. This CANNOT be
            # smaller than the order footprint on the array
            # as it would lead to a set of NaNs in the downsized
            # image and chaos afterward
            binsize = 32 # pixels

            if xtalk_filter:
                # we fetch the mask
                mask = get_mask(mask_file)

                # create the matrix that will become the high-passed image
                im2 = np.array(im)
                im2[mask == 1] = np.nan

                # median-bin and expand back to original size
                tmp = medbin(im2,binsize,binsize)
                tmp[~np.isfinite(tmp)] = 0
                lowf = zoom(tmp,4096//binsize)

                # subtract low-frequency from masked image
                im2 -= lowf

                # remove correlated profile in Y axis
                yprofile1d = np.nanmedian(im2,axis=1)
                yprofile = np.repeat(yprofile1d, 4096).reshape(4096, 4096)
                im -= yprofile

                # first pixel of each amplifier
                first_col_x = (np.append(np.arange(16)*256,np.arange(16)*256-1+256))
                first_col_x.sort()

                # median-filter the first and last ref pixel, which trace the
                # behavior of the first-pixel effect
                amp0 = (median_filter(im[:,0],7) + median_filter(im[:,4095],7))/2

                for i in first_col_x:
                    im[:,i] -= amp0

            im = rot8(im,7) # rotate to match ESPRESSO order disposition

            fits.writeto(outname,im, hdr, overwrite = True)

    return []