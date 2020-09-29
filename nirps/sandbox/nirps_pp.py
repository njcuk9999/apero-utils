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

        #fits.writeto('cube.fits',cube,overwrite =  True)
        #fits.writeto('cube2.fits',cube2,overwrite =  True)

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
    # we need to pass a flat image to get the corresponding mask.


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
    mask = im>5*sig
    mask[np.isfinite(mask) == False] = np.nan
    mask = binary_dilation(mask, iterations = 4)

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

def nirps_pp(files,mask_file = '', doplot = False, force = False):
    # pre-processing of NIRPS images with only left/right and top/bottom pixels
    # if we pass the name of a flat_flat file as 'mask_file', then we also
    # correct for correlated noise between amps using dark pixels. This should
    # not be done on an LED image. You can set plot = True or force = True
    # if you want to see nice plots or force the overwrite of existing pre-process
    # files.

    if mask_file == '':
        xtalk_filter = False
    else:
        xtalk_filter = True

    ref_hdr = fits.getheader('ref_hdr.fits')

    if type(files) == str:
        files = glob.glob(files)

    for file in files:
        outname = '_pp.fits'.join(file.split('.fits'))

        if '_pp.' in file:
            print(file+' is a _pp file')
            continue

        if os.path.isfile(outname) and (force == False):
            print('File : '+outname +' exists')
            continue
        else:
            print('We pre-process '+file)

            hdr = fits.getheader(file)


            im = fits.getdata(file)
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

            if 'MJDEND' not in hdr:
                hdr['MJDEND'] = 0.00
                hdr['EXPTIME'] = 5.57*len(hdr['INTT*'])

            hdr['MJDMID'] = hdr['MJDEND'] - hdr['EXPTIME']/2.0/86400.0

            hdr['INF1000'] = file
            DPRTYPES = ['DARK_DARK','DARK_FP','FLAT_FLAT','DARK_FLAT',
                        'FLAT_DARK','HC_FP','FP_HC','FP_FP','OBJ_DARK',
                        'OBJ_FP','HC_DARK','DARK_HC','HC_HC','DARK_FP', 'FP_DARK','LED_LED']

            if 'STAR_DARK' in file:
                hdr['DPRTYPE'] = 'OBJ_DARK'

            if 'STAR_FP' in file:
                hdr['DPRTYPE'] = 'OBJ_FP'

            hdr['DPRTYPE'] = 'DARK_DARK'

            for DPRTYPE in DPRTYPES:
                if DPRTYPE in file:
                    if DPRTYPE == 'DARK_DARK':
                        hdr['DPRTYPE'] = 'DARK_DARK_TEL'
                    elif DPRTYPE == 'HC_HC':
                        hdr['DPRTYPE'] = 'HCONE_HCONE'
                    elif DPRTYPE == 'FP_HC':
                        hdr['DPRTYPE'] = 'FP_HCONE'
                    elif DPRTYPE == 'HC_FP':
                        hdr['DPRTYPE'] = 'HCONE_FP'
                    elif DPRTYPE == 'DARK_HC':
                        hdr['DPRTYPE'] = 'DARK_HCONE'
                    elif DPRTYPE == 'HC_DARK':
                        hdr['DPRTYPE'] = 'HCONE_DARK'
                    elif DPRTYPE == 'FP_DARK':
                        hdr['DPRTYPE'] = 'FP_DARK'
                    elif DPRTYPE == 'DARK_FP':
                        hdr['DPRTYPE'] = 'DARK_FP'
                    elif DPRTYPE == 'LED_LED':
                        hdr['DPRTYPE'] = 'LED_LED'


            if 'DPRTYPE' not in hdr:
                print('error, with DPRTYPE for ',file)
                return

            if 'OBJECT' not in hdr:
                hdr['OBJECT'] = 'none'

            if 'RDNOISE' not in hdr:
                hdr['RDNOISE']= 10.0,'rdnoise *not* provided, added by _pp'

            if 'GAIN' not in hdr:
                hdr['GAIN']= 1.000,'gain *not* provided, added by _pp'

            if 'SATURATE' not in hdr:
                hdr['SATURATE']= 60000,'saturate *not* provided, added by _pp'

            if 'PVERSION' not in hdr:
                hdr['PVERSION'] = 'NIRPS_SIMU_PP'

            if 'OBSTYPE' not in hdr:
                if hdr['DPRTYPE'][0:4] == 'FLAT':
                    hdr['OBSTYPE'] = 'FLAT'

                if hdr['DPRTYPE'][0:4] == 'DARK':
                    hdr['OBSTYPE'] = 'DARK'

                if hdr['DPRTYPE'][0:2] == 'FP':
                    hdr['OBSTYPE'] = 'ALIGN'

                if hdr['DPRTYPE'][0:2] == 'HC':
                    hdr['OBSTYPE'] = 'COMPARISON'

                if hdr['DPRTYPE'][0:3] == 'OBJ':
                    hdr['OBSTYPE'] = 'OBJECT'

                if hdr['DPRTYPE'][0:3] == 'LED':
                    hdr['OBSTYPE'] = 'LED'

            if hdr['DPRTYPE'][0:3] == 'OBJ':
                hdr['TRG_TYPE'] = 'TARGET'
            else:
                hdr['TRG_TYPE'] = ''

            necessary_kwrd = ['OBSTYPE','TRG_TYPE','OBJECT','OBJRA','OBJDEC','OBJECT','OBJEQUIN','OBJRAPM','OBJDECPM','AIRMASS','RELHUMID','OBJTEMP','GAIA_ID','OBJPLX','OBSRV','GAIN','RDNOISE','FRMTIME','EXPTIME','PI_NAME','CMPLTEXP','NEXP','MJDATE','MJDEND','SBCREF_P','SBCCAS_P','SBCALI_P','SBCDEN_P','DATE-OBS','UTC-OBS','SATURATE','TEMPERAT','SB_POL_T']

            missing = False
            for key in necessary_kwrd:
                if key not in hdr:
                    print('missing keyword : {0}'.format(key))
                    missing = True

                    if key in ref_hdr:
                        hdr[key] = ref_hdr[key]

            if hdr['OBSTYPE'] == 'LED':
                # we cannot correct the capacitive coupling between amplifiers if we have an LED image
                xtalk_filter = False


            if xtalk_filter:
                mask = get_mask(mask_file)

                im2 = np.array(im)
                im2[mask] = np.nan

                # median-bin and expand back to original size
                lowf = zoom(medbin(im2,binsize,binsize),4096//binsize)
                # subtract low-frequency from masked image

                # find the amplifier x-talk map
                xtalk = med32(im2-lowf,doplot = doplot)
                im2 -= xtalk
                # subtract both low-frequency and x-talk from input image
                im -= xtalk

            fits.writeto(outname,im, hdr, overwrite = True)

            return

            """
            # rotates the image so that it matches the order geometry of SPIRou and HARPS
            # redder orders at the bottom and redder wavelength within each order on the left

            # NIRPS = 5
            # SPIROU = 3
            im = rot8(im,5)



            b = fits.getdata(file,ext = 2)
            errslope = fits.getdata(file,ext = 3)
            n = fits.getdata(file,ext = 4)

            b = rot8(b,5)
            errslope = rot8(errslope,5)
            n = rot8(n,5)


            hdu1 = fits.PrimaryHDU()
            hdu1.header = hdr
            hdu1.header['NEXTEND'] = 4
            hdu2 = fits.ImageHDU(im)
            hdu2.header['UNITS'] = ('ADU/S', 'Slope of fit, flux vs time')
            hdu2.header['EXTNAME'] = ('slope', 'Slope of fit, flux vs time')

            hdu3 = fits.ImageHDU(b)
            hdu3.header['UNITS'] = ('ADU', 'Intercept of the pixel/time fit.')
            hdu3.header['EXTNAME'] = ('intercept', 'Intercept of the pixel/time fit.')

            hdu4 = fits.ImageHDU(errslope)
            hdu4.header['UNITS'] = ('ADU/S', 'Formal error on slope fit')
            hdu4.header['EXTNAME'] = ('errslope', 'Formal error on slope fit')

            hdu5 = fits.ImageHDU(n)
            hdu5.header['UNITS'] = ('Nimages', 'N readouts below saturation')
            hdu5.header['EXTNAME'] = ('count', 'N readouts below saturation')

            new_hdul = fits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])

            # just to avoid an error message with writeto
            if os.path.isfile(outname):
                print('file : ' + outname + ' exists, we are overwriting it')
                os.system('rm ' + outname + '')

            new_hdul.writeto(outname, overwrite=True)
            """


    return []