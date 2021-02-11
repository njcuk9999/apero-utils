import numpy as np
from astropy.io import fits
from scipy.ndimage.morphology import binary_dilation
import glob
import os

def xpand_mask(mask1,mask2):
    # find all pixels within mask2 that include a mask1 pixel
    increment = 1
    sum_prev = 0
    while increment != 0:
        mask1 = np.array((mask2) * (binary_dilation(mask1)))
        increment = np.nansum(mask1) - sum_prev
        sum_prev = np.nansum(mask1)
    return mask1

def read_smart_cosmic(file, header = False, variance_cuts = 100, intercept_cuts = 50):

    ron = 30.0 # super-pessimistic noise estimate. Includes uncorrected common noise

    # read header
    hdr = fits.getheader(file)
    im = fits.getdata(file, ext=1) # flux
    intercept = fits.getdata(file, ext=2) # intercept
    errslope = fits.getdata(file, ext=3) # error on slope
    inttime =  fits.getdata(file, ext=4)*hdr['FRMTIME'] # pixel exposure time

    im2 =im*inttime # flux expressed in ADUs, not ADU/s
    errslope*=inttime # same for slope

    variance = errslope**2 # express excursions as variance so we can subtract things
    for i in range(32):
        box = variance[[0,1,2,3,4096-4,4096-3,4096-2,4096-1],i*128:i*128+128]
        variance[:,i*128:i*128+128] -= np.nanmedian(box) # subtract median per-amplifier variance

    im2[im2<0] = 0 # cannot be smaller than zero
    expected = im2+ron**2
    nsig2 = variance/expected
    # number of sigma away from bulk of expected-to-observed variance
    nsig2 /= np.nanpercentile(np.abs(nsig2),68) # 1-sigma


    mask_slope_variance = nsig2>variance_cuts#[0]

    # mask of where variance is bad
    # set to NaN bad pixels
    im[mask_slope_variance] = np.nan

    # remove median per-column intercept
    for i in range(4096):
        intercept[:,i] -= np.nanmedian(intercept[:,i])
    # remove per-region intercept
    for i in range(64):
        for j in range(64):
            intercept[i*64:i*64+64,j*64:j*64+64] -= np.nanmedian(intercept[i*64:i*64+64,j*64:j*64+64])

    # normalize to 1-sigme
    intercept/=np.nanpercentile(np.abs(intercept),68)
    # express as variance
    nsig2 = intercept**2

    mask_intercept_deviation = nsig2>intercept_cuts#[0]
    #mask2 = nsig2>intercept_cuts[1]

    # find spurious pixels
    #mask_intercept_deviation = xpand_mask(mask1,mask2)
    im[mask_intercept_deviation] = np.nan

    #fits.writeto('nsig2_intercept.fits',nsig2, overwrite = True)

    if header == False:
        return im
    else:
        hdr = fits.getheader(file)
        # adding some keywords to quantify the number of bad pixels
        hdr['NBADINTE'] = np.nansum(mask_intercept_deviation), 'Number of bad pix to intercept err'
        hdr['NBADSLOP'] = np.nansum(mask_slope_variance), 'Number of bad pix to slope err'
        hdr['NBADBOTH'] = np.nansum(mask_slope_variance*mask_intercept_deviation), 'Number of bad pix to both slope and interc.'

        return im,hdr


files = glob.glob('*.fits')
for file in files:
    if 'corr' in file:
        continue
    print(file)
    outname = '_corr.'.join(file.split('.'))
    if os.path.isfile(outname):
        print('File {0} exists'.format(outname))
        continue

    im,hdr  = read_smart_cosmic(file,header = True)
    fits.writeto(outname,im,hdr,overwrite = True)
