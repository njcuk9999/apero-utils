import numpy as np
from astropy.io import fits
from astropy.table import Table

def loco2mask(loco_file, xrange = [3200,3800],xyoffset = [4,4], width = 15,outname = 'mask.fits'):
    # Parameters:
    #
    # input file from either ESPRESSO-like or APERO pipeline
    # the 1st extension should be a fits table with the localisation coefficients
    # Does not use the header directly, but propagates it to the mask
    #
    # xrange : x position for the orders to be used in the mask. Should be a list of two integers
    #          xrange = [3200,3800] is reasonable to get the H band
    #
    # xyoffsets : Transforms the science pixels to mask pixels. In both pipelines, we use a cropped 4088x4088
    #             image that removes the reference pixels while the mask is used on the full 4096x4096 image.
    #             For now, the offset is therefore 4 pixels on both axis.
    #
    # widths : widht of the mask relative to the polynomial. 15 pixels for the ESPRESSO pipeline.
    # outname : name of the output mask
    #

    # provide a min and a max position of the position of the center of the orders
    # to be used in the mask creation
    x1, x2 = int(xrange[0]), int(xrange[1])

    # offset between the region used by DRS and full frame. This is to account for the cropping of the image.
    # For both APERO and ESPRESSO codes, this effectively gets rid of the reference pixels
    xoffset, yoffset =  int(xyoffset[0]), int(xyoffset[1])

    print('\n\tReading {0}'.format(loco_file))
    # hdr that gets propagated to the mask
    hdr = fits.getheader(loco_file)

    # we check the first extension header, if it is an ESPRESSO loco, this will an image, if it is APERO, then
    # its the first extension that is the table
    h1 = fits.getheader(loco_file, ext=1)

    if 'IMAGE' in h1['XTENSION']:
        print('\tThis is an *ESPRESSO-like* loco file')
        hdu = 2
    else:
        print('\tThis is an *APERO* loco file')
        hdu = 1

    # read table of coefficiens
    tbl_c = Table.read(loco_file,hdu=hdu)  # centers

    # there may be other keys, just keep keys that include 'COEFF' in the column name
    print('\tConstructing COEFF table')
    keys = np.array(tbl_c.keys())
    keep_key  = np.zeros_like(keys,dtype = bool)
    for i in range(len(keep_key)):
        if 'COEFF' in tbl_c.keys()[i].upper():
            keep_key[i] = True

    coeff_keys = keys[keep_key]

    # map of polynomial coefficients
    fit_c = np.zeros([len(tbl_c),len(coeff_keys)])

    for i,key in  enumerate(coeff_keys):
        fit_c[:,i] = tbl_c[key]

    # mask to hold the expsure
    mask = np.zeros([4096,4096])

    ypixmap, xpixmap = np.indices([4096,4096])

    print('\tMapping coefficients to pixel space')
    ypix = np.arange(4096)
    for i in range(len(fit_c)):
        # get x position of order
        xpix = np.polyval(fit_c[i][::-1],ypix)

        # see if it falls in the authorized pixel range
        if (np.max(xpix)>x1)*(np.max(xpix)<x2):
            # find pixels that are in the range expected
            xpix = np.repeat(xpix+xoffset, 4096).reshape(4096, 4096)

            # add this to the mask
            mask += np.abs(xpix-xpixmap)<width/2.0

    # save mask, be sure that the filename includes .fits as an extension
    if '.fits' not in outname:
        outname = outname+'.fits'

    print('\tWe write "{0}"\n'.format(outname))
    fits.writeto(outname,mask, hdr, overwrite = True)

# Sample demo lines, to be deleted

# APERO file use to create the mask
loco_file = 'NIRPS_GEN_ORDERDEF253_000F1T5_pp_loco_A.fits'
loco2mask(loco_file,outname = 'mask_APERO.fits')

loco_file = 'NIRPS_HA_ORDER_TABLE_A.fits'
loco2mask(loco_file,outname = 'mask_ESPRESSO.fits')

