import numpy as np
from astropy.io import fits
from astropy.table import Table

# provide a min and a max position of the position of the center of the orders
# to be used in the mask creation
x1, x2 = 3200,3800

# offset between the region used by DRS and full frame. This is to account for the cropping of the image.
# For APERO, this effectively gets rid of the reference pixels
yoffset = 4
xoffset = 4


# APERO file use to create the mask
loco_file = 'NIRPS_GEN_ORDERDEF253_000F1T5_pp_loco_A.fits'

hdr = fits.getheader(loco_file)

if hdr['FIBER'] == 'A':
    width = 8

if hdr['FIBER'] == 'B':
    width = 8

tbl_c = Table.read(loco_file,hdu=1)  # centers
fit_c = np.zeros([len(tbl_c),len( tbl_c.keys()[1:])])


for i,key in  enumerate(tbl_c.keys()[1:]):
    fit_c[:,i] = tbl_c[key]

# mask to hold the expsure
mask = np.zeros([4096,4096])

ypixmap, xpixmap = np.indices([4096,4096])

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

# save mask
fits.writeto('mask.fits',mask, overwrite = True)