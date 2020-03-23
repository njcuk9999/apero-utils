import numpy as np
from scipy.signal import medfilt,convolve2d
from astropy.table import Table
from astropy.io import fits

# NIRPS file
file = '20200218111102_HA_DARK_DARK_pp.fits'
outname = 'NIRPS_hotpix.csv'

# SPIRou file
#file = '2401543d.fits'
#outname = 'SPIRou_hotpix.csv'

im = fits.getdata(file)
# set NaNs and inf to zero. NaN pixels will not be flagged as
# hotpix
im[~np.isfinite(im)] = 0

# subtract a DC offset of the image level
im -= np.nanmedian(im)
# express things normalized in terms of sigma
im / np.nanpercentile(np.abs(im),68.2689)

# a hot pixel is a point that is >10-sigma (positive)
# and that has a 5x5 median around it that is within
# +-1 sigma; it is well-behaved and not surrounded by
# bad pixels

medim = medfilt(im,[5,5])

hotpix = np.array(np.array( (np.abs(medim) < 1.0) & (im>10.0)  ),dtype = float)

# find if hot pixels are alone in a 5x5 box
neighbours = convolve2d(hotpix,np.ones([5,5],dtype = float),mode = 'same')

# after the convolution, isolated (within 5x5) hotpixels have
# neighbours = 1
hotpix = hotpix*(neighbours ==1)

# find x/y position of hot pixels

y,x = np.where(hotpix)

# output to a file
tbl = Table()
tbl['nsig'] = im[y,x] #
tbl['xpix'] = x
tbl['ypix'] = y

# write table as csv file
tbl.write(outname,overwrite = True)
# mask to compare with the image. You probably don't want to keep
# that one after debugging
fits.writeto('mask_'+file,hotpix,overwrite = True)