from astropy.io import fits
import numpy as np
from photutils import DAOStarFinder, morphology,data_properties
import matplotlib.pyplot as plt
import os
from astropy.table import Table

plt.ioff()


file_HC = 'NIRPS_GEN_WAVE018_0003.fits'

im = fits.getdata(file_HC)
std = np.nanpercentile(im,[16,84])
std = (std[1]-std[0])/2

# constructing a mask for the edges of the array
mask = np.ones_like(im)
mask[0:10, :] = 0
mask[-10:, :] = 0
mask[:, 0:10] = 0
mask[:, -10:] = 0

daofind = DAOStarFinder(fwhm=3.0, threshold=25. * std)

sources = daofind(im * mask)
keep = sources['sharpness']< 0.6
sources = sources[keep]

w = 10

box = np.zeros([w*8,w*8])

fig,ax = plt.subplots(nrows = 1,ncols = 1)
for xreg in range(4):
    for yreg in range(4):
        gg = (sources['xcentroid']//1024 == xreg)*(sources['ycentroid']//1024 == yreg)
        print(np.sum(gg))
        xx = np.array(np.round(sources['xcentroid'][gg]),dtype = int)
        yy = np.array(np.round(sources['ycentroid'][gg]),dtype = int)

        cube = np.zeros([len(xx),w*2,w*2])
        for j in range(len(xx)):
            cube[j,:,:]=im[yy[j]-w:yy[j]+w,xx[j]-w:xx[j]+w]
            cube[j,:,:] /=np.nansum(im[yy[j]-1:yy[j]+2,xx[j]-1:xx[j]+2])

        box[yreg*w*2:(yreg+1)*w*2,xreg*w*2:(xreg+1)*w*2] = np.nanmedian(cube,axis=0)

        p0 = np.nansum(box,axis = 0)
        p1 = np.nansum(box,axis = 1)

ax.imshow(box)
plt.show()