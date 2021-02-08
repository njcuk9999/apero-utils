from astropy.io import fits
import numpy as np
from photutils import DAOStarFinder, morphology,data_properties
import matplotlib.pyplot as plt
import os
from astropy.table import Table
import etienne_tools as et

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

fig,ax = plt.subplots(nrows = 1,ncols = 2, figsize = [16,8])
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

        tmp = np.nanmedian(cube,axis=0)
        box[yreg*w*2:(yreg+1)*w*2,xreg*w*2:(xreg+1)*w*2] = tmp

        p0 = np.nansum(tmp,axis = 0)
        p1 = np.nansum(tmp,axis = 1)

        #cen, ew, amp, zp, slope
        fit = et.fit_gauss(np.arange(len(p0)),p0,[w,2,np.max(p0),0,0])
        fw0 = fit[1]* np.sqrt(np.log(2)*2)*2

        fit = et.fit_gauss(np.arange(len(p1)),p1,[w,2,np.max(p1),0,0])
        fw1 = fit[1]* np.sqrt(np.log(2)*2)*2

        print(fw0,fw1)

        ax[1].text(xreg*1024+512,yreg*1024+512,'{0:.2f}pix\n{1:.2f}pix'.format(fw1,fw0))

ax[0].imshow(box,origin = 'lower')
ax[0].set(title = 'PSF map')
ax[1].set(xlim = [0,4088],ylim = [0,4088])
plt.show()