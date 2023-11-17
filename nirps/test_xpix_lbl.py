from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
from etienne_tools import odd_ratio_mean, sigma
from astropy.io import fits



#files = glob.glob('/Volumes/courlan/lbl_NIRPS_HE/lblrv/PROXIMA_PCAx_PROXIMA_PCAx/*fits')
files = glob.glob('/Volumes/courlan/lbl_SPIROU/lblrv/GL410-PCA03_GL410-PCA03/*fits')

rjd = np.zeros(len(files))
rv1 = np.zeros(len(files))
rv2 = np.zeros(len(files))
rv3 = np.zeros(len(files))

sig_rv1 = np.zeros(len(files))
sig_rv2 = np.zeros(len(files))
sig_rv3 = np.zeros(len(files))

berv = np.zeros(len(files))

xminmax = 1000,3000
for i in range(len(files)):
    print(i, len(files))
    tbl = Table.read(files[i])
    h = fits.getheader(files[i])

    #rjd[i] = h['MJD-OBS']
    rjd[i] = h['MJDATE']

    berv[i] = h['BERV']
    xpix = np.array(tbl['XPIX'], dtype = float)
    dv = np.array(tbl['dv'], dtype = float)
    sdv = np.array(tbl['sdv'], dtype = float)

    rv1[i],sig_rv1[i] = odd_ratio_mean(dv,sdv)

    g = np.where((xpix<  xminmax[0]))[0]
    rv2[i], sig_rv2[i] = odd_ratio_mean(dv[g],sdv[g])

    g = np.where((xpix>  xminmax[1]))[0]
    rv3[i], sig_rv3[i] = odd_ratio_mean(dv[g],sdv[g])

fig, ax = plt.subplots(nrows = 3, ncols = 1, sharex = True)
ax[0].errorbar(rjd,rv1,yerr=sig_rv1,fmt='.',label = 'no masking')
ax[0].errorbar(rjd,rv2,yerr=sig_rv2,fmt='.', label = 'masking ')
ax[0].legend()
ax[1].errorbar(rjd,rv2-rv3,yerr=np.sqrt(sig_rv2**2+sig_rv3**2),fmt='o')
ax[2].plot(rjd,berv,'.')
plt.show()