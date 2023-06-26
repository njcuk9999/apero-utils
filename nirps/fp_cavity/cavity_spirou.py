import os
import glob
import numpy as np
from astropy.table import Table
from astropy.io import fits
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.constants import c
import george
from george import kernels
from george.modeling import Model
from etienne_tools import odd_ratio_mean

def bin_tbl(tbl):
    tbl2 = Table()

    udate = np.unique(tbl['WAVETIME'])

    tbl2['rjd'] = np.zeros(len(udate),dtype = float)
    tbl2['vrad'] = np.zeros(len(udate),dtype = float)
    tbl2['svrad'] = np.zeros(len(udate),dtype = float)
    tbl2['WAVETIME'] = np.zeros(len(udate),dtype = float)

    for iudate in range(len(udate)):
        idx = tbl['WAVETIME'] == udate[iudate]
        tbl2['vrad'][iudate],tbl2['svrad'][iudate] = odd_ratio_mean(tbl['vrad'][idx],tbl['svrad'][idx])
        tbl2['WAVETIME'][iudate] = udate[iudate]
        tbl2['rjd'][iudate] = np.mean(tbl['rjd'][idx])

    return tbl2

def sigma(y):
    return np.nanpercentile(np.abs(y-np.nanmedian(y)),68)


drift = Table.read('cavity_spirou.csv',format ='ascii')

drift = drift[drift['WAVETIME'] > 59360]
drift = drift[drift['WAVETIME'] < 59490]

drift['CAVITY'] = drift['WCAV_PED'] + drift['WCAV000']
drift['DV'] = (drift['CAVITY']/np.nanmedian(drift['CAVITY'])-1)*c.value
drift = drift[np.argsort(drift['WAVETIME'])]
drift = drift[drift['WAVETIME'] > 58868]

x = drift['WAVETIME']
y = drift['DV']
yerr = np.ones_like(y)*np.std(y)


ll = 30 # best - 20
length = ll**2
sig = sigma(y)**2#np.var(y)

kernel = sig * kernels.Matern32Kernel(length)

gp = george.GP(kernel)
gp.compute(x, yerr)

x_pred = np.linspace(np.min(drift['WAVETIME']), np.max(drift['WAVETIME']), 500)
pred, pred_var = gp.predict(y, x_pred, return_var=True)

pred2, pred_var2 = gp.predict(y, x, return_var=True)

plt.plot(x, y, 'k.',alpha = 0.5)
plt.plot(x_pred, pred, 'r-')
plt.show()

drift['DRIFT'] = pred2 - y

tbl = Table.read('lbl_GL699_GL699_bervzp.rdb')
tbl = tbl[tbl['WAVETIME'] > 59360]
tbl = tbl[tbl['WAVETIME'] < 59490]

tbl2 = bin_tbl(tbl)

tbl2['vrad'] -= odd_ratio_mean(tbl2['vrad'],tbl2['svrad'])[0]

tbl2['CORR'] = 0.0
for i in range(len(tbl2)):
    gg = np.abs(tbl2['WAVETIME'][i] == drift['WAVETIME'])
    if np.sum(gg) == 0:
        continue

    tbl2['CORR'][i] = -drift['DRIFT'][gg][0]


tbl2 = tbl2[tbl2['CORR'] != 0]

print(np.nanstd(tbl2['vrad'] ))
print(np.nanstd(tbl2['vrad'] - tbl2['CORR']))
print(np.nanstd(tbl2['vrad'] + tbl2['CORR']))

#gg = (tbl2['WAVETIME']>58970)*(tbl2['WAVETIME']<59020)


fig, ax = plt.subplots(nrows = 2, ncols =1, sharex = True)
ax[0].errorbar(tbl2['WAVETIME'], tbl2['vrad'], tbl2['svrad'], fmt='k.',alpha = 0.5)
ax[0].plot(tbl2['WAVETIME'], tbl2['CORR'],'ro')
ax[0].errorbar(tbl2['WAVETIME'], tbl2['vrad']-tbl2['CORR'], tbl2['svrad'], fmt='g.',alpha = 0.5)

ax[1].plot(drift['WAVETIME'], drift['DRIFT'],'k.')

plt.show()




tbl = Table.read('lbl_GL699_GL699_bervzp.rdb')
tbl['vrad'] -= odd_ratio_mean(tbl['vrad'],tbl['svrad'])[0]

x = tbl['rjd']
y = tbl['vrad']
yerr = tbl['svrad']


ll = 30 # best - 20
length = ll**2
sig = sigma(y)**2#np.var(y)

kernel = sig * kernels.Matern32Kernel(length)

gp = george.GP(kernel)
gp.compute(x, yerr)

x_pred = np.linspace(np.min(tbl['rjd']), np.max(tbl['rjd']), 500)
pred, pred_var = gp.predict(y, x_pred, return_var=True)

pred2, pred_var2 = gp.predict(y, x, return_var=True)

tbl['vrad'] -= pred2

tbl.write('test.rdb',format = 'rdb',overwrite = True)

plt.errorbar(x, y, yerr = yerr,fmt= 'k.',alpha = 0.5)
plt.errorbar(x, y-pred2, yerr = yerr,fmt= 'r.',alpha = 0.5)
plt.plot(x_pred, pred, 'r-')
plt.show()


#plt.plot(tbl['WAVETIME'], tbl['DV'], 'k.',alpha = 0.5)
#plt.plot(x_pred, pred, "k", lw=1.5, alpha=0.5)
#
#plt.show()