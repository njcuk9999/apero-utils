import glob
from astropy.io import fits
import os
from tqdm import tqdm
import etienne_tools as et
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from compilblrv import *
# DO NOT TOUCH, this is 100% test zone and not to be used by anyone but me!
#
#                                     grrrrrr
#
obj_sci = 'GL699'
obj_template = 'GL699'
doplot = True
force = True
common_weights = True

tbl = compilblrv(obj_sci, obj_template = obj_template, doplot = doplot, force = force, common_weights = common_weights)

tbl2 = Table()
udates =np.unique(tbl['DATE-OBS'])
tbl2['MJDATE'] = np.zeros(len(udates),dtype = float)
tbl2['RV'] = np.zeros(len(udates),dtype = float)
tbl2['ERR'] = np.zeros(len(udates),dtype = float)

for i in range(len(udates)):
    g = (tbl['DATE-OBS'] == udates[i])

    err = tbl['per_epoch_err'][g]
    rv = tbl['per_epoch_mean'][g]

    tbl2['RV'][i] = np.sum( rv/err**2 )/np.sum(1/err**2)
    tbl2['ERR'][i]  = np.sqrt(1/np.sum(1/err**2))
    tbl2['MJDATE'][i] = np.mean(tbl['MJDATE'][g])


if obj_sci == 'GL436':
    tbl2 = tbl2[(tbl2['MJDATE']>59000)*(tbl2['ERR'] < 5)]
    plt.errorbar((tbl2['MJDATE'] / 2.6439) % 1,tbl2['RV'] - np.nanmedian(tbl2['RV']),fmt='.g', yerr=tbl2['ERR'])
    plt.xlabel('Orbital phase')
else:
    #tbl2 = tbl2[(tbl2['MJDATE']>59000)*(tbl2['ERR'] < 5)]
    plt.errorbar(tbl2['MJDATE'],tbl2['RV'] - np.nanmedian(tbl2['RV']),fmt='.g', yerr=tbl2['ERR'])
    plt.xlabel('Date')


plt.ylabel('Velocity [m/s]')
plt.title(obj_sci)
plt.show()

print( np.mean(tbl2['ERR']),np.nanstd(tbl2['RV']))
