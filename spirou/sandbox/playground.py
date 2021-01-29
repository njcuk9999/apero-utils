import glob
from astropy.io import fits
import os
from tqdm import tqdm
import etienne_tools as et
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from compilblrv import *
from astropy.time import Time


if False:
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


if False:
    files = np.array(glob.glob('tbllblrv_G*.csv'))
    #files = np.append(files, 'tbllblrv_FP_FP.csv')



    for file in files:

        if 'GL876' in file:
            continue

        tbl = Table.read(file)
        if 'FP' in file:
            obj = 'FP'
        else:
            obj = tbl['OBJECT'][0]
        #g = (tbl['MJDATE']>0)*(tbl['MJDATE']<1e9)*(tbl['per_epoch_err']<10)
        g = (tbl['MJDATE']>58950)*(tbl['MJDATE']<59137)*(tbl['per_epoch_err']<5)
        if np.nansum(g) <5:
            continue

        tbl = tbl[g]
        tbl['per_epoch_mean'] -= np.nanmedian(tbl['per_epoch_mean'])

        #if 'GL846' in file:
        #    tbl['per_epoch_mean'] +=50

        if obj == 'FP':
            tbl['per_epoch_mean_H'] -= np.nanmedian(tbl['per_epoch_mean_H'])
            plt.errorbar(tbl['MJDATE'], tbl['per_epoch_mean_H'], fmt='k.',yerr=tbl['per_epoch_err_H'],label =obj)
        else:
            plt.errorbar(tbl['MJDATE'], tbl['per_epoch_mean'], fmt='.',alpha = 0.5,yerr=tbl['per_epoch_err'],label =obj)


    plt.xlabel('MJDATE')
    plt.ylabel('dv [m/s]')
    plt.legend()
    plt.show()

# rossiter 58648.526577
if False:
    tbl = Table.read('tbllblrv_TRAPPIST-1_TRAPPIST-1.csv')
    obj = tbl['OBJECT'][0]
    tbl['per_epoch_mean'] -= np.nanmedian(tbl['per_epoch_mean'])
    tbl['per_epoch_mean_H'] -= np.nanmedian(tbl['per_epoch_mean_H'])


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

    t3 = Time(tbl['MJDATE'], format='mjd')
    plt.plot_date(t3.plot_date, tbl['per_epoch_mean_H'], 'r.', alpha=0.5)

    #plt.errorbar(tbl['MJDATE'], tbl['per_epoch_mean_H'], fmt='g.', alpha=0.5,
    #             yerr=tbl['per_epoch_err_H'], label=obj)

    #plt.errorbar(tbl2['MJDATE'],tbl2['RV'],yerr = tbl2['ERR'], fmt = 'r.',)
    plt.xlabel('MJDATE')
    plt.ylabel('dv [m/s]')
    plt.legend()
    plt.show()


if False:
    force = False
    tbl1 =   compilblrv('FP',common_weights = True,force = force)
    tbl2 =   compilblrv('FP',common_weights = False,force = force)

    plt.plot(tbl1['MJDATE'],tbl1['per_epoch_mean_H_2044-4088'],'r.')
    plt.plot(tbl2['MJDATE'],tbl2['per_epoch_mean_H_2044-4088'],'g.')
    plt.show()


if False:
    tbl = Table.read('tbllblrv_FP_FP.csv')

    fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True)

    ax[0].errorbar(tbl['MJDATE'], tbl['per_epoch_mean'], fmt='g.', alpha=0.5,
                 yerr=tbl['per_epoch_err'])
    ax[1].errorbar(tbl['MJDATE'], tbl['per_epoch_DDV'], fmt='g.', alpha=0.5,
                 yerr=tbl['per_epoch_DDVRMS'])
    ax[0].set(xlabel = 'MJDATE', ylabel = 'DV [m/s]',title = 'FP velocity')
    ax[1].set(xlabel = 'MJDATE', ylabel = 'DDV [(m/s)^2]',title= '2nd derivative')
    plt.show()