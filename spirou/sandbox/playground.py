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
if False:
    tbl = compilblrv('FP',force = True)
    tbl['SED'] = np.zeros(len(tbl), dtype = bool)

    for i in range(len(tbl)):
        if 'sed' in tbl['LOCAL_FILE_NAME'][i]:
            tbl['SED'][i] = True

    #g1 = tbl['CDBWAVE'] == '257146F1T2c_pp_e2dsff_C_wave_night_C.fits'
    #g2 = ~g1

    key0 = 'per_epoch_mean_H_1532-2556'
    #key1 = 'per_epoch_mean_Y_0-2044'
    #key2 = 'per_epoch_mean_Y_1532-2556'
    #key3 = 'per_epoch_mean_Y_2044-4088'
    #key4 = 'per_epoch_mean_H_0-2044'
    #key5 = 'per_epoch_mean_H_2044-4088'

    #plt.errorbar(tbl['per_epoch_mean_H'][g1],tbl[key1][g1] - tbl[key2][g1], fmt='.', alpha=0.5,
    #          yerr=tbl['per_epoch_err_H'][g1],color = 'red')
    #plt.errorbar(tbl['per_epoch_mean_H'][g2], tbl[key1][g2] - tbl[key2][g2], fmt='.', alpha=0.5,
    #          yerr=tbl['per_epoch_err_H'][g2],color = 'blue')

    #plt.plot(tbl['MJDATE'], tbl[key1], '.', alpha=0.5,color = 'red',label = key1)#,
    #plt.plot(tbl['MJDATE'], tbl[key2], '.', alpha=0.5, color='green',label = key2)  # ,
    #plt.plot(tbl['MJDATE'], tbl[key3], '.', alpha=0.5, color='blue',label = key3)  # ,
    #plt.plot(tbl['MJDATE'], tbl[key1][tbl['SED']], 'o', alpha=0.5, color='black',label = 'Y band')  # ,
    plt.plot(tbl['MJDATE'][tbl['SED']], tbl[key0][tbl['SED']], 'o', alpha=0.5, color='black',label = 'SED')  # ,
    plt.plot(tbl['MJDATE'][~tbl['SED']], tbl[key0][~tbl['SED']], 'o', alpha=0.5, color='red',label = 'not SED')  # ,
    #plt.plot(tbl['MJDATE'], tbl[key4] - tbl[key5], 'o', alpha=0.5, color='red',label = 'H band')  # ,
    #          yerr=tbl['per_epoch_err_H'][g1],color = 'red')
    #plt.plot(tbl['MJDATE'][g2], tbl[key1][g2], fmt='.', alpha=0.5)#,
    plt.legend()
    plt.xlabel('MJDATE')
    plt.ylabel('dv in H [m/s]')
    plt.show()
if False:
    tbl1 = Table.read('2571513o_pp_e2dsff_C_FP_FP_lbl.fits')
    tbl2 = Table.read('2571514o_pp_e2dsff_C_FP_FP_lbl.fits')

    g = (tbl1['ORDER'] > 35)*(tbl1['ORDER'] < 38)*(tbl1['DVRMS']>1)*(tbl1['DVRMS']<100)
    tbl1 = tbl1[g]
    tbl2 = tbl2[g]


    fig,ax = plt.subplots(nrows = 2, ncols = 1,sharex = True)
    ax[0].plot(tbl1['XPIX'],tbl1['RV'] - tbl2['RV'],'.',alpha = 0.5)#, yerr =tbl1['DVRMS']  )
    ax[0].plot([0,4096],[0,0],color = 'red')
    ax[0].set(xlabel = 'position in pix',ylabel = 'dv [m/s]')


    ax[1].plot(tbl1['XPIX'],tbl1['DDV'] - tbl2['DDV'],'.',alpha = 0.5)#, yerr =tbl1['DVRMS']  )
    ax[1].plot([0,4096],[0,0],color = 'red')
    ax[1].set(xlabel = 'position in pix',ylabel = 'ddv ')

    plt.show()

    plt.xlabel('MJDATE')
    plt.ylabel('dv in H left/right [m/s]')
    plt.show()
if False:
    tbl1 = Table.read('2571513o_pp_e2dsff_C_FP_FP_lbl.fits')
    tbl2 = Table.read('2571514o_pp_e2dsff_C_FP_FP_lbl.fits')

    g = (tbl1['ORDER'] > 40)*(tbl1['ORDER'] < 45)*(tbl1['DVRMS']>1)*(tbl1['DVRMS']<100)
    tbl1 = tbl1[g]
    tbl2 = tbl2[g]


    fig,ax = plt.subplots(nrows = 2, ncols = 1,sharex = True)
    ax[0].plot(tbl1['XPIX'],tbl1['RV'] - tbl2['RV'],'.',alpha = 0.5)#, yerr =tbl1['DVRMS']  )
    ax[0].plot([0,4096],[0,0],color = 'red')
    ax[0].set(xlabel = 'position in pix',ylabel = 'dv [m/s]')


    ax[1].plot(tbl1['XPIX'],tbl1['DDV'] - tbl2['DDV'],'.',alpha = 0.5)#, yerr =tbl1['DVRMS']  )
    ax[1].plot([0,4096],[0,0],color = 'red')
    ax[1].set(xlabel = 'position in pix',ylabel = 'ddv ')

    plt.show()


if True:
    obj_sci = 'TRAPPIST-1'
    obj_template = 'TRAPPIST-1'
    doplot = False
    force = True


    tbl1 = compilblrv(obj_sci, obj_template = obj_template, doplot = doplot, force = force, common_weights = False,
                   get_cumul_plot = False)
    #bl2 = compilblrv(obj_sci, obj_template = obj_template, doplot = doplot, force = force, common_weights = False,
    #               get_cumul_plot = False, fcut = 0.8)

    keep = tbl1['per_epoch_err'] < 20
    tbl1 = tbl1[keep]
    #tbl2 = tbl2[keep]

    fig,ax = plt.subplots(nrows = 3, ncols = 1,sharex = True)

    rv_date1= []
    err_date1 = []
    mjdate_date1 = []

    #rv_date2= []
    #err_date2 = []
    #mjdate_date2 = []

    for date in np.unique(tbl1['DATE-OBS']):
        gg = tbl1['DATE-OBS'] == date
        rv = tbl1['per_epoch_mean'][gg]
        err = tbl1['per_epoch_err'][gg]

        mjdate_date1 = np.append(mjdate_date1,np.nanmean(tbl1['MJDATE'][gg]))
        rv_date1 = np.append(rv_date1, np.nansum(rv/err**2)/np.nansum(1/err**2))
        err_date1 = np.append(err_date1,np.sqrt(1/np.nansum(1/err**2)))

        #gg = tbl2['DATE-OBS'] == date
        #rv = tbl2['per_epoch_mean'][gg]
        #err = tbl2['per_epoch_err'][gg]

        #mjdate_date2 = np.append(mjdate_date2,np.nanmean(tbl2['MJDATE'][gg]))
        #rv_date2 = np.append(rv_date2, np.nansum(rv/err**2)/np.nansum(1/err**2))
        #err_date2 = np.append(err_date2,np.sqrt(1/np.nansum(1/err**2)))

    #rv_date1 -= np.nanmedian(rv_date1)
    #rv_date2 -= np.nanmedian(rv_date2)

    #ax.errorbar(mjdate_date2+.2,rv_date2, yerr=err_date2, color='red', fmt='.', alpha=0.5)
    #ax.plot(mjdate_date2,rv_date1-rv_date2,'o', color='black',  alpha=0.5)

    ax[0].errorbar(tbl1['MJDATE'],tbl1['per_epoch_mean'], yerr = tbl1['per_epoch_err'], fmt='.', alpha=0.5)
    ax[0].errorbar(mjdate_date1,rv_date1, yerr=err_date1, color='red', fmt='.', alpha=0.9)
    ax[1].errorbar(tbl1['MJDATE'],tbl1['per_epoch_DDV'], yerr = tbl1['per_epoch_DDVRMS'], fmt='.', alpha=0.5)
    ax[2].errorbar(tbl1['MJDATE'],tbl1['per_epoch_DDDV'], yerr = tbl1['per_epoch_DDDVRMS'], fmt='.', alpha=0.5)
    ax[0].set(xlabel = 'Date', ylabel = 'RV [m/s]')
    ax[1].set(xlabel = 'Date', ylabel = '2nd deriv [m/s]')
    ax[2].set(xlabel = 'Date', ylabel = '3rd deriv [m/s]')
    plt.show()