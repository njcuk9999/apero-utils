import numpy as np
#from lin_mini import lin_mini # please don't ask Etienne, it's somewhere on this github repository!
import matplotlib.pyplot as plt
from astropy.table import Table
from bisector import *
from astropy.time import Time
from ccf2rv import *


exclude_orders = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,40,41,42,43,44,45,46,47,48]#10,11,12,13,14,16,17,18,20,21,22,48]
# number of median-absolute deviations within an epoch to consider a point discrepant
nMAD_cut = 5
tbl = get_object_rv('Gl699',mask = 'gl699_neg',method = 'bisector_40_60',force = True,exclude_orders = exclude_orders,
                    snr_min = 120.0,weight_type = '')


fig,ax = plt.subplots(nrows = 2, ncols=1)

ax[0].plot(tbl['RV'],'g.')
ax[1].plot(tbl['BIS_SLOPE'],'g.')
In [23]: plt.show()


tbl['RV'] -= np.nanmedian(tbl['RV'])

plt.plot(tbl['MJDATE'],tbl['RV'],'g.')
#plt.plot(tbl['MJDATE'],tbl['RV_WAVFP']-tbl['RV_SIMFP'],'r.')
plt.show()

#tbl = get_object_rv('TOI-1452',method = 'template',force = True,exclude_orders = exclude_orders)
#period = 11.064093

udates = np.unique(tbl['DATE-OBS'])

# Table with epoch-binned data
tbl_bin = Table()

tbl_bin['RV'] = np.zeros(len(udates), dtype = float)
tbl_bin['RV_SIG'] = np.zeros(len(udates), dtype = float)
tbl_bin['FORMAL_SIG'] = np.zeros(len(udates), dtype = float)

tbl_bin['RV_N'] = np.zeros(len(udates), dtype = int)
tbl_bin['MJDATE_MEAN'] = np.zeros(len(udates), dtype = float)
tbl_bin['MJDATE_MIN'] = np.zeros(len(udates), dtype = float)
tbl_bin['MJDATE_MAX'] = np.zeros(len(udates), dtype = float)

tbl['RV'] -= np.nanmean(tbl['RV'])
tbl['RV'] = 1000*tbl['RV']

tbl['KEEP'] = np.ones_like(tbl,dtype = bool)
# Table per night
for ite in range(1):
    tbl2 = Table()
    for i in range(len(udates)):
        g = (udates[i] == tbl['DATE-OBS'])

        nMAD = (tbl['RV'][g] - np.nanmedian(tbl['RV'][g]))/np.nanmedian(np.abs(tbl['RV'][g] - np.nanmedian(tbl['RV'][g])))
        nMAD = np.abs(nMAD)
        print(np.max(nMAD))

        tbl['KEEP'][g] = nMAD<nMAD_cut

        tbl_bin['RV'][i] = np.mean(tbl['RV'][g])
        tbl_bin['RV_SIG'][i] = np.std(tbl['RV'][g])
        tbl_bin['RV_N'][i] = np.sum(g)


        tbl_bin['MJDATE_MEAN'][i] = np.mean(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MIN'][i] = np.min(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MAX'][i] = np.max(tbl['MJDATE'][g])

        tbl_bin['FORMAL_SIG'][i] =  tbl_bin['RV_SIG'][i]/np.sqrt(tbl_bin['RV_N'][i])

    print(len(tbl))
    tbl = tbl[tbl['KEEP']]


stop
t2 = Time(tbl_bin['MJDATE_MEAN'], format = 'mjd')
t3 = Time(tbl['MJDATE'], format = 'mjd')

phase_bin = 2*np.pi*tbl_bin['MJDATE_MEAN']/period
phase = 2*np.pi*tbl['MJDATE']/period


fit, pcov = curve_fit(sinusoidal, phase_bin, tbl_bin['RV'])


dt = np.max(tbl_bin['MJDATE_MEAN']) - np.min(tbl_bin['MJDATE_MEAN'])
time_plot = np.arange(np.min(tbl_bin['MJDATE_MEAN'])-dt/10,np.max(tbl_bin['MJDATE_MEAN'])+dt/10,dt/1000)
phase_plot = 2*np.pi*time_plot/period

#model_plot = #amps[0]+amps[1]*np.sin(2*np.pi*time_plot/period)+amps[2]*np.cos(2*np.pi*time_plot/period)
model_bin =  sinusoidal(phase_bin,*fit)  #amps[0]+amps[1]*np.sin(phase)+amps[2]*np.cos(phase)
model=  sinusoidal(phase,*fit)  #amps[0]+amps[1]*np.sin(phase)+amps[2]*np.cos(phase)
model_plot =  sinusoidal(phase_plot,*fit)  #amps[0]+amps[1]*np.sin(phase)+amps[2]*np.cos(phase)

print('Amplitude of the sinusoidal at {0} days: {1} m/s'.format(period,fit[1]))
print('Mean/Median per-epoch STDDEV {0}/{1} m/s'.format(np.mean(tbl_bin["RV_SIG"]),np.median(tbl_bin["RV_SIG"])))


fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True)

for i in range(len(t2)):
    ax[0].plot_date(t2.plot_date,tbl_bin['RV'],'g.')
    ax[0].plot_date([t2[i].plot_date,t2[i].plot_date],[tbl_bin['RV'][i]-tbl_bin['FORMAL_SIG'][i],tbl_bin['RV'][i]+tbl_bin['FORMAL_SIG'][i]],'g')

ax[0].plot_date(t3.plot_date,tbl['RV'],'r.',alpha = 0.5)
ax[1].plot_date(t3.plot_date,tbl['RV'] - model,'r.',alpha = 0.5)


ax[0].plot(Time(time_plot, format = 'mjd').plot_date,model_plot,'r:')
ax[0].set(ylabel = 'Residuals [m/s]')

for i in range(len(t2)):
    ax[1].plot_date(t2.plot_date,tbl_bin['RV']-model_bin,'g.')
    ax[1].plot_date([t2[i].plot_date,t2[i].plot_date],[tbl_bin['RV'][i]-tbl_bin['FORMAL_SIG'][i]-model_bin[i],tbl_bin['RV'][i]+tbl_bin['FORMAL_SIG'][i]-model_bin[i]],'g')

ax[1].plot(Time(time_plot, format = 'mjd').plot_date,np.zeros(len(time_plot)),'r:')
ax[1].set(xlabel = 'Date', ylabel = 'Residuals [m/s]')
plt.tight_layout()
plt.show()

