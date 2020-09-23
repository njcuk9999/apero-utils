import numpy as np
#from lin_mini import lin_mini # please don't ask Etienne, it's somewhere on this github repository!
import matplotlib.pyplot as plt
from astropy.table import Table
from bisector import *
from astropy.time import Time
from ccf2rv import *

def sinusoidal(phase,dphase,amp,zp):
    return np.sin( (phase+dphase))*amp+zp


if True:  # best setup for TOI-1278
    exclude_orders = [10,11,12,13,14,16,17,18,20,21,22,48]
    # number of median-absolute deviations within an epoch to consider a point discrepant
    nMAD_cut = 5
    tbl = get_object_rv('TOI-1278',mask = 'gl846_neg',method = 'template',force = True,exclude_orders = exclude_orders,
                        snr_min = 20.0)
    period = 14.4

if False:  # best setup for TOI-1278
    exclude_orders = [10,11,12,13,14,16,17,18,20,21,22,48]
    # number of median-absolute deviations within an epoch to consider a point discrepant
    nMAD_cut = 5
    tbl = get_object_rv('Gl699',mask = 'gl699_neg',method = 'template',force = True,exclude_orders = exclude_orders,
                        snr_min = 100.0)
    period = 14.4



#tbl = get_object_rv('TOI-1452',method = 'template',force = True,exclude_orders = exclude_orders)
#period = 11.064093

tbl_bin =

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

