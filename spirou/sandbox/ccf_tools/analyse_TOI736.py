import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from bisector import *
from astropy.time import Time
from ccf2rv import *
from per_epoch_table import per_epoch_table

def sinusoidal(phase,dphase,amp,zp):
    return np.sin( (phase+dphase))*amp+zp

# do not *formally* exclude an order, but this is done later with the bandpass keyword
exclude_orders = [28,47,48]

object = 'TOI-736'
mask =  'gl699_neg'
method = 'all'
sanitize = True
# number of median-absolute deviations within an epoch to consider a point discrepant
tbl,dico = get_object_rv(object,mask =mask,
                    method = method,force = True,
                    exclude_orders = exclude_orders,
                    snr_min = 20.0, velocity_window = 20, sanitize = sanitize,
                    dvmax_per_order = 500.0, bandpass = 'H',
                    doplot = True, do_blacklist = True,
                    detailed_output = True,
                    sed_match = False)


rv = np.array(tbl['RV'])
rv -= np.mean(rv)

ccf = np.array(dico['MEAN_CCF'])

ccf2 = np.array(ccf)
for i in range(34):
    ccf2[:,i] = np.roll(ccf2[:,i],int(-rv[i]*10))

moy = np.mean(ccf2,axis=1)
for i in range(34):
    ccf2[:,i] -= moy

for i in range(34):
    ccf2[:,i] = np.roll(ccf2[:,i],int(rv[i]*10))

damps = np.arange(10,55,0.1)

all_ccs = np.zeros([ccf2.shape[0],len(damps)])

for ite in range(len(damps)):
    print(ite)
    ccf3 = np.zeros_like(ccf2)
    for i in range(34):
        ccf3[:,i] = np.roll(ccf2[:,i],int(damps[ite]*rv[i]*10))

    all_ccs[:,ite] = np.nanmean(ccf3,axis=1)

plt.plot(dico['ccf_RV'],moy)
plt.show()

plt.plot(damps,all_ccs[np.argmin(moy),:])
plt.show()


plt.imshow(all_ccs/np.std(all_ccs),aspect = 'auto',extent = [np.min(damps),np.max(damps),np.min(dico['ccf_RV']),np.max(dico['ccf_RV'])])
plt.show()

# period for the sinusoidal currve
period = 14.4

# create the table with bis per epoch
tbl_bin = per_epoch_table(tbl,nMAD_cut = 5)

# get time stamps friendly for plotting
t2 = Time(tbl_bin['MJDATE_MEAN'], format = 'mjd')
t3 = Time(tbl['MJDATE'], format = 'mjd')

# get phase for sine fitting
phase_bin = 2*np.pi*tbl_bin['MJDATE_MEAN']/period
phase = 2*np.pi*tbl['MJDATE']/period

# fit sinusoid
fit, pcov = curve_fit(sinusoidal, phase_bin, tbl_bin['RV'])

# some plotting fiddling
dt = np.max(tbl_bin['MJDATE_MEAN']) - np.min(tbl_bin['MJDATE_MEAN'])
time_plot = np.arange(np.min(tbl_bin['MJDATE_MEAN'])-dt/10,np.max(tbl_bin['MJDATE_MEAN'])+dt/10,dt/1000)
phase_plot = 2*np.pi*time_plot/period

model_bin =  sinusoidal(phase_bin,*fit)
model=  sinusoidal(phase,*fit)
model_plot =  sinusoidal(phase_plot,*fit)

print('Amplitude of the sinusoidal at {0} days: {1:.2f} m/s'.format(period, 1000*fit[1]))
print('Mean velocity: {1:.2f} m/s'.format(period, 1000*fit[2]))

print('Mean/Median per-epoch STDDEV {0}/{1} km/s'.format(np.mean(tbl_bin["ERROR_RV"])
                                                        ,np.median(tbl_bin["ERROR_RV"])))

fig, ax = plt.subplots(nrows = 2, ncols = 1,sharex = True, figsize = (14,8))

for i in range(len(t2)):
    ax[0].plot_date(t2.plot_date,tbl_bin['RV'],'g.')
    ax[0].plot_date([t2[i].plot_date,t2[i].plot_date],[tbl_bin['RV'][i]-tbl_bin['ERROR_RV'][i],
                                                       tbl_bin['RV'][i]+tbl_bin['ERROR_RV'][i]],'g')
ax[0].plot_date(t3.plot_date,tbl['RV'],'r.',alpha = 0.5)
ax[1].errorbar(t3.plot_date,tbl['RV'] - model,yerr=tbl['ERROR_RV'], linestyle="None",
               fmt='o',color = 'green', alpha = 0.2, label = 'Individual measurements')

ax[0].plot(Time(time_plot, format = 'mjd').plot_date,model_plot,'r:')
ax[0].set(ylabel = 'Velocity [km/s]',title = object)


ax[1].errorbar(t2.plot_date, tbl_bin['RV'] - model_bin, yerr=tbl_bin['ERROR_RV'],
               linestyle="None", fmt='o',
               alpha = 0.5, capsize = 2, color = 'black',label = 'Epoch mean')

ax[1].legend()
ax[1].plot(Time(time_plot, format = 'mjd').plot_date,np.zeros(len(time_plot)),'r:')
ax[1].set(xlabel = 'Date', ylabel = 'Residuals [km/s]',ylim = [-.15,0.15],
          xlim = [np.min(Time(time_plot, format = 'mjd').plot_date),
          np.max(Time(time_plot, format = 'mjd').plot_date)]
          )

for label in ax[1].get_xticklabels():
  label.set_rotation(25)
  label.set_ha('right')

plt.tight_layout()
plt.savefig(object+'.pdf')
plt.show()

sigma = np.std((tbl_bin['RV'] - model_bin))
mean_error = np.mean(tbl_bin['ERROR_RV'])
median_error = np.nanmedian(tbl_bin['ERROR_RV'])
reduced_chi2 = np.std((tbl_bin['RV'] - model_bin)/tbl_bin['ERROR_RV'])
print('\n--- values for the per-night weighted-mean points ---\n')
print(' mean ERROR_RV {0:.2f} m/s, median ERROR_RV {1:.2f} m/s, '
      'reduced chi2 {2:.2f} '.format(mean_error*1e3, median_error*1e3, reduced_chi2))

mean_error = np.mean(tbl['ERROR_RV'])
median_error = np.nanmedian(tbl['ERROR_RV'])
print('\n--- values for the individual points ---\n')
print(' mean ERROR_RV {0:.2f} m/s, median ERROR_RV {1:.2f} m/s'.format( mean_error*1e3,median_error*1e3))


f = open('TOI1278_obslog.tex','w')
# create an observation log in tex format
# Nice when you want to write a paper in the end, hey, that's the point of all these observations!
for i in range(len(tbl)):
    f.write('{0:.4f} & ${1:.3f} \pm {2:.3f}$ & {3:.3f} \\\\ \n'.format(tbl['MJDATE'][i],tbl['RV'][i], tbl['ERROR_RV'][i],tbl['D2_RESIDUAL_CCF'][i]))
f.close()