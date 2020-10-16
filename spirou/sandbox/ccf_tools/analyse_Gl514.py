import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.time import Time
from ccf2rv import get_object_rv


exclude_orders = [0, 11, 12, 13, 15, 16, 20, 21, 22, 47, 48]
# number of median-absolute devs within an epoch to consider a point discrepant
nMAD_cut = 5


tbl = get_object_rv('Gl514', mask='gl514_neg', method='all',
                    force=True, exclude_orders=exclude_orders,
                    sanitize=True,
                    snr_min=120, weight_type='', bandpass='HK',
                    velocity_window=10, do_blacklist=True)

fig, ax = plt.subplots(nrows=2, ncols=1)

ax[0].plot(tbl['RV'], 'g.')
ax[1].plot(tbl['BIS_SLOPE'], 'g.')
tbl['RV'] -= np.nanmedian(tbl['RV'])

plt.plot(tbl['MJDATE'], tbl['RV'], 'g.')
# plt.plot(tbl['MJDATE'],tbl['RV_WAVFP']-tbl['RV_SIMFP'],'r.')
plt.show()

udates = np.unique(tbl['DATE-OBS'])

# Table with epoch-binned data
tbl_bin = Table()

tbl_bin['RV'] = np.zeros(len(udates), dtype=float)
tbl_bin['RV_SIG'] = np.zeros(len(udates), dtype=float)
tbl_bin['FORMAL_SIG'] = np.zeros(len(udates), dtype=float)

tbl_bin['RV_N'] = np.zeros(len(udates), dtype=int)
tbl_bin['MJDATE_MEAN'] = np.zeros(len(udates), dtype=float)
tbl_bin['MJDATE_MIN'] = np.zeros(len(udates), dtype=float)
tbl_bin['MJDATE_MAX'] = np.zeros(len(udates), dtype=float)

tbl['RV'] -= np.nanmean(tbl['RV'])
tbl['RV'] = 1000*tbl['RV']

tbl['KEEP'] = np.ones_like(tbl, dtype=bool)
# Table per night
for ite in range(1):
    tbl2 = Table()
    for i in range(len(udates)):
        g = (udates[i] == tbl['DATE-OBS'])

        nMAD = ((tbl['RV'][g] - np.nanmedian(tbl['RV'][g]))
                / np.nanmedian(
                    np.abs(tbl['RV'][g] - np.nanmedian(tbl['RV'][g]))
                    )
                )
        nMAD = np.abs(nMAD)
        print(np.max(nMAD))

        tbl['KEEP'][g] = nMAD < nMAD_cut

        tbl_bin['RV'][i] = np.mean(tbl['RV'][g])
        tbl_bin['RV_SIG'][i] = np.std(tbl['RV'][g])
        tbl_bin['RV_N'][i] = np.sum(g)

        tbl_bin['MJDATE_MEAN'][i] = np.mean(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MIN'][i] = np.min(tbl['MJDATE'][g])
        tbl_bin['MJDATE_MAX'][i] = np.max(tbl['MJDATE'][g])

        tbl_bin['FORMAL_SIG'][i] = (tbl_bin['RV_SIG'][i]
                                    / np.sqrt(tbl_bin['RV_N'][i])
                                    )

    print(len(tbl))
    tbl = tbl[tbl['KEEP']]

t2 = Time(tbl_bin['MJDATE_MEAN'], format='mjd')
t3 = Time(tbl['MJDATE'], format='mjd')

dt = np.max(tbl_bin['MJDATE_MEAN']) - np.min(tbl_bin['MJDATE_MEAN'])
time_plot = np.arange(
        np.min(tbl_bin['MJDATE_MEAN'])-dt/10,
        np.max(tbl_bin['MJDATE_MEAN'])+dt/10,
        dt/1000)

print('Mean/Median per-epoch STDDEV {0}/{1} m/s'
      .format(np.mean(tbl_bin["RV_SIG"]), np.median(tbl_bin["RV_SIG"])))


fig, ax = plt.subplots(nrows=1, ncols=1, sharex=True)

for i in range(len(t2)):
    ax.plot_date(t2.plot_date, tbl_bin['RV'], 'g.')
    ax.plot_date([t2[i].plot_date, t2[i].plot_date],
                 [tbl_bin['RV'][i]-tbl_bin['FORMAL_SIG'][i],
                  tbl_bin['RV'][i]+tbl_bin['FORMAL_SIG'][i]],
                 'g')

ax.plot_date(t3.plot_date, tbl['RV'], 'r.', alpha=0.5)

ax.set(ylabel='RV [m/s]')

plt.tight_layout()
plt.show()
