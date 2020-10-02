import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from bisector import *
from astropy.time import Time
from ccf2rv import *
from per_epoch_table import *

object = 'K2-33'
exclude_orders = [47]
# number of median-absolute deviations within an epoch to consider a point discrepant
nMAD_cut = 5
method = 'template'
mask = 'gl846_neg'

tbl1 = get_object_rv(object,
                    bandpass='HK',
                    method=method,
                    dvmax_per_order = 10.0, # large because it's a somewhat fast rotator
                    exclude_orders = exclude_orders,
                    force = True,
                    snr_min = 25.0,
                    weight_type = '',
                    sanitize = True,
                    doplot = True,
                    mask = mask,
                    do_blacklist = True)


tbl_bin1 = per_epoch_table(tbl1)

fig, ax = plt.subplots(nrows = 2, ncols= 1,sharex = True)

ax[0].errorbar(tbl1['MJDATE'], tbl1['RV'],yerr=tbl1['ERROR_RV'], linestyle="None", fmt='o', alpha = 0.2)
ax[0].errorbar(tbl_bin1['MJDATE_MEAN'], tbl_bin1['RV'],yerr=tbl_bin1['ERROR_RV'], linestyle="None", fmt='o',
               alpha = 0.5, capsize = 2, color = 'black')

ax[1].set(xlabel = 'MJD')
ax[0].set(ylabel = 'Velocity [km/s]')
ax[1].plot(tbl1['MJDATE'],tbl1['D3_RESIDUAL_CCF'],'g.')
ax[1].set(ylabel = 'd2 activity parameter')
ax[0].set(title = 'Object = {0}, sanitized = {1}'.format(object, True))
plt.savefig('K2-33_velocity.pdf')
plt.show()
