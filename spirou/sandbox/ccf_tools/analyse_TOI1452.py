import numpy as np
#from lin_mini import lin_mini # please don't ask Etienne, it's somewhere on this github repository!
import matplotlib.pyplot as plt
from astropy.table import Table
from bisector import *
from astropy.time import Time
from ccf2rv import *
from per_epoch_table import *

def sinusoidal(phase,dphase,amp,zp):
    return np.sin( (phase+dphase))*amp+zp

object = 'TOI-1452'
exclude_orders = [-1]
# number of median-absolute deviations within an epoch to consider a point discrepant
nMAD_cut = 5

# compare with and without sanitizing
tbl1 = get_object_rv(object,mask = 'gl846_neg',method = 'template',force = True,exclude_orders = exclude_orders,
                    snr_min = 20.0,sanitize = True, weight_type = 'ccf_depth', bandpass = 'YJHK')

tbl2 = get_object_rv(object,mask = 'gl846_neg',method = 'template',force = True,exclude_orders = exclude_orders,
                    snr_min = 20.0,sanitize = False, weight_type = 'ccf_depth', bandpass = 'YJHK')

fig, ax = plt.subplots(nrows = 2, ncols = 1)
ax[0].plot(tbl1['ERROR_RV'],'o', color = 'blue',label = 'error with sanitize')
ax[0].plot(tbl2['ERROR_RV'],'o', color = 'red',label = 'error without sanitize')
ax[0].set(xlabel = 'nth frame', ylabel ='error in m/s',title = object)
ax[1].plot(tbl2['ERROR_RV']/tbl1['ERROR_RV'],'ro')
ax[1].set(xlabel = 'nth frame', ylabel ='ratio not sanitize/sanitize')
ax[0].legend()
plt.show()

tbl_bin1 = per_epoch_table(tbl1)

plt.errorbar(tbl1['MJDATE'], tbl1['RV'],yerr=tbl1['ERROR_RV'], linestyle="None", fmt='o', alpha = 0.2)
plt.errorbar(tbl_bin1['MJDATE_MEAN'], tbl_bin1['RV'],yerr=tbl_bin1['ERROR_RV'], linestyle="None", fmt='o',
               alpha = 0.5, capsize = 2, color = 'black')
plt.xlabel('MJD')
plt.ylabel('Velocity [km/s]')
plt.title('Object = {0}, sanitized = {1}'.format(object, True))
plt.show()
