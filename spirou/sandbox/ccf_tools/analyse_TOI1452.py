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


exclude_orders = [0,11,12,13,15,16,20,21,22,47,48]
# number of median-absolute deviations within an epoch to consider a point discrepant
nMAD_cut = 5
tbl1 = get_object_rv('TOI-1452',mask = 'gl846_neg',method = 'template',force = True,exclude_orders = exclude_orders,
                    snr_min = 20.0,sanitize = True, weight_type = 'ccf_depth')

#tbl2 = get_object_rv('TOI-1452',mask = 'gl846_neg',method = 'template',force = True,exclude_orders = exclude_orders,
#                    snr_min = 20.0,sanitize = False, weight_type = 'ccf_depth')

tbl_bin1 = per_epoch_table(tbl1)

tbl['RV'] -= np.nanmedian(tbl['RV'])
tbl = tbl[tbl['CCF_RESIDUAL_RMS'] <0.008]

plt.errorbar(tbl1['MJDATE'], tbl1['RV'],yerr=tbl1['ERROR_RV'], linestyle="None", fmt='o')
plt.show()

print(np.nanstd(tbl['RV']  )  ,np.nanmedian(np.abs(tbl['RV'])))
print(np.nanmedian(np.abs(tbl['RV'][1:]-tbl['RV'][:-1]))/np.sqrt(2) )