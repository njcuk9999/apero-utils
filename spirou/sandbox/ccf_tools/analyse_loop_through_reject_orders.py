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
# number of median-absolute deviations within an epoch to consider a point discrepant
nMAD_cut = 5


rms = np.zeros(49)
for iord in range(49):
    exclude_orders = [iord]
    # compare with and without sanitizing
    tbl1 = get_object_rv(object,mask = 'gl846_neg',method = 'template',force = True,exclude_orders = exclude_orders,
                        snr_min = 20.0,sanitize = True, weight_type = 'ccf_depth', bandpass = 'YJHK', doplot = False)

    rms[iord] = np.mean(tbl1['ERROR_RV'])
    print(iord,rms[iord])