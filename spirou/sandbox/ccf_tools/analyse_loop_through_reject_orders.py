import numpy as np
import matplotlib.pyplot as plt
from ccf2rv import get_object_rv


def sinusoidal(phase, dphase, amp, zp):
    return np.sin(phase+dphase)*amp+zp


obj = 'TOI-1452'
# number of median-absolute devs within an epoch to consider a point discrepant
nMAD_cut = 5

# you can give some known bad orders
known_bad = [36]

# get the RMS of errors for the case where not additional order is rejected
tbl1 = get_object_rv(
        obj,
        mask='gl846_neg',
        method='template',
        force=True,
        exclude_orders=known_bad,
        snr_min=20.0,
        sanitize=False,
        weight_type='ccf_depth',
        bandpass='YJHK',
        doplot=False,
        )
ref_rms = np.mean(tbl1['ERROR_RV'])

ERROR_RV = np.zeros(49)

for iord in range(0, 49):
    exclude_orders = np.array(np.unique(np.append(iord, known_bad)), dtype=int)

    print('\nWe try rejecting order[s] : {0}\n'.format(exclude_orders))
    tbl1 = get_object_rv(
            obj,
            mask='gl846_neg',
            method='template',
            force=True,
            exclude_orders=exclude_orders,
            snr_min=20.0,
            sanitize=False,
            weight_type='ccf_depth',
            bandpass='YJHK',
            doplot=False,
            )

    ERROR_RV[iord] = np.mean(tbl1['ERROR_RV'])/ref_rms
    print(iord, ERROR_RV[iord])

plt.plot(ERROR_RV, 'ko')
plt.xlabel('Nth order rejected')
plt.ylabel('RMS relative to no rejection')
plt.title('Testing rejection of orders')
plt.show()
