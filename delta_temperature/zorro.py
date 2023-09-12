import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

# tbl = Table.read('lbl_PROXIMA_PROXIMA.rdb', format='ascii.rdb')
# title = 'Proxima'
# tbl = Table.read('lbl_TOI756_GL514.rdb', format='ascii.rdb')
# title = 'TOI756'
tbl1 = Table.read('/Volumes/courlan/lbl_SPIROU/lblrdb/lbl_TRAPPIST1PC_tc_TRAPPIST1PC_tc.rdb', format='ascii.rdb')
tbl2 = Table.read('/Users/eartigau/Downloads/lbl_TRAPPIST1_TRAPPIST1.rdb', format='ascii.rdb')
title = 'T1'

tbl1['phase'] = ((tbl1['rjd'] - 57322.51736) / 1.51087081 + 0.5) % 1
tbl2['phase'] = ((tbl2['rjd'] - 57322.51736) / 1.51087081 + 0.5) % 1

for b in np.unique(tbl1['DATE-OBS']):
    g = np.where(tbl1['DATE-OBS'] == b)[0]
    if len(g) < 12:
        tbl1['DTEMP'][g] = np.nan
    else:
        # tbl1['DTEMP'][g] -= np.mean(tbl1['DTEMP'][g])
        tbl1['DTEMP'][g] -= np.polyval(np.polyfit(tbl1['phase'][g], tbl1['DTEMP'][g], 3), tbl1['phase'][g])

for b in np.unique(tbl2['DATE-OBS']):
    g = np.where(tbl2['DATE-OBS'] == b)[0]
    if len(g) < 12:
        tbl2['DTEMP'][g] = np.nan
    else:
        # tbl2['DTEMP'][g] -= np.mean(tbl2['DTEMP'][g])
        tbl2['DTEMP'][g] -= np.polyval(np.polyfit(tbl2['phase'][g], tbl2['DTEMP'][g], 3), tbl2['phase'][g])

tbl1 = tbl1[np.isfinite(tbl1['DTEMP'])]
tbl2 = tbl2[np.isfinite(tbl2['DTEMP'])]

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 8), sharex=True, sharey=True)

for b in np.unique(tbl1['DATE-OBS']):
    g = np.where(tbl1['DATE-OBS'] == b)[0]
    ax[0].errorbar(tbl1['phase'][g], tbl1['DTEMP'][g], yerr=tbl1['sDTEMP'][g], fmt='o-', alpha=0.5, label=b)
for b in np.unique(tbl2['DATE-OBS']):
    g = np.where(tbl2['DATE-OBS'] == b)[0]
    ax[1].errorbar(tbl2['phase'][g], tbl2['DTEMP'][g], yerr=tbl2['sDTEMP'][g], fmt='o-', alpha=0.5)

ax[0].set(xlabel='Phase', ylabel='Delta Temperature [K]', title='Avec correction de persistance')
ax[1].set(xlabel='Phase', title='Sans correction de persistance')
# plt.title(title)
ax[0].legend()
ax[0].grid()
ax[1].grid()
plt.tight_layout()
plt.savefig('{}_dtemp.png'.format(title))
plt.show()
