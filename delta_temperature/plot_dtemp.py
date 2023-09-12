import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table

#tbl = Table.read('lbl_PROXIMA_PROXIMA.rdb', format='ascii.rdb')
#title = 'Proxima'
#tbl = Table.read('lbl_TOI756_GL514.rdb', format='ascii.rdb')
#title = 'TOI756'
tbl = Table.read('/Users/eartigau/Desktop/lbl_Gl873_tc_Gl873_tc.rdb', format='ascii.rdb')
title = 'T1'


for k in tbl.keys():
    if 'nm' not in k:
        continue
    if k.startswith('svrad'):
        continue

    plt.errorbar(tbl['plot_date'], tbl[k] - np.nanmedian(tbl[k]), yerr=tbl['s'+k], fmt='.', alpha=0.5, label=k)
plt.legend()
plt.show()


fig, ax = plt.subplots(nrows = 5, ncols=1,figsize=(15,8),sharex = True  )

#tbl['plot_date'] = tbl['plot_date'] % 5.1624

period = np.inf

ax[0].plot_date(tbl['plot_date'], tbl['vrad_605nm'],'k.',alpha = 0.5)
ax[1].plot_date(tbl['plot_date'], tbl['DTEMP'],'k.',alpha = 0.5)
ax[2].plot_date(tbl['plot_date'], tbl['fwhm'],'k.',alpha = 0.5)
ax[3].plot_date(tbl['plot_date'], tbl['contrast'],'k.',alpha = 0.5)
ax[4].plot_date(tbl['plot_date'], tbl['vrad_chromatic_slope'],'k.',alpha = 0.5)

ax[0].errorbar(tbl['plot_date'], tbl['vrad_605nm'],yerr = tbl['svrad_605nm'],fmt= 'k.',alpha = 0.5)
ax[1].errorbar(tbl['plot_date'], tbl['DTEMP'], yerr=tbl['sDTEMP'], fmt='k.', ecolor='k', alpha=0.5)
ax[2].errorbar(tbl['plot_date'], tbl['fwhm'],yerr = tbl['sig_fwhm'],fmt= 'k.',alpha = 0.5)
ax[3].errorbar(tbl['plot_date'], tbl['contrast'],yerr = tbl['sig_contrast'],fmt= 'k.',alpha = 0.5)
ax[4].errorbar(tbl['plot_date'], tbl['vrad_chromatic_slope'],yerr = tbl['svrad_chromatic_slope'],fmt= 'k.',alpha = 0.5)

ax[0].set(ylabel = 'RV [m/s]',title = title)
ax[1].set(ylabel = 'Delta Temperature [K]')
ax[2].set(ylabel = 'FWHM [m/s]',xlabel = 'Date')
ax[3].set(ylabel = 'Contrast',xlabel = 'Date')
ax[4].set(ylabel = 'Chromatic slope m/s/µm',xlabel = 'Date')

ax[0].grid()
ax[1].grid()
ax[2].grid()
ax[3].grid()
ax[4].grid()
plt.tight_layout()
plt.savefig('{}_rv_dtemp_fwhm.png'.format(title))
plt.show()