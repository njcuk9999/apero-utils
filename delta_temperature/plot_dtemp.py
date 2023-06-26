import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table

tbl = Table.read('lbl_TOI756_GL514.rdb.txt', format='ascii.rdb')
title = 'TOI-756 [Template: GL514]'

fig, ax = plt.subplots(nrows = 3, ncols=1,figsize=(15,8),sharex = True  )

period = np.inf

ax[0].plot_date(tbl['plot_date'], tbl['vrad'],'k.',alpha = 0.5)
ax[1].plot_date(tbl['plot_date'], tbl['DTEMP'],'k.',alpha = 0.5)
ax[2].plot_date(tbl['plot_date'], tbl['fwhm'],'k.',alpha = 0.5)

ax[0].errorbar(tbl['plot_date'], tbl['vrad'],yerr = tbl['svrad'],fmt= 'k.',alpha = 0.5)
ax[1].errorbar(tbl['plot_date'], tbl['DTEMP'], yerr=tbl['sDTEMP'], fmt='k.', ecolor='k', alpha=0.5)
ax[2].errorbar(tbl['plot_date'], tbl['fwhm'],yerr = tbl['sig_fwhm'],fmt= 'k.',alpha = 0.5)

ax[0].set(ylabel = 'RV [m/s]',title = title)
ax[1].set(ylabel = 'Delta Temperature [K]')
ax[2].set(ylabel = 'FWHM [m/s]',xlabel = 'Date')

ax[0].grid()
ax[1].grid()
ax[2].grid()
plt.tight_layout()
plt.savefig('{}_rv_dtemp_fwhm.png'.format(title))
plt.show()