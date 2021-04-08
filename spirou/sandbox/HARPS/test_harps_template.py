from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table

tbl1 = Table.read('templates/Template_GL699_HARPS.fits')
tbl2 = Table.read('templates/Template_GL699_tcorr_HARPS.fits')


plt.plot(tbl1['wavelength'],tbl1['flux'],alpha =0.5, label = 'original')
plt.plot(tbl2['wavelength'],tbl2['flux'],alpha =0.5, label = 'corr')
plt.plot(tbl2['wavelength'],tbl2['flux']/tbl1['flux'],alpha =0.5, label = 'original/corr')
plt.ylim([0,2])
plt.legend()
plt.show()