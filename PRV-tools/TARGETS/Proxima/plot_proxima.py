import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

# TODO change the path to match your local path
file_path = '/Users/eartigau/Downloads/lbl_PROXIMA_PROXIMA.rdb'

# read into a table
tbl = Table.read(file_path, format='ascii.rdb')
tbl = tbl[tbl['vrad'] < (-21215)] # removing outliers

plt.plot_date(tbl['plot_date'], tbl['vrad'], 'k.',alpha = 0.1)
plt.errorbar(tbl['plot_date'], tbl['vrad'], yerr=tbl['svrad'], fmt='none', ecolor='k', alpha=0.5)
plt.show()