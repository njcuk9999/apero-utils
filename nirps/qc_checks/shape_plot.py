import os
import numpy as np
from astropy.io import fits
from astropy.time import Time
import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

path = '/cosmos99/nirps/apero-data/nirps_he_online/calib/'

files = list(glob.glob(os.path.join(path, '*shapel.fits')))


mjd = np.zeros(len(files))
dx = np.zeros(len(files))
dy = np.zeros(len(files))
a = np.zeros(len(files))
b = np.zeros(len(files))
c = np.zeros(len(files))
d = np.zeros(len(files))



for it in tqdm(range(len(files))):
    filename = files[it]
    # open header
    hdr = fits.getheader(filename)
    # push the values into numpy arrays
    mjd[it] = float(hdr['MJDMID'])
    dx[it] = float(hdr['SHAPE_DX'])
    dy[it] = float(hdr['SHAPE_DY'])
    a[it] = float(hdr['SHAPE_A'])
    b[it] = float(hdr['SHAPE_B'])
    c[it] = float(hdr['SHAPE_C'])
    d[it] = float(hdr['SHAPE_D'])


# make mjd a time array
mjd = Time(mjd, format='mjd').to_datetime()

# Plot
plt.close()
fig, frames = plt.subplots(ncols=2, nrows=3, sharex='all')
# Invisible proxy for the legend (no handle, just text)
invisible_proxy = Line2D([0], [0], color='none')
# flatten frames
fframes = frames.flatten()
pkwargs = dict(marker='.', ls='None')
variables = [dx, dy, 1-a, c, b, 1-d]
labels = ['dx', 'dy', '1-a', 'b', 'c', '1-d']

for it, frame in enumerate(fframes):
    # plot the data
    frame.plot(mjd, variables[it], **pkwargs, label=labels[it])
    # set the labels
    frame.set(xlabel='Date', ylabel=labels[it])
    # set the xlabels to vertical
    frame.tick_params(axis='x', rotation=90)
    # Set the y-axis formatter to show scientific notation
    frame.yaxis.set_major_formatter(ScalarFormatter())
    frame.yaxis.get_major_formatter().set_powerlimits((-10, 10))
    # Set date formatting to YYYY-mm-dd HH:MM:SS
    frame.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
    # Move y-axis to the right for odd columns (1, 3, 5 in 0-based indexing)
    if it % 2 == 1:
        frame.yaxis.tick_right()
        frame.yaxis.set_label_position("right")  # Move the y-axis label as well
    # plot legend
    frame.legend([invisible_proxy], [labels[it]], handletextpad=0,
                 handlelength=0, loc=0, frameon=True)


plt.subplots_adjust(bottom=0.2, hspace=0.01, wspace=0.01)
plt.show()


