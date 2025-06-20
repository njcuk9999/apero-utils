from astropy.table import Table
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
from etienne_tools import snail, sigma
import numpy as np
import matplotlib.cm as cm
from matplotlib.colors import Normalize

all_fps = np.array(glob.glob('/Users/eartigau/glitch_fp/data_NIRPS_HE/*wave_fplines_A.fits'))

mjd = np.array([hdr['MJD-OBS'] for hdr in [fits.getheader(f) for f in snail(all_fps)]])

ord = np.argsort(mjd)
all_fps = all_fps[ord]
mjd = mjd[ord]

# we set colors scaling with time using a heat palette. The heat
# is proportional to the mjd - the first mjd
colors = cm.rainbow((mjd - mjd.min()) / (mjd.max() - mjd.min()))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))
for i in snail(range(len(all_fps))):
    f = all_fps[i]
    tbl = Table.read(f)
    w = tbl['WAVE_MEAS'].data.data
    keep = (w > 1600) * (w < 1605)
    tbl = tbl[keep]

    wave = tbl['WAVE_MEAS'].data.data
    nsig = tbl['NSIG'].data.data

    off = (tbl['WAVE_MEAS'] / tbl['WAVE_REF'] - 1).data.data * 3e8
    sig1 = sigma(off * nsig)
    sig = sig1 / nsig

    ord = np.argsort(wave)
    wave = wave[ord]
    off = off[ord]
    sig = sig[ord]

    x = wave
    y = off
    yerr = sig
    nsig = y / yerr

    valid = np.isfinite(x + y + yerr) & (np.abs(nsig) < 8)
    x = x[valid]
    y = y[valid]
    yerr = yerr[valid]

    ax.errorbar(x, y, yerr=yerr, fmt='.-', alpha=.2, label=f'{mjd[i]:.2f}', color=colors[i])

# add a colorbar with the correct range
norm = Normalize(vmin=mjd.min(), vmax=mjd.max())
fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.rainbow), ax=ax, label='MJD')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('FP position error (m/s)')
plt.show()