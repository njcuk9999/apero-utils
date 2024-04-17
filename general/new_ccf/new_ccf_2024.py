# ius spline
import glob

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.constants import c
# ius spline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from tqdm import tqdm


def gauss_slope(x, cen, w, depth, slope, amp):
    """

    Parameters
    ----------
    x --> velocity grid
    cen --> center of gaussian
    w --> e-width of the gaussian
    depth --> normalized depth to continuum
    slope --> slope of the continuum
    amp --> amplitude of the continuum. Normally very close to 1.

    Returns
    -------
    a gaussian with the properties listed on the x grid
    """
    xpix = (x - cen)
    gauss = (1 - depth * np.exp(-2 * xpix ** 2 / w ** 2))
    return gauss * (1 - slope * xpix) * amp


def doppler(wave0, v):
    # velocity expressed in m/s
    # relativistic calculation

    v = np.array(v)
    return wave0 * np.sqrt((1 - v / c) / (1 + v / c))


linelist = '/Users/eartigau/ccf/PROXIMA_neg_ccf.fits'
tbl = dict(Table.read(linelist))

# find all files
search_path = '/Users/eartigau/ccf/data/*fits'
files = np.array(glob.glob(search_path))
files.sort()

v0 = -5e4  # m/s
v1 = 0.0  # m/s
velostep = 500  # m/s
dvs = np.arange(v0, v1, velostep)

# placeholders for the results
RVS = np.zeros(len(files))
sig_RVS = np.zeros(len(files))
rjd = np.zeros(len(files))
bisspan = np.zeros(len(files))
contrasts = np.zeros(len(files))
fwhm = np.zeros(len(files))

for ifile, f in enumerate(files):
    print('Processing file {}/{}'.format(ifile, len(files)))
    sp = fits.getdata(f)
    h = fits.getheader(f)
    rjd[ifile] = h['MJD-OBS']
    berv = h['BERV'] * 1e3
    wave = fits.getdata('/Users/eartigau/ccf/wavesol/' + h['WAVEFILE'])
    wave = doppler(wave, -berv)

    ccfs = np.zeros((len(dvs), sp.shape[0]))
    for iord in tqdm(range(sp.shape[0]), leave=False):
        ccf = np.zeros_like(dvs)
        w1 = np.min(wave[iord, :])
        w2 = np.max(wave[iord, :])
        g = (tbl['ll_mask_s'] > w1) & (tbl['ll_mask_s'] < w2)

        if np.sum(g) < 20:
            # if <20 lines with mask in this order, skip
            continue

        ww = tbl['ll_mask_s'][g]
        wmask = tbl['w_mask'][g]

        # wave of the order considered
        sp2 = sp[iord]
        w2 = wave[iord]
        # valid domain
        valid = np.isfinite(sp2)
        # spline to interpolate the spectrum
        spl = InterpolatedUnivariateSpline(w2[valid], sp2[valid], k=2, ext=1)
        # spline to monitor if line falls on a NaN
        spl_valid = InterpolatedUnivariateSpline(w2, np.array(valid, dtype=float), k=1, ext=1)
        valid2 = []
        for i, dv in enumerate(dvs):
            valid2.append(spl_valid(doppler(ww, -dv)))
        valid2 = np.array(valid2)
        # must be >99% valid
        g = np.nanmin(valid2, axis=0) > 0.99
        ww = ww[g]
        wmask = wmask[g]
        # loop on dvs and interpolate the spectrum
        for i, dv in enumerate(dvs):
            w2 = doppler(ww, -dv)
            ccfs[i, iord] = np.sum(spl(w2) * wmask)

    # sum of all per-order CCFs
    ccf = np.nansum(ccfs, axis=1)
    sig = np.sqrt(ccf)

    # oversampling ratio of the wavelength grid
    oversampling_ratio = (c / np.nanmedian(wave / np.gradient(wave, axis=1))) / np.nanmedian(np.gradient(dvs))

    # error on the RV from the photon noise
    grad = np.gradient(ccf) / np.gradient(dvs)
    # Bouchy 2001 formula
    sig_RVS[ifile] = 1 / np.sqrt(np.sum((grad / sig) ** 2)) * np.sqrt(oversampling_ratio)

    # normalize to about 1. Doesn't need to be very accurate
    ccf /= np.nanmedian(ccf)

    # valuves in the p0 and CCF fit.
    # cen, w, depth, amp, slope, amp
    p0 = [dvs[np.argmin(ccf)], 4e3, 1 - np.min(ccf), 0, 1]
    fit = curve_fit(gauss_slope, dvs, ccf, p0=p0, sigma=sig)
    print('\tRV = {:.2f} m/s'.format(fit[0][0]))
    print('\tfwhm = {:.2f} m/s'.format(fit[0][1] * 2.354))
    print('\tamp = {:.2f}'.format(fit[0][2]))
    print('\tslope = {:.2e} '.format(fit[0][3]))
    print('\tzp = 1 {} {:.2e} '.format(['-', '+'][fit[0][4] > 1], np.abs(fit[0][4] - 1)))

    # get the 0th term of the fit, this is the velocity
    RVS[ifile] = fit[0][0]
    # fit but without de gaussian. Just to normalize the CCF to unity
    fit2 = np.array(fit[0])
    # set depth to 0
    fit2[2] = 0
    # normalized CCF without the gaussian. Sets the continuum flat to 1
    ccf2 = (1 - ccf / gauss_slope(dvs, *fit2)) / fit[0][2]
    # should optionally be other values than 80% and 30% percentile of depth
    cut1 = 0.8
    # find the two points where the CCF is above the cut and interpolate with a linear fit
    lims1 = np.where(ccf2 > cut1)[0][[0, -1]]
    v1 = np.polyval(np.polyfit(ccf2[lims1[0] - 1:lims1[0] + 1], dvs[lims1[0] - 1:lims1[0] + 1], 1), cut1)
    v2 = np.polyval(np.polyfit(ccf2[lims1[1]:lims1[1] + 2], dvs[lims1[1]:lims1[1] + 2], 1), cut1)
    span80 = (v1 + v2) / 2.0
    # same at the other bisector cut
    cut1 = 0.3
    lims1 = np.where(ccf2 > cut1)[0][[0, -1]]
    v1 = np.polyval(np.polyfit(ccf2[lims1[0] - 1:lims1[0] + 1], dvs[lims1[0] - 1:lims1[0] + 1], 1), cut1)
    v2 = np.polyval(np.polyfit(ccf2[lims1[1]:lims1[1] + 2], dvs[lims1[1]:lims1[1] + 2], 1), cut1)
    span30 = (v1 + v2) / 2.0
    # difference between the two bisector cuts
    bisspan[ifile] = span80 - span30
    #
    contrasts[ifile] = fit[0][2]
    fwhm[ifile] = fit[0][1] * np.sqrt(np.log(2) * 2) * 2

plt.show()

plt.xlabel('RJD')
plt.ylabel('BIS [m/s]')
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.errorbar(rjd, bisspan, yerr=np.nanstd(bisspan), fmt='o')

plt.xlabel('RJD')
plt.ylabel('RV [m/s]')
plt.title('Proxima Centauri')
plt.tight_layout()
plt.savefig('/Users/eartigau/ccf/PROXIMA_BIS.pdf')
plt.savefig('/Users/eartigau/ccf/PROXIMA_BIS.png')
plt.show()

print(np.nanstd(RVS), np.nanmean(sig_RVS))
plt.xlabel('RJD')
plt.ylabel('RV [m/s]')
plt.grid(color='gray', linestyle='--', linewidth=0.5)
plt.errorbar(rjd, RVS, yerr=sig_RVS, fmt='o')

plt.xlabel('RJD')
plt.ylabel('RV [m/s]')
plt.title('Proxima Centauri')
plt.tight_layout()
plt.savefig('/Users/eartigau/ccf/PROXIMA_ccf.pdf')
plt.savefig('/Users/eartigau/ccf/PROXIMA_ccf.png')
plt.show()