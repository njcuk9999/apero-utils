import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.optimize import curve_fit
import glob
from astropy.io import fits
from tqdm import tqdm

def gauss(v,v0,ew,zp,amp):
    # gaussian with a constant offset. As we know that the ccfs are negative structures, amp will be negative
    return zp+amp*np.exp( -0.5*(v-v0)**2/ew**2)

"""
All masks are compared at an SNR of 100 in the bandpass of interest. Provide a template and a number of mask files
In the final plots, it is better to have <10 masks to keep things readable
"""

# template file
file = 'Template_s1d_GL699_sc1d_v_file_AB.fits'

# provide a list of masks
masks = np.array(glob.glob('GL699*neg.mas'))
masks = np.append(masks, np.array(glob.glob('masque*.mas')))

# which photometric bandpass is used for the comparison.
# We *strongly* suggest to perform this analysis one bandpass at a time. This will
# avoid differences in masks in terms of number of lines
bandpass = 'H'

# parameters for the comparison that you probably don't need to change
c = 299792.458 # speed of light
pix_scale = 2.2# km/s approx velocity step in SPIRou spectra
step = 1000 # in m/s, step in the CCF. Doesn't change anything once you are well below the pix scale
drv0 = 10 # width in km/s for the CCF. Unless you have a super weird CCF, leave at ~30 km/s

# bandpass definitions
if 'Y' == bandpass:
    wave1,wave2 = 938,1113
if 'J' == bandpass:
    wave1,wave2 = 1153, 1354
if 'H' == bandpass:
    wave1,wave2 = 1462, 1808
if 'K' == bandpass:
    wave1,wave2 = 1957, 2400

# getting the systemic velocity from header
hdr = fits.getheader(file,ext=1)
dv0 = hdr['OBJRV']

# get the template
tbl = Table.read(file)
wave_template = np.array(tbl['wavelength'])
flux_template = np.array(tbl['flux'])

# find the scaling of the spline such that we remain at an SNR~100
sampling_ratio = pix_scale/(step/1000)

# clip the template to the wavelength domain, leave 1% margin at edges. This speeds up the spline afterward
keep = (wave_template > wave1 * .99) & (wave_template < wave2 * 1.01) & np.isfinite(flux_template)
wave_template = wave_template[keep]
flux_template = flux_template[keep]

# create a spline for the CCF
spline = ius(wave_template, 10000*flux_template/np.nanmedian(flux_template)/sampling_ratio, k=5)


# create a dictionary (later transformed into a table) to hold all sorts of outputs
tbl = dict()

# mask names
tbl['MASKS'] = np.array(masks)

# trim extensions if .mass
for i in range(len(masks)):
    tbl['MASKS'][i] = np.array(tbl['MASKS'][i].split('.mas')[0])

# object's systemic velocity for that mask
tbl['SYSTEMIC_VELOCITY'] = np.zeros_like(masks, dtype = float)
# pixel width of the Gaussian fit
tbl['GAUSSIAN_WIDTH'] = np.zeros_like(masks, dtype = float)
# pixel FWHM of the Gaussian fit
tbl['GAUSSIAN_FWHM'] = np.zeros_like(masks, dtype = float)
# continuum level of the Gaussian fit
tbl['GAUSSIAN_CONTINUUM'] = np.zeros_like(masks, dtype = float)
# amplitude of the Gaussian fit
tbl['GAUSSIAN_AMPLITUDE'] = np.zeros_like(masks, dtype = float)
# depth of the CCF
tbl['CCF_DEPTH'] = np.zeros_like(masks, dtype = float)
# RMS of the residuals to the gaussian fit
tbl['GAUSSIAN_FIT_RMS'] = np.zeros_like(masks, dtype = float)
# same as above, normalized expressed as a fraction of continuum level
tbl['FRACT_GAUSSIAN_FIT_RMS'] = np.zeros_like(masks, dtype = float)

# RV content of CCF
tbl['Q'] = np.zeros_like(masks, dtype = float)
# Ratio of CCF SNR to spectrum SNR
tbl['CCF_SNR_BOOST'] = np.zeros_like(masks, dtype = float)
# Mean SNR of CCF
tbl['CCF_SNR'] = np.zeros_like(masks, dtype = float)
# RV accuracy from photon statistics alone
tbl['DVRMS'] = np.zeros_like(masks, dtype = float)

# create velocity steps for CCF
dvs = np.arange(dv0 - drv0, dv0 + drv0, step / 1000.0)



it = 0
mask_file = masks[0]
print('({0}/{1}) {2}'.format(it,len(masks),mask_file))

# load mask
mask = np.loadtxt(mask_file)

# extract columns
wave_mask = mask[:,0]
weight_mask = mask[:,2]

# clip the mask to the wavelength domain
keep = (wave_mask> wave1) & (wave_mask<wave2)
wave_mask = wave_mask[keep]
weight_mask = weight_mask[keep]

flux = spline(wave_mask)
sigma = np.sqrt(flux)

iline = 500

for iline in range(100):
    ccfs = np.zeros([len(dvs),len(weight_mask)], dtype=float)

    # construct CCF
    for i,dv in enumerate(dvs):
        ccfs[i,:] = spline(wave_mask*(1+dv/c))*weight_mask

    epsilons = []
    dvrms_all = []
    for epsilon in np.arange(-weight_mask[iline]*.1,weight_mask[iline]*.1,weight_mask[iline]/10.0):

        weight_mask2 = np.array(weight_mask)

        weight_mask2[iline] += epsilon

        ccf = np.nansum(ccfs,axis=1)+ccfs[:,iline]*epsilon

        ccf /= np.nansum(weight_mask2)

        # assume photon-counting statistics

        # point-to-point uncertainfy in the CCF from photon noise in the weighted sum
        ccf_noise = np.sqrt(np.nansum((weight_mask2 * sigma) ** 2)) / np.nansum(np.abs(weight_mask2))

        # expressed in m/s, we find the Q value
        Q = np.sqrt(np.nansum((np.gradient(ccf) / step) ** 2))

        # we find the RV accuracy in m/s
        dvrms = ccf_noise / Q

        dvrms_all.append(dvrms)
        epsilons.append(epsilon)

        weight_mask[iline]+=epsilons[np.argmin(dvrms_all)]

    #dvrms_all -= np.min(dvrms_all)
    #dvrms_all/=np.max(dvrms_all)

    #print(dvrms_all)
    print(np.min(dvrms_all))