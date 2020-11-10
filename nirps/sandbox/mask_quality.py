import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.optimize import curve_fit
import glob
from astropy.io import fits

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
masks = np.array(glob.glob('GL*neg.mas'))
masks = np.append(masks, np.array(glob.glob('masque*.mas')))

# which photometric bandpass is used for the comparison.
# We *strongly* suggest to perform this analysis one bandpass at a time. This will
# avoid differences in masks in terms of number of lines
bandpass = 'K'

# parameters for the comparison that you probably don't need to change
c = 299792.458 # speed of light
pix_scale = 2.2# km/s approx velocity step in SPIRou spectra
step = 500 # in m/s, step in the CCF. Doesn't change anything once you are well below the pix scale
drv0 = 30 # width in km/s for the CCF. Unless you have a super weird CCF, leave at ~30 km/s

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
ccf = np.zeros_like(dvs, dtype=float)


# loop through masks
for it,mask_file in enumerate(masks):
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

    # sum of weights for CCF normalization
    norm = np.sum(weight_mask)

    # construct CCF
    for i,dv in enumerate(dvs):
        ccf[i] = np.sum(spline(wave_mask*(1+dv/c))*weight_mask)/norm

    # get mean flux from mask
    flux = spline(wave_mask)

    # assume photon-counting statistics
    sigma = np.sqrt(flux)

    # point-to-point uncertainfy in the CCF from photon noise in the weighted sum
    ccf_noise =1/np.nansum(np.abs(weight_mask))*np.sqrt( np.nansum( (weight_mask*sigma)**2 )  )

    # update table
    tbl['CCF_SNR'][it] = np.nanmedian(ccf)/ccf_noise
    tbl['CCF_SNR_BOOST'][it] = 100/ccf_noise/np.sqrt(sampling_ratio)

    # expressed in m/s, we find the Q value
    Q = np.sqrt(np.nansum((np.gradient(ccf)/step)**2))

    # we find the RV accuracy in m/s
    tbl['DVRMS'][it] =  ccf_noise/Q

    # fit a gaussian to the CCF
    p0 = [dv0,1,np.nanmedian(ccf),np.min(ccf)-np.nanmedian(ccf)]
    fit, pcov = curve_fit(gauss,dvs,ccf,p0 = p0)
    fit[1] = np.abs(fit[1])

    # update the table
    tbl['SYSTEMIC_VELOCITY'][it] = fit[0]
    tbl['GAUSSIAN_WIDTH'][it] = fit[1]
    tbl['GAUSSIAN_FWHM'][it] = fit[1]*np.sqrt(np.log(2)*2)*2
    tbl['GAUSSIAN_CONTINUUM'][it] = fit[2]
    tbl['GAUSSIAN_AMPLITUDE'][it] = fit[3]
    tbl['GAUSSIAN_FIT_RMS'][it] =  np.nanstd(ccf-gauss(dvs,*fit))
    tbl['FRACT_GAUSSIAN_FIT_RMS'][it] =  np.nanstd(ccf-gauss(dvs,*fit))/fit[2]

    tbl['CCF_DEPTH'][it] = -fit[3]/fit[2]
    tbl['Q'][it] = Q

# transform dictionnary into a table
tbl = Table(tbl)

# sort from best to worst masks
tbl = tbl[np.argsort(-tbl['DVRMS'])]

plt.clf()
# some plots just for fun
fig, ax = plt.subplots(nrows = 1,ncols=2,sharey = True)

ax[0].barh(np.arange(len(tbl)),tbl['DVRMS'],height = 0.8,label = 'DVRMS',alpha = 0.7)
ax[0].set_yticks(np.arange(len(tbl)))
ax[0].set_yticklabels(tbl['MASKS'])#,step=1)
ax[0].set_xlabel('DVRMS [m/s] for median(SNR) ~ 100')
ax[0].set_ylabel('Mask')

ax[1].barh(np.arange(len(tbl)),tbl['CCF_SNR_BOOST'],height = 0.8,label = 'CCF_SNR_BOOST',alpha = 0.7)
ax[1].set_yticks(np.arange(len(tbl)))
ax[1].set_yticklabels(tbl['MASKS'])#,step=1)
ax[1].set_xlabel('CCF_SNR_BOOST')
ax[1].set_ylabel('Mask')

plt.tight_layout()
plt.savefig( file.split('.fits')[0]+'_dvrms.png')
plt.show()
plt.clf()

fig, ax = plt.subplots(nrows = 1,ncols=2,sharey = True)
ax[0].barh(np.arange(len(tbl)),tbl['GAUSSIAN_FWHM'],height = 0.8,label = 'GAUSSIAN_FWHM',alpha = 0.7)
ax[0].set_yticks(np.arange(len(tbl)))
ax[0].set_yticklabels(tbl['MASKS'])#,step=1)
ax[0].set_xlabel('GAUSSIAN_FWHM')
ax[0].set_ylabel('Mask')


ax[1].barh(np.arange(len(tbl)),tbl['CCF_DEPTH'],height = 0.8,label = 'CCF_DEPTH',alpha = 0.7)
ax[1].set_yticks(np.arange(len(tbl)))
ax[1].set_yticklabels(tbl['MASKS'])#,step=1)
ax[1].set_xlabel('CCF_DEPTH')
ax[1].set_ylabel('Mask')

plt.tight_layout()
plt.savefig( file.split('.fits')[0]+'_gaussian.png')
plt.show()
plt.clf()


syms = ['o','*','x','+']
for i in range(len(tbl)):
    plt.plot(tbl['GAUSSIAN_FWHM'][i],tbl['CCF_DEPTH'][i],syms[(i % len(syms))],label = tbl['MASKS'][i])
plt.xlabel('GAUSSIAN FWHM [pix]')
plt.ylabel('CCF DEPTH\n[fraction of continuum)')
plt.legend()
plt.tight_layout()
plt.savefig( file.split('.fits')[0]+'_gaussian2.png')

plt.show()

tbl.write(file.split('.fits')[0]+'_compare.csv',overwrite = True)
