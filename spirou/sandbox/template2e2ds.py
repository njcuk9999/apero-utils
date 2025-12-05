import numpy as np
import glob
from astropy.io import fits
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import etienne_tools as et

# create a fake e2ds from the template

obj_sci = 'GL699'
obj_template = 'GL699'
science_path = 'tellurics/'
ref_blaze_file = '2498F798T802f_pp_blaze_AB.fits'
template_path = 'templates/'

template_file = template_path+'Template_s1d_'+obj_template+'_sc1d_v_file_AB.fits'

template = fits.getdata(template_file)

scifiles = glob.glob(science_path+obj_sci+'/2*tcorr*AB.fits')

# sort files by name so that they are consecutive in time
scifiles.sort()

all_durations = []

# we create the spline of the template to be used everywhere further down
valid = np.isfinite(template['flux'])
bl = fits.getdata(ref_blaze_file)
for i in range(bl.shape[0]):
    bl[i]/=np.nanpercentile(bl[i],90)

# get info on template systvel for splining correctly

print('defining all the splines required later')
flux = np.array(template['flux'])

# template removed from its systemic velocity, the spline is redefined for good
#
# Flux, 1st and 2nd derivatives. We could add more derivatives if need be one day
# Be careful to update the 'valid' mask above
spline = ius(template['wavelength'][valid],flux[valid],k=3,ext=1)
spline_mask = ius(template['wavelength'], np.isfinite(template['flux']), k=1, ext=1)

# flag to take a completely new RV measurement
failed_convergence = True
for ifile in range(len(scifiles)):

    # get the science file info
    sp, hdr = fits.getdata(scifiles[ifile], header=True)
    wave = et.fits2wave(hdr)
    model = np.zeros_like(sp)
    model_mask = np.zeros_like(sp)

    for ii in range(sp.shape[0]):
        # RV shift the spline and give it the shape of the model
        model[ii] = spline(et.doppler(wave[ii],-hdr['BERV']*1000) ) * bl[ii]
        model_mask[ii] = spline_mask(et.doppler(wave[ii], -hdr['BERV']*1000))
    model[model_mask<0.99] = np.nan

    for ii in range(sp.shape[0]):
        model[ii] = model[ii] * np.nanmedian(sp[ii]) / np.nanmedian(model[ii])

    outname = 'fake'.join(scifiles[ifile].split('tcorr'))
    print('we write {0}'.format(outname))
    model[~np.isfinite(sp)] = np.nan
    fits.writeto(outname, model, hdr, overwrite = True)