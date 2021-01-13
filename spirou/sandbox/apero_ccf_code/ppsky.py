#import numpy as np
#from astropy.io import fits
#from scipy.signal import medfilt
import glob
from fits2wave import *
from sigma import *
#import matplotlib.pyplot as plt
#from scipy.ndimage.morphology import binary_dilation
from sklearn.decomposition import PCA
#from scipy.interpolate import InterpolatedUnivariateSpline
from wave2wave import *
from low_pass_filter import *
from tqdm import tqdm

# use the master wavelength grid as a reference, all spectra are shifted back to
# this grid for the construction of principal components
master_wave = fits.getdata('MASTER_WAVE_SPIROU.fits')

# list of OH lines
oh_lines = np.genfromtxt('list_OH_v2.0.dat')
intensities = oh_lines[:,1]
oh_lines = oh_lines[:,0]/1e1

oh_lines = oh_lines[intensities>1e-3*np.nanmax(intensities)]

# sky files to be used for the pca construction
files = glob.glob('sky_files/2*o_pp_e2dsff_AB.fits')

cube = np.zeros([49*4088,len(files)])

# constructing an enormous cube of all files
for i in tqdm(np.arange(len(files))):
    im, hdr = fits.getdata(files[i], header = True)

    im = wave2wave(im, fits2wave(hdr), master_wave)

    for iord in range(49):
        im[iord] -= lowpassfilter(im[iord],width = 101)

    cube[:, i] = im.ravel()

# replacing NaNs with zeros
cube[np.isfinite(cube) == False]=0

# finding the RMS to the median to know if it is  true line
med = np.nanmedian(cube,axis=1)
tmp = np.array(cube)
for i in range(len(files)):
    tmp[:,i] -= med
rms = np.nanmedian(np.abs(tmp))

# find likelihood of a point being a true sky emission
nsig = med/rms
p1 = np.exp(-0.5*(nsig)**2)
p2 = 1-p1

pp = p2/(p1+p2)

# set to zero out-of-line regions
for i in np.arange(len(files)):
    tmp = cube[:,i]
    tmp = tmp*pp
    tmp[tmp<0] = 0
    cube[:, i] = tmp

# constructs components for each photometric bandpass
# Y, J, H and K limits
cuts = [0,1150,1400,1900,9999]

# Number of components
n_components = 9

# puts the components for each band end-to-end
out = np.zeros( [cube.shape[0], n_components*4+1])
out[:,0] = master_wave.ravel()

# loop through bands and fill the big PCA cube
for ite in range(len(cuts)-1):
    print(ite)
    band = (master_wave.ravel() > cuts[ite]) * (master_wave.ravel() < cuts[ite+1])


    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(cube[band,:])

    for i in range(n_components):
        principalComponents[:,i] -= principalComponents[0,i]

    out[band,1+ite*n_components:1+(ite+1)*n_components] = principalComponents


# save sky files
fits.writeto('sky_PCs.fits',out,overwrite = True)