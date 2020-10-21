import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from fits2wave import fits2wave
from scipy.interpolate import InterpolatedUnivariateSpline as ius
import os
import warnings
#from sigma import *
from scipy.signal import medfilt
from tqdm import tqdm

def lowpassfilter(v,w):
    # low-pass filtering of vector v using a running median of width w. This code differs from medfilt in proprely
    # handling NaNs
    v3 = []
    i3 = []
    #
    for i in range(0,len(v),w//3):

        i0 = i-w//2
        i1 = i+w//2

        if i0 <0:
            i0 = 0
        if i1 > (len(v)):
            i1 = len(v)
        vv =v[i0:i1]

        if np.max(np.isfinite(vv)):
            v3.append(np.nanmedian(vv))
            i3.append(i)

    v3 = np.array(v3)
    i3 = np.array(i3)

    if len(i3)<5:
        return np.zeros_like(v)+np.nan

    spline = ius(i3,v3, k=2, ext=3)
    return spline(np.arange(len(v)))


def template_rv(rv,template_spline,nanmask):
    v = np.arange(len(nanmask))

    tmp = template_spline(v+rv)
    tmp[np.roll(nanmask,-int(rv))==False] = np.nan

    return tmp

# velocity of light, could be retrieved from DRS constants
c = 299792.4580
# object name
DRSOBJN = 'GL699'

if False:
    # where to search for s1d files
    path_to_s1d = '/Volumes/courlan/xdrv/s1d'

    # get s1d files
    files = np.array(glob.glob(path_to_s1d + '/2[4-5]*s1d_v*.fits'))

    os.system('mkdir {0}'.format(DRSOBJN))
    # checking all files to find the proper ones
    for i in tqdm(range(len(files))):
        hdr = fits.getheader(files[i],ext=1)
        if hdr['DRSOBJN'] == DRSOBJN:
            cmd = 'ln -s  {1} {0}/'.format(DRSOBJN,files[i])
            os.system(cmd)

files = np.array(glob.glob('{0}/*s1d_v*.fits'.format(DRSOBJN)))
files = files[np.argsort(files)]

mjdates = np.zeros_like(files,dtype = float)
for i in tqdm(range(len(files))):
    h = fits.getheader(files[i], ext=1)
    mjdates[i] = h['MJDATE']



# used to shift the mask to the right position in the spectrum
systemic = -110.0


# get template, should be capitalized as for the DRSOBJN
path_to_template = '/Volumes/courlan/xdrv/templates/Template_s1d_{0}_sc1d_v_file_AB.fits'.format(DRSOBJN)

# proper mjdates
mjd1 = 58550
mjd2 = 59020

# CCF mask of lines in the spectrum
mask = (np.loadtxt('/Volumes/courlan/xdrv/masks/{0}_neg.mas'.format(DRSOBJN)))[:,0]

# keep only YJHK bands
keep = np.zeros_like(mask,dtype = bool)
keep[(mask >938.600)*(mask<1113.400)] = True # Y band
keep[(mask >1153.586)*(mask<1354.422)] = True # J band
keep[(mask >1462.897)*(mask<1808.544)] = True # H band
keep[(mask >1957.792)*(mask<2343.105)] = True # K band

mask = mask[keep]


# read template
template = Table.read(path_to_template)
w = template['wavelength']
f = template['flux']

# find pixels that are valid for the spline
nanmask = np.isfinite(f)
v = np.arange(len(template))
template_spline = ius(v[nanmask],f[nanmask],k=5,ext=1)



# removing files with the wrong object
files = files[files != '']

# keep only the files within the proper date range
keep = (mjdates>mjd1)*(mjdates<mjd2)
files = files[keep]
mjdates = mjdates[keep]


# matrix ton construct velocities
all_velocities = np.zeros([len(files),len(mask)])

for ifile in tqdm(range(len(files))):
    tbl = Table.read(files[ifile])

    h = fits.getheader(files[ifile],ext = 1)
    rv = h['BERV']
    f = np.array(tbl['flux'])
    f/=np.nanmedian(f)

    f/=lowpassfilter(f/template_rv(rv,template_spline,nanmask),101)

    residual = f - template_rv(rv,template_spline,nanmask)

    dd = template_rv(rv,template_spline,nanmask) - template_rv(rv+1e-3,template_spline,nanmask)
    mask2 = mask*(1+systemic/c)*(1-rv/c)
    ipix = np.array(ius(w, np.arange(len(w)))(mask2), dtype=int)

    ipix_min = np.zeros_like(ipix)
    for i in range(1,len(ipix)-1):
        ipix_min[i] = ipix[i] - (ipix[i]-ipix[i-1])//2

    for i in range(len(mask)-1):
        dd_segment = dd[ipix_min[i]:ipix_min[i+1]+1]
        residual_segment = residual[ipix_min[i]:ipix_min[i+1]+1]

        all_velocities[ifile,i] = np.sum(residual_segment*dd_segment)/np.sum(dd_segment**2)

mjdates = np.array(mjdates)

keep =  np.isfinite(np.sum(all_velocities,axis=0))
all_velocities2 = all_velocities[:,keep]
mask2 = mask[keep]

rms = np.nanmedian(np.abs(all_velocities2),axis=0)

rms[np.isfinite(rms)==0] = np.inf
rms[rms<=10] = np.inf
w = 1/rms**2

w[np.isfinite(w) == 0] = 0
w/=np.nansum(w)

w2d = np.zeros_like(all_velocities2)

for i in range(len(mask2)):
    w2d[:,i]=w[i]

plt.plot(mjdates,np.nansum(all_velocities2*w2d,axis=1),'g.')
plt.show()
