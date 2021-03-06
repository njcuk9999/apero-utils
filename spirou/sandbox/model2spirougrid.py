import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import rotBroad
from PyAstronomy import funcFit as fuf
from scipy.signal import convolve

c = 2.99792458e5

def epsilon(wave):
    # for a wavelength in nm, return the Epsilon value
    # of limb darkening. Here we have a place-holder law
    # It is 1.0 for 900 nm and 0.1 for 2500 nm
    #
    # it has to be between 0 and 1 for all wavelength
    

    eps = np.polyval(np.polyfit([900,2500],[1.0,0.1],1), wave)

    return eps

def get_broad_lsf(wave0, epsilon, vsini):
    # create a rotation profile
    wave_tmp = (1+np.arange(-np.int(vsini)-1,np.int(vsini)+1)/c)*wave0*10
    f_tmp = np.zeros_like(wave_tmp)
    f_tmp[np.int(vsini)] = 1

    lsf = rotBroad(wave_tmp,f_tmp,epsilon,vsini)
    lsf = lsf[lsf!=0]
    #print(lsf)
    #print(wave0)
    return lsf

def get_lsf(beta = 2.24, ew = 2.20):
    """
    wave and flux of a spectrum

    we assume that the input spectrum is a SPIRou s1d_v
    the

    LSF shape : np.exp( -0.5 * np.abs( ( wave/wave0 - 1)*c/ew  )**beta )

    """

    dvs = np.arange(-np.ceil(ew*5),np.ceil(ew*5)+1)

    lsf = np.exp(-0.5 * np.abs(dvs/ew ) ** beta)

    lsf /= np.sum(lsf)

    return dvs, lsf

"""
 please download a model and the corresponding wavelength grid from here:
 
http://phoenix.astro.physik.uni-goettingen.de/?page_id=15

    these are expressed in energy while SPIRou spectra are expressed in
    photons, do not forget to multiply the model by lambda to get a spectrum
    proportional to photons
 
"""



# read a spirou template to get the wavelength grid. This has to be a _v_ file
wave_spirou = np.array(Table.read('Template_s1d_Gl699_sc1d_v_file_AB.fits')['wavelength'])

# we read the model wavelength grid
# Careful with the various definition of wavelength, SPIRou uses
# nm and many models use Angs. Here we added a /10 factor to convert to
# nm.
wave_model = fits.getdata('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')/10.0
flux_model = fits.getdata('lte03000-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')

# velocity parameters
dv0 = 1.0 # velocity difference in km/s
vsini = 13.0 # in km/s

# LSF parameters
beta = 2.24 # shape factor exponent
ew = 2.20 # e width factor

# number of wavelength samplings you want for your limb darkening law
nlimb = 10 # 10 seems like a fair number

s1d = interp1d(wave_model,flux_model)

flux_interpol = np.zeros_like(wave_spirou)

dvs, lsf = get_lsf(beta = beta, ew = ew)


wave_limb = np.arange(np.min(wave_spirou), np.max(wave_spirou),np.max(wave_spirou)/nlimb)
wave_limb = np.append(wave_limb, np.max(wave_spirou))

lsfs = []
for i in range(len(wave_limb)):
    # get the stellar broadening LSF
    lsf2 = get_broad_lsf(wave_limb[i], epsilon(wave_limb[i]), vsini)

    # get the combined instrumental and broadening LSF
    lsf_full = np.convolve(lsf,lsf2,mode = 'full')
    lsf_full /= np.nansum(lsf_full)

    lsfs.append(lsf_full)


lsfs = np.array(lsfs)

# get the corresponding velocity grid
dvs = np.arange(len(lsf_full), dtype=float)
dvs -= np.mean(dvs)

# perform the convolution
for i in range(len(dvs)):
    print('Dv bin {0} in {1} of lsf, normalisation = {2}'.format(i+1,len(dvs),lsf_full[i]))

    lsf_int = interp1d(wave_limb, lsfs[:,i])

    flux_interpol += s1d(wave_spirou*(1+(dvs[i]-dv0)/c))*lsf_int(wave_spirou)


xrange = [1800,1805]
# plot the resulting convolved spectrum
plt.plot(wave_model,flux_model,'k-',label = 'input spectrum')
plt.xlim(xrange)
g = (wave_model > xrange[0])*(wave_model < xrange[1])
plt.ylim([np.min(flux_model[g]),np.max(flux_model[g])])

plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux (arbitrary)')
plt.plot(wave_spirou,flux_interpol,'r-', label = 'convolved spectrum')
plt.title('star : vsini={0}km/s, epsilon={1}, dv={2}km/s\nLSF: beta={3}, ew={4}'.format(vsini, epsilon, dv0, beta, ew ))
plt.legend()
plt.savefig('comp_convolve_sp.pdf')
plt.show()
