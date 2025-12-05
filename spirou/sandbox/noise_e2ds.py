from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def get_noise_model(file, width_per_fiber = 12, RON_floor = 6, RON_zp = 20, dark_current_per_s = 0.015):
    # we expect a file that is TCORR and we will need
    # to have an e2dsff_AB file corresponding in the same
    # folder as well as the recon. They are needed to
    # determine the noise.
    #
    # optional params :
    #
    # RON_zp is the zero point assuming a sqrt(N) decrease, N is the number of readouts
    # RON_floor is the per-pixel noise limit at infinite N
    # dark_current_per_s : in e-/s/pixel, the mean noise level
    # width_per_fiber : the pixel width expressed in pixel per fiber, for the AB, this is times 2.

    fic = ''.join(file.split('tcorr_'))

    print('We read file {0}'.format(fic))

    sp,hdr = fits.getdata(fic,header = True)
    sp[sp<0] = 0
    sp[~np.isfinite(sp)] = 0

    recon_file =  'e2dsff_recon_'.join(fic.split('e2dsff_'))
    print('We read file {0}'.format(recon_file))

    recon = fits.getdata(recon_file)

    dark_current = dark_current_per_s*hdr['EXPTIME']
    nread = hdr['NREADS']

    # width of trace in pixels
    # is 10 pix if A, B or C and 20 pix if AB
    width = len(hdr['FIBER'])*width_per_fiber

    # readout noise is modelled as a 20e-/sqrt(N readouts) plus a 9 e- variance
    # added in quadrature to have the observed asymptotic behavior
    ron = np.sqrt(RON_zp/np.sqrt(nread)+RON_floor**2)

    # effective RON over the width of trace
    ron *= np.sqrt(width)

    # effective dark current of the width of trace
    dark_current*=width

    noise = np.sqrt( sp+dark_current+ron**2  )
    # express in the tcorr frame, scaling the noise level
    noise/=recon

    #return noise model
    return noise

noise = get_noise_model('2574599o_pp_e2dsff_AB.fits')
sp1 = fits.getdata('2574599o_pp_e2dsff_tcorr_AB.fits')
sp2 = fits.getdata('2574598o_pp_e2dsff_tcorr_AB.fits')

plt.plot(sp1[35])#-sp2[35])
plt.plot(noise[35])
plt.show()
