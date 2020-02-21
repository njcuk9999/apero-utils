import numpy as np
import glob
import os
from scipy.signal import medfilt
import matplotlib.pyplot as plt
from astropy.io import fits

"""

Combine a set of DARK_FP files into master calibrations to correct the leakage between
the calibration and science channels. The DARP_FPs are read, we normalize each order (both AB and C) by the
C fiber FP value and median-combine all files into a super-DARK_FP. These files are then used to subtract the
FP leakage onto the AB channels

"""

# find files to be combined into a master DARK_FP calibration
fics = np.array(glob.glob('calib_files/2*_AB.fits'))


# will contain all the AB and C fiber images
all_AB = []
all_C = []

# kernel to smooth the spectrum in the AB fiber to that only high-frequencies are kept
# this is necessary to avoid contaminating the correction frame with the thermal background
# of individual dark_fp files

w_smooth = 15 # e-width of the smoothing kernel
# gaussian that goes from -3 to +3 e-width
ker = np.exp(-.5*(np.arange(-3*w_smooth,3*w_smooth+1)/(w_smooth) )**2)
ker = ker/np.sum(ker)

# loop  through input files
for i in range(len(fics)):
    print(i,len(fics))

    C_file = 'C'.join(fics[i].split('AB'))
    if os.path.isfile(C_file) == False:
        print('file {0} does not exist, we won''t include {1} in the combination!'.format(C_file, fics[i]))
        continue
    # AB file receiving the leakage
    cal_AB = fits.getdata(fics[i])
    # C file from which the leakage comes
    cal_C = fits.getdata( C_file)

    # loop around orders
    for iord  in range(cal_AB.shape[0]):

        # remove the pedestal from the FP to avoid an offset from thermal background
        cal_C[iord] -= np.nanpercentile(cal_C[iord],5)


        # median filtering has to be an odd number
        tmp =  medfilt(cal_AB[iord],2*(w_smooth//2)+1)
        tmp[np.isfinite(tmp) == False] = 0

        # find a proxy for the low-frequency in the AB channel
        mask = np.ones_like(tmp)
        mask[tmp==0]=0
        low_f = np.convolve(tmp, ker, mode='same') / np.convolve(mask, ker, mode='same')

        """
        # just for tests
        if iord == 30:
            plt.plot(cal_AB[iord])
            plt.plot(low_f)
            plt.show(block=True)            
        """

        # keep only high-frequencies
        cal_AB[iord] -= low_f

        # normalize both channels by the 90th percentile of the leaking fiber
        amp = np.nanpercentile(cal_C[iord],90)
        cal_AB[iord]/=amp
        cal_C[iord]/=amp

    # append onto all other AB anc C spectra
    all_AB.append(cal_AB)
    all_C.append(cal_C)

# make np arrays for the median
all_AB = np.array(all_AB)
all_C = np.array(all_C)

# save the master frames for correction
fits.writeto('master_C.fits',np.nanmedian(all_C,axis= 0),overwrite = True)
fits.writeto('master_AB.fits',np.nanmedian(all_AB,axis= 0),overwrite = True)

