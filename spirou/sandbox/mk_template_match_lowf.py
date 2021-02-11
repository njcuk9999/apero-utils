import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import etienne_tools as et
import glob

files = np.array(glob.glob('tellurics/TRAPPIST-1/2??????o_pp_e2dsff_tcorr_AB.fits'))#[::5]
files = files[np.argsort(files)]

i1 = 0

iord = 45

sp1 = fits.getdata(files[i1])[iord]
sp2 = np.zeros([len(files),4088])
wave = et.fits2wave(fits.getheader(files[i1]))

for i2 in range(len(files)):
    print(i2,len(files))
    sp2[i2] = et.wave2wave(files[i2],files[i1])[iord]

sp2_hf = np.array(sp2)
for i2 in range(len(files)):
    print(i2,len(files))
    sp2_hf[i2] -= et.lowpassfilter(sp2_hf[i2] )

for ite in range(3):
    ref = np.nanmedian(sp2_hf,axis=0)
    for i2 in range(len(files)):
        g = np.isfinite(sp2_hf[i2] * ref)
        amp = np.nansum(sp2_hf[i2][g] * ref[g]) / np.nansum(ref[g] ** 2)
        print(amp)
        sp2_hf[i2]/=amp
        sp2[i2]/=amp

fits.writeto('test.fits',sp2_hf,overwrite = True)
fits.writeto('test2.fits',sp2,overwrite = True)
fits.writeto('test3.fits',sp2 - sp2_hf,overwrite = True)



#g = np.isfinite(sp1_hf*sp2_hf)
#amp = np.nansum(sp2_hf[g]*sp1_hf[g])/np.nansum(sp2_hf[g]**2)
#print(amp)

i1 = 30
i2 = 103

plt.plot(wave[iord],sp2_hf[i2],alpha = 0.5,label = 'sp[{0}]'.format(i2))
plt.plot(wave[iord],ref,alpha = 0.5,label = 'model HP')
plt.plot(wave[iord],sp2_hf[i2] - ref,label = 'difference')
plt.legend()
plt.show()

plt.plot(wave[iord],sp2[i2],alpha = 0.5,label = 'sp[{0}]'.format(i2))
plt.plot(wave[iord],sp2[i1],alpha = 0.5,label = 'sp[{0}]'.format(i1))
plt.plot(wave[iord],sp2[i2] - sp2[i1],alpha = 0.5,label = 'diff = sp[{0}] - sp[{1}]'.format(i2,i1))
plt.legend()

plt.show()
