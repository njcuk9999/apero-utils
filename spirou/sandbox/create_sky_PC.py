import numpy as np
from astropy.io import fits
from scipy.signal import medfilt
import glob
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from tqdm import tqdm
import etienne_tools as et
import os
# use the master wavelength grid as a reference, all spectra are shifted back to

master_wave = fits.getdata('MASTER_WAVE_SPIROU.fits')

oh_lines = np.genfromtxt('list_OH_v2.0.dat')
intensities = oh_lines[:,1]
oh_lines = oh_lines[:,0]/1e1

oh_lines = oh_lines[intensities>1e-3*np.nanmax(intensities)]

files = np.array(glob.glob('all_sky/2[3-9]*o_pp_e2dsff_AB.fits'))
for i in tqdm(range(len(files))):
    hdr = fits.getheader(files[i])
    if hdr['EXTSN035'] > 2:
        files[i] = ''

    if hdr['EXPTIME'] < 180: # at least 3 minutes of integration
        files[i] = ''

files = files[files != '']

cube = np.zeros([49*4088,len(files)])

rms = np.zeros(len(files))
for i in tqdm(np.arange(len(files))):

    name_hp =  '_hp.'.join(files[i].split('.'))
    import numpy as np
    from astropy.io import fits
    from scipy.signal import medfilt
    import glob
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from tqdm import tqdm
    import etienne_tools as et
    import os

    # use the master wavelength grid as a reference, all spectra are shifted back to

    master_wave = fits.getdata('MASTER_WAVE_SPIROU.fits')

    oh_lines = np.genfromtxt('list_OH_v2.0.dat')
    intensities = oh_lines[:, 1]
    oh_lines = oh_lines[:, 0] / 1e1

    oh_lines = oh_lines[intensities > 1e-3 * np.nanmax(intensities)]

    files = np.array(glob.glob('all_sky/2[3-9]*o_pp_e2dsff_AB.fits'))
    for i in tqdm(range(len(files))):
        hdr = fits.getheader(files[i])
        if hdr['EXTSN035'] > 2:
            files[i] = ''

        if hdr['EXPTIME'] < 180:  # at least 3 minutes of integration
            files[i] = ''

    files = files[files != '']

    cube = np.zeros([49 * 4088, len(files)])

    rms = np.zeros(len(files))
    for i in tqdm(np.arange(len(files))):

        name_hp = '_hp.'.join(files[i].split('.'))

        if os.path.isfile(name_hp) == False:
            im, hdr = fits.getdata(files[i], header=True)

            im = et.wave2wave(np.array(im), et.fits2wave(hdr), master_wave)

            for iord in range(49):
                if np.nanmean(master_wave[iord]) < 2000:
                    im[iord] -= et.lowpassfilter(im[iord], width=301)
                else:
                    im[iord] -= et.lowpassfilter(im[iord], width=101)
            fits.writeto(name_hp, im, hdr)

        im, hdr = fits.getdata(name_hp, header=True)

        cube[:, i] = im.ravel()

    for ite in range(2):
        med = np.nanmedian(cube, axis=1)

        for i in range(cube.shape[1]):
            amp = np.nansum(med * cube[:, i]) / np.nansum(med ** 2)
            cube[:, i] /= amp
            print(amp)
    med = np.nanmedian(cube, axis=1)

    rms = np.zeros(cube.shape[1])
    for i in range(cube.shape[1]):
        rms[i] = et.sigma(cube[:, i])

    keep = rms < (1.5 * np.nanmedian(rms))
    cube = cube[:, keep]
    rms = rms[keep]
    files = files[keep]

    fits.writeto('cube1.fits', cube, overwrite=True)

    p1 = np.nanpercentile(cube, [16, 84], axis=1)
    sig = (p1[1] - p1[0]) / 2

    fits.writeto('cube1b.fits', cube, overwrite=True)

    med[np.isfinite(med) == 0] = 0
    for i in range(cube.shape[1]):
        bad = cube[:, i] > 5 * sig
        cube[bad, i] = med[bad]
    fits.writeto('cube1c.fits', cube, overwrite=True)

    med = np.nanmedian(cube, axis=1)
    med[np.isfinite(med) == False] = 0
    p_global = 1 - np.exp(-0.5 * (np.exp(med / et.sigma(med))) ** 2)

    cube[np.isfinite(cube) == False] = 0

    for i in np.arange(len(files)):
        nsig = cube[:, i] / et.sigma(cube[:, i])
        nsig[np.isfinite(nsig) == False] = 0
        pp = 1 - np.exp(-.5 * nsig ** 2 / 3 ** 2)

        # local probability
        p = np.ones_like(pp)
        p[nsig < 0] = 1 - pp[[nsig < 0]]
        cube[:, i] *= p
        cube[:, i] *= p_global

    for i in range(0, len(files), 5):
        plt.plot(master_wave.ravel(), cube[:, i], alpha=0.6)
    plt.plot(oh_lines, np.zeros_like(oh_lines), 'r+')
    plt.show()

    fits.writeto('cube.fits', cube, overwrite=True)

    cuts = [0, 1150, 1400, 1900, 9999]

    n_components = 9
    out = np.zeros([cube.shape[0], n_components * (len(cuts) - 1) + 1])
    out[:, 0] = master_wave.ravel()

    for ite in range(len(cuts) - 1):
        print(ite)
        band = (master_wave.ravel() > cuts[ite]) * (master_wave.ravel() < cuts[ite + 1])

        cube2 = np.zeros_like(cube)
        cube2[band, :] = cube[band, :]

        pca = PCA(n_components=n_components)
        principalComponents = pca.fit_transform(cube2)

        for i in range(n_components):
            principalComponents[:, i] -= principalComponents[0, i]
        out[band, 1 + ite * n_components:1 + (ite + 1) * n_components] = principalComponents[band, :]

    fits.writeto('sky_PCs.fits', out, overwrite=True)
    if os.path.isfile(name_hp) == False:
        im, hdr = fits.getdata(files[i], header = True)

        im = et.wave2wave(np.array(im),et.fits2wave(hdr) , master_wave )

        for iord in range(49):
            if np.nanmean(master_wave[iord])<2000:
                im[iord] -= et.lowpassfilter(im[iord],width = 301)
            else:
                im[iord] -= et.lowpassfilter(im[iord],width = 101)
        fits.writeto(name_hp, im, hdr)

    im, hdr = fits.getdata(name_hp, header=True)

    cube[:, i] = im.ravel()



for ite in range(2):
    med = np.nanmedian(cube, axis=1)

    for i in range(cube.shape[1]):
        amp = np.nansum(med*cube[:,i])/np.nansum(med**2)
        cube[:,i]/=amp
        print(amp)
med = np.nanmedian(cube,axis=1)

rms = np.zeros(cube.shape[1])
for i in range(cube.shape[1]):
    rms[i] = et.sigma(cube[:,i])


keep = rms<(1.5*np.nanmedian(rms))
cube = cube[:,keep]
rms = rms[keep]
files = files[keep]

fits.writeto('cube1.fits',cube,overwrite=True)


p1 = np.nanpercentile(cube,[16,84],axis=1)
sig = (p1[1]-p1[0])/2

fits.writeto('cube1b.fits',cube,overwrite=True)

med[np.isfinite(med) == 0] = 0
for i in range(cube.shape[1]):
    bad = cube[:,i] > 5*sig
    cube[bad,i] = med[bad]
fits.writeto('cube1c.fits',cube,overwrite=True)


med = np.nanmedian(cube, axis =1)
med[np.isfinite(med) == False] = 0
p_global = 1-np.exp(-0.5*(np.exp(med/et.sigma(med)))**2)


cube[np.isfinite(cube) == False]=0

for i in np.arange(len(files)):
    nsig = cube[:, i] / et.sigma(cube[:, i])
    nsig[np.isfinite(nsig) == False] = 0
    pp = 1-np.exp(-.5 * nsig ** 2 / 3 ** 2)

    # local probability
    p = np.ones_like(pp)
    p[nsig<0] = 1-pp[[nsig<0]]
    cube[:,i]*=p
    cube[:,i]*=p_global

for i in range(0,len(files),5):
    plt.plot(master_wave.ravel(),cube[:,i],alpha = 0.6)
plt.plot(oh_lines,np.zeros_like(oh_lines),'r+')
plt.show()

fits.writeto('cube.fits',cube,overwrite=True)

cuts = [0,1150,1400,1900,9999]

n_components = 9
out = np.zeros( [cube.shape[0], n_components*(len(cuts)-1)+1])
out[:,0] = master_wave.ravel()

for ite in range(len(cuts)-1):
    print(ite)
    band = (master_wave.ravel() > cuts[ite]) * (master_wave.ravel() < cuts[ite+1])

    cube2 = np.zeros_like(cube)
    cube2[band,:] = cube[band,:]

    pca = PCA(n_components=n_components)
    principalComponents = pca.fit_transform(cube2)

    for i in range(n_components):
        principalComponents[:,i] -= principalComponents[0,i]
    out[band,1+ite*n_components:1+(ite+1)*n_components] = principalComponents[band,:]

fits.writeto('sky_PCs.fits',out,overwrite = True)
os.system('rsync -av -e "ssh  -oPort=5822" sky_PCs.fits artigau@venus.astro.umontreal.ca:/home/artigau/www/')