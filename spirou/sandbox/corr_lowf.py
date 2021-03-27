import numpy as np
from astropy.io import fits
import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import etienne_tools as et

# get the files for TRAPPIST-1
files = np.array(glob.glob('e2ds/24*o_pp_e2dsff_tcorr_AB.fits'))

# create the cube that holds the high frequencies
highf = np.zeros([4088*49,len(files)])
# create the cube that holds the low frequencies
lowf = np.zeros([4088*49,len(files)])


dates = np.zeros(len(files),dtype = 'U99')

for i in tqdm(range(len(files))):
    # read info from file
    im,hdr = fits.getdata(files[i],header = True)
    dates[i] = hdr['DATE'].split('T')[0]
    # doppler shift the data
    im2 = et.doppler_shift( et.fits2wave(hdr),im, hdr['BERV']*1000)

    # put low and high frequencies in their respective cubes
    lowf[:, i] = et.lowpassfilter(im2.ravel())
    highf[:, i] = im2.ravel() - lowf[:,i]

# get median HP spectrum
med_highf = np.nanmedian(highf,axis=1)

# find the respective amplitudes of high-frequencies
amps = np.zeros(highf.shape[1])
prev_date =''
color = 'green'
for i in np.arange(highf.shape[1]):
    amps[i] = et.get_ratio(highf[:, i], med_highf)
    lowf[:,i]/=amps[i]
    prev_date = dates[i]

    hdr = fits.getheader(files[i])
    plt.plot(hdr['MJDATE'],amps[i],'o',color=color)

plt.xlabel('MJDATE')
plt.ylabel('Amplitude high-f')
plt.show()

# keep in a box the mean low frequency spectrum
lowf_mean = np.reshape(np.nanmean(lowf,axis=1),im.shape)

fig,ax = plt.subplots(nrows = 2, ncols=1,sharex = True,sharey = True)
for i in tqdm(range(len(files))):
    # we read again the image
    im,hdr = fits.getdata(files[i],header = True)
    # we set the amplitude of the image right
    im/=amps[i]

    wave = et.fits2wave(hdr)
    ax[0].plot(wave[35],im[35],alpha =0.2)


    # we keep only the high frequencies
    im = np.reshape(im.ravel()-et.lowpassfilter(im.ravel()),im.shape)
    # we add back the low frequencies but shift  backward as they are
    # expressed at BERV =0
    im += et.doppler_shift( et.fits2wave(hdr),lowf_mean, -hdr['BERV']*1000)

    ax[1].plot(wave[35],im[35],alpha =0.2)

    outname = files[i].split('/')
    outname[-1]= '_corrlowf.'.join(outname[-1].split('.'))
    fits.writeto('/'.join(outname),im,hdr,overwrite = True)

ax[0].set(xlabel = 'Wavelength',title='Before low-f match')
ax[1].set(xlabel = 'Wavelength',title='After low-f match')
plt.tight_layout()
plt.show()