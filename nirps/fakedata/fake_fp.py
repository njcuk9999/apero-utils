from astropy.io import fits
import glob
import numpy as np

object_key = 'WAVE,FP,FP'
#object_key = 'CONTAM,DARK,FP'
path = '/cosmos99/nirps/raw-data/nirps_ha/2023-0[4:5]-*/*.fits'

files = glob.glob(path)

mjd_min = 60000
mjd_max = 60073

p99_min = 20

good_files = []
bad_files = []

good_im = np.zeros([4096,4096])
for i in range(len(files)):

    if 'FAKE' in files[i]:
        continue

    h = fits.getheader(files[i])
    if h['OBJECT'] != object_key:
        continue
    mjd = h['MJD-OBS']
    if mjd < mjd_min or mjd > mjd_max:
        continue

    im = np.array(fits.getdata(files[i],ext=1))
    nread = np.array(fits.getdata(files[i],ext=3))

    p99 = np.nanpercentile(im,99)

    if p99 > p99_min:
        good_files.append(files[i])
        good_im+=im
        status = 'GOOD\t'
    else:
        bad_files.append(files[i])
        status = 'BAD\t'

    p99txt = '{:5.1f}'.format(p99)
    print(status,files[i],p99txt,h['MJD-OBS'],np.max(nread))

print('*'*80)
print('BAD FILES')
print('*'*80)
for bad_file in bad_files:
    print(bad_file.split('/')[-1])
print('*'*80)

good_im/=len(good_files)

for i in range(len(bad_files)):

    with fits.open(bad_files[i]) as hdul:
        hdul[1].data = good_im
        hdul[0].header['FAKE'] = True,'Fake data because of FPs that have no flux'
        outname = bad_files[i].replace('.fits','_FAKE.fits')
        print('We write {}'.format(outname))
        hdul.writeto(outname,overwrite=True)