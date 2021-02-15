import glob
from astropy.io import fits
import numpy as np

# code to get the best LBL mask reference for any file

all_masks = np.array(glob.glob('masks/*pos.fits'))

NSPTEMPL = np.zeros(len(all_masks))
OBJTEMP = np.zeros(len(all_masks))
OBJ = np.zeros_like(all_masks)
REFS = np.zeros_like(all_masks,dtype = bool)

bad = ['GL876']

for i in range(len(all_masks)):
    OBJ[i] = all_masks[i].split('_')[0].split('/')[1]

    hdr = fits.getheader(all_masks[i])
    if 'NSPTEMPL' not in hdr.keys():
        continue
    if 'OBJTEMP' not in hdr.keys():
        continue

    OBJTEMP[i] = hdr['OBJTEMP']
    NSPTEMPL[i] = hdr['NSPTEMPL']

    REFS[i] = (NSPTEMPL[i] > 200) and ('GL' in OBJ[i]) and (OBJ[i] not in bad)


for i in range(len(all_masks)):
    if 'FP' in OBJ[i]:
        continue
    if len(OBJ[i])<3:
        continue

    print(OBJ[i]+','+OBJ[i])
    print(OBJ[i]+','+OBJ[REFS][np.argmin(np.abs(OBJTEMP[i] - OBJTEMP[REFS]) + 1e9*(OBJ[i] == OBJ[REFS]))])