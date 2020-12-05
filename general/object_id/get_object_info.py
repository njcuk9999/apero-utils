import os
import glob
import numpy as np
import pandas as pd
from astropy.io import fits
import tqdm

def add_val(adict, kwd, hdul):
    try:
        adict[kwd].append(hdul[kwd])
    except KeyError:
        adict[kwd].append(np.nan)

keys = ['OBJECT',
        'OBJRV',
        'OBJTEMP',
        'MJDEND',
        'RA_DEG',
        'DEC_DEG',
	 ]
outdict = dict()
for k in keys:
    outdict[k] = []
# fnames = []
for filepath in tqdm.tqdm(glob.glob('/spirou/cfht_nights/common/raw/20*/*o.fits')):
    #hdul = fits.open(filepath)
    hdr = fits.getheader(filepath)
    for k in keys:
        #add_val(outdict, k, hdul)
        add_val(outdict, k, hdr)
    # fnames.append(os.path.basename(filepath))
df = pd.DataFrame([])
for k in keys:
    df[k] = outdict[k]
# df.to_csv('object_info.csv', index=False)
