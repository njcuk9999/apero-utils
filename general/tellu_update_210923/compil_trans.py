import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import etienne_tools as et
from tqdm import tqdm

# all trans files ever observed
files = np.array(glob.glob('*trans_AB.fits'))

# you may update TELLUP_TRANS_THRESH to about -2 As
# the detrending is done on a per-pixel basis. A strong
# local bias will not affect of other bits of the
# domain.

# create a cube with all e2ds from transfiles
trans_cube = np.zeros([49,4088,len(files)])

# water and dry exponents for all trans files
TLPEH2O = np.zeros(len(files))
TLPEOTR = np.zeros(len(files))
# construct the trans cube
for i in tqdm(range(len(files))):
    trans = fits.getdata(files[i])
    hdr_trans = fits.getheader(files[i])
    TLPEH2O[i] = hdr_trans['TLPEH2O']
    TLPEOTR[i] = hdr_trans['TLPEOTR']
    trans_cube[:,:,i] = np.log(trans)

# sample vectors for the reconstruction
sample = np.zeros([3,len(TLPEH2O)])
# bias level of the residual
sample[0] = 1
# water abso
sample[1] = TLPEH2O
# dry abso
sample[2] = TLPEOTR

# sigma_cut for the worst_offender
sigma_cut = 5.0
# min number of trans
min_trans_files = len(files)//10

# create the reference vectors
zero_residual = np.full_like(trans, np.nan)
TLPEH2O_residual = np.full_like(trans, np.nan)
TLPEOTR_residual = np.full_like(trans, np.nan)

for ix in tqdm(range(trans.shape[1])):
    for iord in range(trans.shape[0]):
        # get one pixel of the trans_cube for all observations
        tmp = trans_cube[iord,ix,:]

        # construct a linear model with offset and water+dry components
        if np.sum(np.isfinite(tmp)) < min_trans_files:
            continue

        # loop until no point is an outlier beyond 3 sigma
        worst_offender = np.inf

        while worst_offender >sigma_cut:
            amp, recon = et.lin_mini(tmp, sample)

            sig = (tmp-recon)/et.sigma(tmp-recon)
            worst_offender = np.nanmax(sig)

            if worst_offender>sigma_cut:
                tmp[np.nanargmax(sig)] = np.nan
            else:
                zero_residual[iord,ix] = amp[0]
                TLPEH2O_residual[iord,ix] = amp[1]
                TLPEOTR_residual[iord,ix] = amp[2]

            if np.sum(np.isfinite(tmp)) < min_trans_files:
                break

plt.plot(zero_residual.ravel(), label = 'DC component')
plt.plot(TLPEH2O_residual.ravel(), label = 'Water component')
plt.plot(TLPEOTR_residual.ravel(), label = 'Dry component')
plt.ylim([-5,5])
plt.legend()
plt.show()

"""

save 3 arrays that have the shape of the e2ds : 
zero_residual [dc offset of residuals], 
TLPEH2O_residual [linear dependency with water absorption]
TLPEOTR_residual [linear dependency with dry absorption]

table : file names, TLPEOTR, TLPEH2O, SNR035, ob
"""


"""
**** end of the bit that you do only once
"""

# random spectrum
isp = 100

# predicted residuals
model = np.exp(zero_residual+TLPEH2O_residual*TLPEH2O[isp]+TLPEOTR_residual*TLPEOTR[isp])

# just for fun, we plot the residual and model residual
model = model.ravel()
tmp = np.exp(trans_cube[:,:,isp].ravel())
model[~np.isfinite(tmp)] = np.nan

output = tmp/model
# recon is the product of model and the pre-clean transmission.

plt.plot(model,alpha = 0.5,color = 'red', label = 'model')
plt.plot(tmp,alpha = 0.5, color = 'green',label = 'residual')
plt.plot(output,alpha = 0.9,color = 'black', label = 'ratio')
plt.ylim([0,2])
plt.legend()
plt.show()




