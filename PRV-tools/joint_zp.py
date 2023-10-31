from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import glob
from etienne_tools import odd_ratio_mean,sigma

# import GP package George
import george
from george import kernels
from etienne_tools import sigma

files = np.array(glob.glob('/Users/eartigau/all_rdbs/spirou/lbl2_G*PCA03_drift.rdb'))

sigs = np.zeros(len(files))
ns = np.zeros(len(files))
for i in range(len(files)):
    tbl = Table.read(files[i], format='rdb')
    tbl['vrad'] -= np.nanmedian(tbl['vrad'])
    sigs[i] = sigma(tbl['vrad'] )
    ns[i] = sigma(tbl['vrad']/tbl['svrad'] )


files = files[(sigs<5)*(ns<3)]

t = np.array([])
v = np.array([])
sv = np.array([])
for f in files:
    print(f)
    tbl = Table.read(f, format='rdb')
    tbl['vrad'] -= odd_ratio_mean(tbl['vrad'],tbl['svrad'])[0]

    plt.errorbar(tbl['rjd'],tbl['vrad'],yerr = tbl['svrad'],fmt = '.',label = f.split('/')[-1],alpha = 0.2)
    sig = sigma(tbl['vrad'] - np.roll(tbl['vrad'], 1))/np.sqrt(2)
    sig2 = np.nanmedian(tbl['svrad'])

    if sig>sig2:
        sigg = np.sqrt(sig**2-sig2**2)
    else:
        sigg = 0

    t = np.append(t,tbl['rjd'])
    v = np.append(v,tbl['vrad'])
    sv = np.append(sv,np.sqrt(tbl['svrad']**2+sigg**2))


valid = np.where(np.isfinite(v))
t = t[valid]
v = v[valid]
sv = sv[valid]

kernel = np.var(v) * kernels.ExpSquaredKernel(50)
gp = george.GP(kernel)
gp.compute(t, sv)
x_pred = np.linspace(np.min(t), np.max(t), 1200)
pred, pred_var = gp.predict(v, x_pred, return_var=True)


pred2, pred_var2 = gp.predict(v, t, return_var=True)

plt.fill_between(x_pred, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var),
                color="k", alpha=0.2)
plt.legend()
plt.show()
