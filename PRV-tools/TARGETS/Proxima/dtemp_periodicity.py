from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# fitting a quasiperiodic GP to the residuals
from celerite import terms, GP
from scipy.optimize import minimize


class CustomTerm(terms.Term):
    parameter_names = ("log_a", "log_b", "log_c", "log_P")

    def get_real_coefficients(self, params):
        log_a, log_b, log_c, log_P = params
        b = np.exp(log_b)
        return (
            np.exp(log_a) * (1.0 + b) / (2.0 + b), np.exp(log_c),
        )

    def get_complex_coefficients(self, params):
        log_a, log_b, log_c, log_P = params
        b = np.exp(log_b)
        return (
            np.exp(log_a) / (2.0 + b), 0.0,
            np.exp(log_c), 2*np.pi*np.exp(-log_P),
        )


def p2(t, period, t0, t1,t2,amp0,amp1,amp2):
    return amp0*np.sin(2*np.pi*(t-t0)/period)+amp1*np.sin(2*np.pi*(t-t1)/(period/2.0))+amp2*np.sin(2*np.pi*(t-t2)/(
            period/3.0))+zp

# TODO change the path to match your local path
file_path = '/Users/eartigau/lbl_NIRPS_HE/lblrdb/lbl_PROXIMA_PROXIMA.rdb'
#file_path = '/Users/eartigau/lbl_NIRPS_HE/lblrdb/lbl_GJ581_GJ581.rdb'

tbl = Table.read(file_path, format='ascii.rdb')

period = 90
t0 = 0
t1 = 0
t2=0
amp0 = 3
amp1 = 3
amp2=3
zp = 0

p0 = [period, t0, t1, t2, amp0, amp1,amp2]

fit,recon = curve_fit(p2, tbl['plot_date'], tbl['DTEMP'], p0=p0, sigma=tbl['sDTEMP'], absolute_sigma=True)

p0 = [period, t0, t1, t2, np.nanstd(tbl['d2v']), 0,0]
fit2,recon = curve_fit(p2, tbl['plot_date'], tbl['d2v'], p0=p0, sigma=tbl['sd2v'], absolute_sigma=True)

t = np.arange(tbl['plot_date'][0]-10,tbl['plot_date'][-1]+10,1)

tbl['DTEMP_residu'] = tbl['DTEMP'] - p2(tbl['plot_date'], *fit)

fig, ax = plt.subplots(nrows = 2, ncols=2,figsize=(10,5),sharex = True  )
ax[0,0].plot_date(tbl['plot_date'], tbl['DTEMP'], 'r.',alpha = 0.1)
ax[0,0].errorbar(tbl['plot_date'], tbl['DTEMP'], yerr=tbl['sDTEMP'], fmt='r.', ecolor='r', alpha=0.5)
recon2 = p2(t, *fit)
ax[0,0].plot(t, recon2, 'k--',alpha = 0.9)
ax[0,0].set(xlabel = 'Date',ylabel='DTEMP (K)',title='Proxima Centauri')
ax[0,0].grid()


ax[0,1].plot_date(tbl['plot_date'], tbl['d2v'], 'r.',alpha = 0.1)
ax[0,1].errorbar(tbl['plot_date'], tbl['d2v'], yerr=tbl['sd2v'], fmt='r.', ecolor='r', alpha=0.5)
recon3 = p2(t, *fit2)

tbl['d2v_residu'] = tbl['d2v'] - p2(tbl['plot_date'], *fit2)

ax[0,1].plot(t, recon3, 'k--',alpha = 0.9)
ax[0,1].set(xlabel = 'Date',ylabel='d2v (m2/s2)',title='Proxima Centauri')
ax[0,1].grid()

ax[1,0].plot_date(tbl['plot_date'], tbl['DTEMP_residu'], 'r.',alpha = 0.1)
ax[1,0].errorbar(tbl['plot_date'], tbl['DTEMP_residu'], yerr=tbl['sDTEMP'], fmt='r.', ecolor='r', alpha=0.5)
ax[1,0].grid()

ax[1,1].plot_date(tbl['plot_date'], tbl['d2v_residu'], 'r.',alpha = 0.1)
ax[1,1].errorbar(tbl['plot_date'], tbl['d2v_residu'], yerr=tbl['sd2v'], fmt='r.', ecolor='r', alpha=0.5)
ax[1,1].grid()


plt.show()

