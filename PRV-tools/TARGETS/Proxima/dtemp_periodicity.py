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

t = np.arange(tbl['plot_date'][0]-10,tbl['plot_date'][-1]+10,1)

tbl['DTEMP_residu'] = tbl['DTEMP'] - p2(tbl['plot_date'], *fit)

fig, ax = plt.subplots(nrows = 2, ncols=1,figsize=(10,5),sharex = True  )
ax[0].plot_date(tbl['plot_date'], tbl['DTEMP'], 'r.',alpha = 0.1)
ax[0].errorbar(tbl['plot_date'], tbl['DTEMP'], yerr=tbl['sDTEMP'], fmt='r.', ecolor='r', alpha=0.5)
recon2 = p2(t, *fit)
ax[0].plot(t, recon2, 'k--',alpha = 0.9)
ax[0].set(xlabel = 'Date',ylabel='DTEMP (K)',title='Proxima Centauri')
ax[0].grid()



"""

# Define the kernel
bounds = dict(log_a=(None, None), log_b=(None, 5.0), log_c=(3,10),
              log_P=(1.0,3.0))
kernel = CustomTerm(log_a=0.1, log_b=0.5, log_c=5, log_P=2.0,
                    bounds=bounds)

# Define the GP
gp = GP(kernel, mean=np.mean(tbl['DTEMP_residu']))
gp.compute(tbl['plot_date'], tbl['sDTEMP'])

# Define the objective function (negative log-likelihood in this case).
def neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)

# And the gradient of the objective function.
def grad_neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.grad_log_likelihood(y)[1]

# You need to compute the GP once before starting the optimization.
gp.compute(tbl['plot_date'], tbl['sDTEMP'])

# Print the initial ln-likelihood.
print(gp.log_likelihood(tbl['DTEMP_residu']))

# Run the optimization routine.
p0 = gp.get_parameter_vector()

results = minimize(neg_log_like, p0, jac=grad_neg_log_like, method='Nelder-Mead', args=(tbl['DTEMP_residu'], gp))

# Update the kernel and print the final log-likelihood.
gp.set_parameter_vector(results.x)
print(gp.log_likelihood(tbl['DTEMP_residu']))

# Plot the data.
x = np.linspace(tbl['plot_date'][0], tbl['plot_date'][-1], 500)
mu, var = gp.predict(tbl['DTEMP_residu'], x, return_var=True)
std = np.sqrt(var)

plt.fill_between(x, mu+std, mu-std, color="g", alpha=0.5)
plt.plot(x, mu, color="g", alpha=0.5)
plt.errorbar(tbl['plot_date'], tbl['DTEMP_residu'], yerr=tbl['sDTEMP'], fmt='r.', ecolor='r', alpha=0.5)
plt.xlabel('Date')
plt.ylabel('DTEMP (K)')
plt.title('Proxima Centauri')
plt.grid()
plt.show()


"""
ax[1].plot_date(tbl['plot_date'], tbl['DTEMP_residu'], 'r.',alpha = 0.1)
ax[1].errorbar(tbl['plot_date'], tbl['DTEMP_residu'], yerr=tbl['sDTEMP'], fmt='r.', ecolor='r', alpha=0.5)
ax[1].grid()

plt.show()

