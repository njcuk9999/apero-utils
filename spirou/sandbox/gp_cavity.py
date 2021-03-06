"""
Run simple GP on FP cavity
"""
import celerite
import corner
import emcee
import george
import matplotlib.pyplot as plt
import numpy as np
from celerite import terms
from george import kernels
from pandas import read_csv
from scipy.optimize import minimize

cavity_data_path = "/home/vandal/Documents/spirou/cavity_x.csv"

df = read_csv(cavity_data_path)
df.sort_values("mjdates", inplace=True, ignore_index=True)
t = df.mjdates.values
cav_len = df.cavity.values

# TODO: Loop over segments and fit each of them
# Keep only last part for test
inds = t > 58950
t = t[inds]
cav_len = cav_len[inds]

plt.plot(t, cav_len, "ko")
plt.ylabel("Cavity length [nm]")
plt.xlabel("MJD")
plt.show()

ALL_KERNELS = [
    "ExpTerm",  # Celerite Exponential (real) term
    "SqExpKernel",  # George exponential term
    "Matern52Kernel",  # George Matern 5/2 kernel
    "Matern32Kernel",  # George Matern 3/2 kernel
    "Matern32Term",  # Celerite Matern 3/2 (fast)
]

white_noise = True
george_fast = False  # Not super stable
kernel_name = "Matern32Term"
mval = np.mean(cav_len)

# Kernel parameter guess
amp = 4.5
tscale = 8.0

# Instatiate kernel
if kernel_name == "ExpTerm":
    bounds = dict(log_a=(-10, 10), log_c=(-10, 10))
    ker = terms.RealTerm(log_a=np.log(amp), log_c=np.log(tscale), bounds=bounds)
elif kernel_name == "Matern32Term":
    bounds = dict(log_rho=(-10, 10), log_sigma=(-10, 10))
    ker = terms.Matern32Term(
        log_sigma=np.log(amp), log_rho=np.log(tscale), bounds=bounds
    )
elif kernel_name == "SqExpKernel":
    ker = amp * kernels.ExpSquaredKernel(tscale)
elif kernel_name == "Matern52Kernel":
    ker = amp * kernels.Matern52Kernel(tscale)
elif kernel_name == "Matern32Kernel":
    ker = amp * kernels.Matern32Kernel(tscale)
else:
    raise ValueError("Invalid kernel")


# Setup GP model
if issubclass(type(ker), kernels.Kernel):
    gp = george.GP(
        ker,
        white_noise=np.log(2.0 ** 2) if white_noise else None,
        fit_white_noise=white_noise,
        fit_mean=True,
        solver=george.HODLRSolver if george_fast else george.BasicSolver,
    )
    wn_ind = 0
elif issubclass(type(ker), terms.Term):
    gp = celerite.GP(
        ker,
        log_white_noise=np.log(2.0 ** 2),
        fit_white_noise=True,
        mean=mval,
        fit_mean=True,
    )
    wn_ind = 2
gp.compute(t)


# Define the objective function (negative log-likelihood in this case).
def nll(p):
    gp.set_parameter_vector(p)
    ll = gp.log_likelihood(cav_len, quiet=True)
    return -ll if np.isfinite(ll) else 1e25


# And the gradient of the objective function.
def grad_nll(p):
    gp.set_parameter_vector(p)
    if isinstance(gp, celerite.GP):
        return -gp.grad_log_likelihood(cav_len, quiet=True)[1]
    else:
        return -gp.grad_log_likelihood(cav_len, quiet=True)


# And the posterior probability for emcee
def log_prob(p):
    try:
        gp.set_parameter_vector(p)
        lprior = gp.log_prior()
        if not np.isfinite(lprior):
            return -np.inf
        return gp.log_likelihood(cav_len) + lprior
    except celerite.solver.LinAlgError:
        print(f"Skipping Linalg Error with parameters {p}")
        return -np.inf


p0 = gp.get_parameter_vector()
result = minimize(nll, p0, jac=grad_nll, method="L-BFGS-B")
gp.set_parameter_vector(result.x)

# Setup mcmc
init_params = gp.get_parameter_vector()
ndim, nwalk = init_params.size, 32
p0 = init_params + 1e-5 * np.random.randn(nwalk, ndim)
sampler = emcee.EnsembleSampler(nwalk, ndim, log_prob)

p0, lp, _ = sampler.run_mcmc(p0, 500, progress=True)

sampler.reset()
sampler.run_mcmc(p0, 2000, progress=True)

fchain = sampler.get_chain(flat=True)

gp.set_parameter_vector(np.median(fchain, axis=0))

# GP prediction
x_pred = np.linspace(t.min(), t.max(), num=500)
x_pred = x_pred[~np.isin(x_pred, t)]
x_pred = np.sort(np.append(x_pred, t))
pred, pred_var = gp.predict(cav_len, x_pred, return_var=True)
pred_std = np.sqrt(pred_var)

# from astropy.stats import sigma_clip
residuals = cav_len - pred[np.isin(x_pred, t)]
print("RMS of residuals", np.std(residuals))
print("Median GP enveloppe:", np.median(pred_std))
# plt.hist(residuals)
# plt.show()

# TODO: iteratively remove worse point as long as Nsigma (84th-16th)/2
# res_rms = np.percentile

# TODO: explore GP mean error formalism ?

wn = np.exp(gp.get_parameter_vector()[wn_ind])
fig, axs = plt.subplots(2, 1, sharex="col")
axp = axs[0]
axres = axs[1]

# axp.plot(t, wl, 'k.')
# axp.errorbar(t, wl, yerr=wn, fmt="k.")
axp.plot(t, cav_len, "k.")
axp.plot(x_pred, pred, "r")
# env = np.sqrt(np.abs(pred_std ** 2 - wn ** 2))
env = pred_std
axp.fill_between(x_pred, pred - env, pred + env, color="r", alpha=0.5)
# for i in range(50):
#     mu = gp.sample_conditional(wl, x_pred)
#     axp.plot(x_pred, mu)
axp.set_ylabel("Cavity length [nm]")

axres.plot(t, cav_len - pred[np.isin(x_pred, t)], "k.")
axres.axhline(0.0, linestyle="--", color="r")
axres.set_xlabel("MJD")
axres.set_ylabel("Residuals [nm]")

plt.tight_layout()
plt.show()

names = gp.get_parameter_names()
corner.corner(sampler.flatchain, labels=names)
plt.show()

fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex="col")
samples = sampler.get_chain()
labels = gp.get_parameter_vector()
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")
plt.show()
