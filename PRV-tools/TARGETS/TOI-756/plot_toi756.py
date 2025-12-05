from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import radvel
import radvel.likelihood
from radvel.plot import orbit_plots, mcmc_plots
import yaml

def load_yaml(param_file):
    with open(param_file, "r") as yamlfile:
        params = yaml.load(yamlfile, Loader=yaml.FullLoader)

    return dict(params)

def plot_results(like,alpha = 0.9):
    plt.errorbar(
        like.x, like.model(t)+like.residuals(),
        yerr=like.yerr, fmt='o',alpha=alpha
        )
    plt.plot(ti, like.model(ti),alpha=alpha)
    plt.xlabel('Time')
    plt.ylabel('RV')
    plt.draw()

def initialize_model(inputs, tbl):

    time_base = 0

    nplanet =  len(inputs['planet_params'])

    params = radvel.Parameters(nplanet,basis=inputs['basis']) # number of planets

    nplanet = 1
    for idplanet in inputs['planet_params'].keys():
        print(idplanet)

        params['per{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['per'][0]) # period 1
        params['tc{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['tc'][0])
        params['secosw{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['secosw'][0])
        params['sesinw{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['sesinw'][0])
        params['k{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['k'][0])
        nplanet += 1

    if inputs['keplerian_fit']['gpfit']:
        for key in inputs['gp_params'].keys():
            params['gp_{}'.format(key)] = radvel.Parameter(value=inputs['gp_params'][key][0])

    mod = radvel.RVModel(params, time_base=time_base)

    mod.params['dvdt'] = radvel.Parameter(value=0, vary=inputs['keplerian_fit']['dvdt'])
    mod.params['curv'] = radvel.Parameter(value=0, vary=inputs['keplerian_fit']['curv'])

    t = np.array(tbl['rjd'])
    vel = np.array(tbl['vrad'])
    errvel = np.array(tbl['svrad'])

    if inputs['keplerian_fit']['gpfit']:
        hnames = ['gp_'+p for p in inputs['gp_params'].keys()]

        like = radvel.likelihood.GPLikelihood(mod,t, vel,
                                              errvel, hnames,
                                              kernel_name=inputs['keplerian_fit']['kernel_name']
                                              )

    else:
        like = radvel.likelihood.RVLikelihood(mod, t, vel, errvel)

    value, vary = inputs['keplerian_fit']['gamma']
    like.params['gamma'] = radvel.Parameter(value=value, vary=vary, linear=True)
    value, vary = inputs['keplerian_fit']['jit']
    like.params['jit'] = radvel.Parameter(value=value, vary=vary, linear=True)

    nplanet = 1
    for idplanet in inputs['planet_params'].keys():
        print(idplanet)

        like.params['secosw{}'.format(nplanet)].vary = False
        like.params['sesinw{}'.format(nplanet)].vary = False
        like.params['k{}'.format(nplanet)].vary = True
        like.params['per{}'.format(nplanet)].vary = False
        like.params['tc{}'.format(nplanet)].vary = True
        nplanet += 1


    like.params['dvdt'].vary = inputs['keplerian_fit']['dvdt']
    like.params['curv'].vary = inputs['keplerian_fit']['curv']

    post = radvel.posterior.Posterior(like)

    nplanet = 1
    for idplanet in inputs['planet_params'].keys():
        print(idplanet)
        for key in inputs['planet_params'][idplanet].keys():
            print(key)
            key2 = '{}{}'.format(key,nplanet)
            value, priortype, range, vary = inputs['planet_params'][idplanet][key]

            if priortype == 'uniform':
                post.priors += [radvel.prior.HardBounds(key2, range[0], range[1])]
            if priortype == 'Jeffreys':
                post.priors += [radvel.prior.Jeffreys(key2, range[0], range[1])]
            if priortype == 'Gaussian':
                post.priors += [radvel.prior.Gaussian(key2, value, range)]

        nplanet += 1

    if inputs['keplerian_fit']['gpfit']:
        for key in inputs['gp_params'].keys():
            value, priortype, range, vary = inputs['gp_params'][key]
            if vary:
                if priortype == 'Gaussian':
                    post.priors += [radvel.prior.Gaussian('gp_{}'.format(key), value, range)]
                if priortype == 'Jeffreys':
                    post.priors += [radvel.prior.Jeffreys('gp_{}'.format(key), range[0], range[1])]

    return mod, like, post

inputs = load_yaml('/Users/eartigau/apero-utils/PRV-tools/TARGETS/proxima.yaml')

tbl = Table.read(inputs['rdbfile'], format='ascii.rdb')
tbl['vrad'] -= np.nanmedian(tbl['vrad'])
tbl = tbl[np.abs(tbl['vrad'])<10]

t = tbl['rjd']+2400000.5

vel = tbl['vrad']
errvel = tbl['svrad']
margin = 2.0
ti = np.linspace(t[0]-margin,t[-1]+margin,500)

tt = tbl['plot_date']
ti_plot = np.linspace(tt[0]-margin,tt[-1]+margin,500)

mod, like, post = initialize_model(inputs, tbl)

print(like)
print(post)

res  = optimize.minimize(
    post.neglogprob_array,     # objective function is negative log likelihood
    post.get_vary_params(),    # initial variable parameters
    method='Nelder-Mead',           # Powell also works
    )

print('BIC : ',post.bic())

plt.errorbar(tbl['rjd'],post.residuals(),yerr = tbl['svrad'],fmt ='g.')
plt.show()

like.model(ti)

like.params['gamma'].value
plt.figure(figsize=(10,5))
plt.plot_date(tbl['plot_date'],tbl['vrad'], 'k.',alpha = 0.01)
plt.errorbar(tbl['plot_date'], tbl['vrad'], yerr=tbl['svrad'], fmt='o',color = 'r', ecolor='k', alpha=0.9)
plt.plot(ti_plot, (like.model(ti)+like.params['gamma'].value),alpha=0.9)

plt.xlabel('Date')
plt.grid(color = 'black',alpha = 0.2)
plt.ylabel('Velocity [m/s]')
#plot_results(like)             # plot best fit model
print(post)
plt.tight_layout()
plt.show()




if not inputs['keplerian_fit']['gpfit']:
    RVPlot = orbit_plots.MultipanelPlot(post, legend=True)
else:
    RVPlot = orbit_plots.GPMultipanelPlot(
        post,
        subtract_gp_mean_model=False,
        plot_likelihoods_separately=False,
        subtract_orbit_model=False
    )

RVPlot.plot_multipanel()
plt.show()
stop
#df = radvel.mcmc(post,nwalkers=20,nrun=400,savename='rawchains.h5')
chains = radvel.mcmc(post,nwalkers=inputs['mcmc_params']['nwalkers'],
                     nrun=inputs['mcmc_params']['nwalkers'],
                     ensembles=3,
                     savename=inputs['mcmc_params']['savename'])

Corner = mcmc_plots.CornerPlot(post, chains)
Corner.plot()


quants = chains.quantile([0.159, 0.5, 0.841]) # median & 1sigma limits of posterior distributions

for par in post.params.keys():
    if post.params[par].vary:
        med = quants[par][0.5]
        high = quants[par][0.841] - med
        low = med - quants[par][0.159]
        err = np.mean([high,low])
        err = radvel.utils.round_sig(err)
        med, err, errhigh = radvel.utils.sigfig(med, err)
        print('{} : {} +/- {}'.format(par, med, err))