from astropy.table import Table
from scipy import optimize
from radvel.plot import orbit_plots, mcmc_plots
import corner
import matplotlib.pyplot as plt
import numpy as np
import os
import radvel
import radvel.likelihood
import yaml
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from tqdm import tqdm
from astropy.table import Table

def unit(key):
    if key.upper()=='VRAD':
        return 'm/s'
    if key.upper()=='SVRAD':
        return 'm/s'
    if key.upper()=='DTEMP':
        return 'K'
    if key.upper()=='SVRAD':
        return 'K'
    if key.upper()=='D2V':
        return 'm$^2$/s$^2$'
    if key.upper()=='SD2V':
        return 'm$^2$/s$^2$'
    if key.upper()=='FWHM':
        return 'm/s'
    if key.upper()=='SIG_FWHM':
        return 'm/s'
    if key.upper()=='FWHM':
        return 'm/s'
    if key.upper()=='SFWHM':
        return 'm/s'
    if key.upper()=='BERV':
        return 'km/s'
    if key.upper()=='CONTRAST':
        return 'normalized depth'
    return ''

def get_plot_time(t,timestep = 0.1):
    margin = (t[-1]-t[0])/20.0
    ti = np.linspace(t[0]-margin,t[-1]+margin,int(np.ceil(t[-1]-t[0])/timestep))

    return ti

def bin_table(tbl):
    # Pass a table and bin per day
    # returns a table with the binned values
    # for non-float values, it takes the first value
    # for float values, it takes the weighted average
    tbl2 = Table()
    udate = np.unique(np.floor(tbl['rjd']))

    for key in tbl.keys():
        tbl2[key] = np.zeros_like(udate,dtype = tbl[key].dtype)

    for i in tqdm(np.arange(len(udate)), leave = False):
        g = np.where(np.floor(tbl['rjd']) == udate[i])[0]

        w = 1./tbl['svrad'][g]**2

        for key in tbl.keys():
            if tbl[key].dtype == float:
                tbl2[key][i] = np.average(tbl[key][g], weights = w)
            else:
                tbl2[key][i] = tbl[key][g[0]]

        tbl2['svrad'][i] = np.sqrt(1./np.sum(w))

    return tbl2

def load_yaml(param_file):
    # load the yaml file
    # return a dictionary with the inputs

    with open(param_file, "r") as yamlfile:
        inputs = yaml.load(yamlfile, Loader=yaml.FullLoader)

    inputs = dict(inputs)

    # adding all keys that are absent but needed as False
    emptypar = {'dvdt':False,'curv':False,'jit':False,'gpfit':False,'fit_berv_gp':False}
    for key in emptypar.keys():
        if key not in inputs['keplerian_fit'].keys():
            inputs['keplerian_fit'][key] = emptypar[key]

    return inputs

def initialize_model(inputs, tbl2, key_quantity = 'vrad', key_x_value = 'jd'):
    # initialize the model
    # inputs: dictionary with the inputs
    # tbl: table with the data

    tbl = Table(tbl2)

    time_base = 0

    if 'planet_params' in inputs.keys():
        # find the number of planets to initialize the model
        nplanet =  len(inputs['planet_params'])-1

        params = radvel.Parameters(nplanet,basis=inputs['planet_params']['basis']) # number of planets

        nplanet = 1 # counter of planet number
        for idplanet in inputs['planet_params'].keys():
            if idplanet == 'basis':
                continue

            print('ID : {}'.format(idplanet))
            # initialize relevant parameters
            params['per{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['per'][0]) # period 1
            params['tc{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['tc'][0])
            params['secosw{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['secosw'][0])
            params['sesinw{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['sesinw'][0])
            params['k{}'.format(nplanet)] = radvel.Parameter(value= inputs['planet_params'][idplanet]['k'][0])
            nplanet += 1
    else:
        # no planet in model
        params = radvel.Parameters(0)

    keplerian_only = False
    # conditions to fit only keplerians
    if inputs['activity_correction']['type'] == 'None':
        keplerian_only = True

    if inputs['activity_correction']['type'] == 'Detrend' and key_quantity=='vrad':
        keplerian_only = True

    if inputs['activity_correction']['type'] == 'BERV' and key_x_value!='BERV':
        keplerian_only = True


    # initialize the GP parameters
    if not keplerian_only:
        for key in  inputs['activity_correction']['gp_params'].keys():
            if key == 'kernel':
                continue
            params['gp_{}'.format(key)] = radvel.Parameter(value=inputs['activity_correction']['gp_params'][key][0])

    mod = radvel.RVModel(params, time_base=time_base)

    tbl = tbl[np.argsort(tbl[key_x_value])]

    t_gp = np.array(tbl[key_x_value])
    vel = np.array(tbl[key_quantity])
    errvel = np.array(tbl['s'+key_quantity])


    if inputs['keplerian_fit']['dvdt']:
        fit = np.polyfit(t_gp-np.mean(t_gp),vel,1,w=1/errvel**2)
        print('polyfit : {} '.format(fit))
        dvdt,_ = fit
        curv = 0.0
    else:
        curv, dvdt = 0.0,0.0

    mod.params['dvdt'] = radvel.Parameter(value=dvdt, vary=inputs['keplerian_fit']['dvdt'])
    mod.params['curv'] = radvel.Parameter(value=curv, vary=inputs['keplerian_fit']['curv'])

    if keplerian_only:
        print('We only fit keplerians')
        like = radvel.likelihood.RVLikelihood(mod, t_gp, vel, errvel)
    else:
        hnames = ['gp_'+p for p in inputs['activity_correction']['gp_params'].keys()]
        hnames.remove('gp_kernel')

        like = radvel.likelihood.GPLikelihood(mod,t_gp, vel,
                                              errvel, hnames,
                                              kernel_name=inputs['activity_correction']['gp_params']['kernel']
                                              )

    # fit the gamma term
    value, vary = inputs['keplerian_fit']['gamma']
    value = np.nanmedian(vel)
    print('gamma = {}, allowed to vary : {}'.format(value,vary))
    like.params['gamma'] = radvel.Parameter(value=value, vary=vary, linear=True)


    value, vary = inputs['keplerian_fit']['jit']
    print('jitter = {}, allowed to vary : {}'.format(value,vary))
    like.params['jit'] = radvel.Parameter(value=value, vary=vary, linear=True)

    # addition of planets in the merit function
    if 'planet_params' in inputs.keys():
        nplanet = 1
        for idplanet in inputs['planet_params'].keys():
            if idplanet == 'basis':
                continue

            print(idplanet)

            for key in inputs['planet_params'][idplanet].keys():
                print(key, inputs['planet_params'][idplanet][key][3])
                like.params['{}{}'.format(key,nplanet)].vary = inputs['planet_params'][idplanet][key][3]
            nplanet += 1

    if not keplerian_only:
        for key in inputs['activity_correction']['gp_params'].keys():
            if key == 'kernel':
                continue
            like.params['gp_{}'.format(key)].vary = inputs['activity_correction']['gp_params'][key][3]

    # allow for variation in dvdt and curv if requested
    like.params['dvdt'].vary = inputs['keplerian_fit']['dvdt']
    like.params['curv'].vary = inputs['keplerian_fit']['curv']

    #
    post = radvel.posterior.Posterior(like)

    value, vary = inputs['keplerian_fit']['jit']
    if vary:
        post.priors += [radvel.prior.Jeffreys('jit', 0.1*value, 10.*value)]

    if 'planet_params' in inputs.keys():
        nplanet = 1
        for idplanet in inputs['planet_params'].keys():
            if idplanet == 'basis':
                continue

            print(idplanet)
            for key in inputs['planet_params'][idplanet].keys():
                print(key)
                key2 = '{}{}'.format(key,nplanet)
                value, priortype, valrange, vary = inputs['planet_params'][idplanet][key]

                if priortype == 'Uniform':
                    post.priors += [radvel.prior.HardBounds(key2, valrange[0], valrange[1])]
                if priortype == 'Jeffreys':
                    post.priors += [radvel.prior.Jeffreys(key2, valrange[0], valrange[1])]
                if priortype == 'Gaussian':
                    post.priors += [radvel.prior.Gaussian(key2, value, valrange)]

            nplanet += 1

    if not keplerian_only:
        for key in inputs['activity_correction']['gp_params'].keys():
            if key == 'kernel':
                continue

            value, priortype, valrange, vary = inputs['activity_correction']['gp_params'][key]
            print('gp_{}'.format(key), value, priortype, valrange, vary)

            if vary:
                if priortype == 'Gaussian':
                    post.priors += [radvel.prior.Gaussian('gp_{}'.format(key), value, valrange)]
                if priortype == 'Jeffreys':
                    post.priors += [radvel.prior.Jeffreys('gp_{}'.format(key), valrange[0], valrange[1])]

    return mod, like, post

def bervfit(time, time1, amp1, time2, amp2, zp):
    lenght_of_year = 365.2425
    return np.sin( 2*np.pi*(time-time1)/lenght_of_year ) * amp1 + np.sin( 4*np.pi*(time-time2)/lenght_of_year ) * amp2+zp

def approx_berv(t,tbl,timekey = 'jd'):

    from scipy.optimize import curve_fit

    tt = tbl[timekey]
    fit,_ = curve_fit(bervfit, tt, tbl['BERV'], p0=[0,0,0,0,0])

    return bervfit(t,*fit)

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.timeseries import LombScargle
from astropy.stats import sigma_clip

def find_nearest(array,value):
    index=(np.abs(array-value)).argmin()
    return index

def BIC(k, n, log_L):
    return k * np.log(n) - 2 * log_L

def sigma(tmp):
    # return a robust estimate of 1 sigma
    sig1 = 0.682689492137086
    p1 = (1-(1-sig1)/2)*100
    return (np.nanpercentile(tmp,p1) - np.nanpercentile(tmp,100-p1))/2.0

def periodogram(tbl):
    # Remove RV outliers
    rv = sigma_clip(tbl['vrad'].data, sigma = 3)
    index_mask = rv.mask == True
    tbl = tbl[~index_mask]
    print('BAD RV = ', sum(index_mask))

    bjd, rv, erv, dTEMP, edTEMP = tbl['rjd'].data + 2400000, tbl['vrad'].data, tbl['svrad'].data, tbl['DTEMP'].data, tbl['sDTEMP'].data

    # Inspect RV curve and second derivative

    fig, ax = plt.subplots(2, 2, figsize = (16,8))

    ax[0,0].set_title('RV', fontsize = 14)
    ax[0,0].errorbar(bjd - 2400000, rv - np.nanmedian(rv), yerr = erv, color = 'k', linestyle = "none", marker = "o", alpha = 0.75, capsize = 2, zorder = 1, label = r'$\sigma_{\rm RV}$ = ' + '{0:.2f} m/s\nRMS = {1:.2f} m/s'.format(np.median(erv), sigma(rv)))
    ax[0,0].set_xlabel("BJD - 2400000", fontsize = 14)
    ax[0,0].set_ylabel("RV (m/s)", fontsize = 14)
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(axis="both", labelsize=12, direction = 'in', length = 5)
    ax[0,0].tick_params(axis="both", which = 'minor', direction = 'in', length = 3)
    ax[0,0].yaxis.set_ticks_position('both')
    ax[0,0].xaxis.set_ticks_position('both')
    ax[0,0].legend(loc = 'best', fontsize = 12)

    ls = LombScargle(bjd, rv, erv)
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[0,1].plot(1 / frequency, power, color="k", label = "Lomb-Scargle periodogram")
    ax[0,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(bjd, np.ones(len(bjd)), fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[0,1].plot(1/frequency, power, color="r", label = "Window function")
    ax[0,1].set_xlabel("Period (days)", fontsize=14)
    ax[0,1].set_ylabel("Power", fontsize=14)
    ax[0,1].set_xlim([1,1000])
    ax[0,1].set_ylim([0,1])
    ax[0,1].set_xscale('log')
    ax[0,1].tick_params(axis="both", labelsize=12, direction = 'out', length = 5)
    ax[0,1].tick_params(axis="both", which = 'minor', direction = 'out', length = 3)
    ax[0,1].yaxis.set_ticks_position('both')
    ax[0,1].xaxis.set_ticks_position('both')
    ax[0,1].legend(loc = 'upper right', fontsize = 10)


    ax[1,0].set_title('$d$TEMP', fontsize = 14)
    ax[1,0].errorbar(bjd - 2400000, dTEMP, yerr = edTEMP, color = 'k', linestyle = "none", marker = "o", alpha = 0.75, capsize = 2, zorder = 1)
    ax[1,0].set_xlabel("Time (BJD - 2400000)", fontsize = 14)
    ax[1,0].set_ylabel(r"$d$TEMP (K)", fontsize = 14)
    ax[1,0].minorticks_on()
    ax[1,0].tick_params(axis="both", labelsize=12, direction = 'in', length = 5)
    ax[1,0].tick_params(axis="both", which = 'minor', direction = 'in', length = 3)
    ax[1,0].yaxis.set_ticks_position('both')
    ax[1,0].xaxis.set_ticks_position('both')

    ls = LombScargle(bjd, dTEMP, edTEMP)
    frequency, power = ls.autopower()
    fap = ls.false_alarm_level(0.001)
    ax[1,1].plot(1 / frequency, power, color="k", label = "Lomb-Scargle periodogram")
    ax[1,1].axhline(y=fap,linestyle="--",color="k",label="99.9% confidence")
    wf = LombScargle(bjd, np.ones(len(bjd)), fit_mean = False, center_data = False)
    frequency, power = wf.autopower()
    ax[1,1].plot(1/frequency, power, color="r", label = "Window function")
    ax[1,1].set_xlabel("Period (days)", fontsize=14)
    ax[1,1].set_ylabel("Power", fontsize=14)
    ax[1,1].set_xlim([1,1000])
    ax[1,1].set_ylim([0,1])
    ax[1,1].set_xscale('log')
    ax[1,1].tick_params(axis="both", labelsize=12, direction = 'out', length = 5)
    ax[1,1].tick_params(axis="both", which = 'minor', direction = 'out', length = 3)
    ax[1,1].yaxis.set_ticks_position('both')
    ax[1,1].xaxis.set_ticks_position('both')

    plt.tight_layout()
    plt.show()

def load_table(inputs):
    tbl_file = inputs['rdbfile']

    tbl = Table.read(tbl_file, format='ascii.rdb')
    if inputs['bin_per_epoch']:
        tbl = bin_table(tbl)

    if 'jd' not in tbl.colnames:
        tbl['jd'] = 2400000.5 + tbl['rjd']

    if 'sig_fwhm' in tbl.colnames:
        tbl['sfwhm'] = tbl['sig_fwhm']

    if 'valid_rv_range' in inputs['keplerian_fit'].keys():
        print('\tRemoving outliers outside of valid_rv_range = {} m/s'.format(inputs['keplerian_fit'][
                                                                                  'valid_rv_range']))
        dv = tbl['vrad'] - np.nanmedian(tbl['vrad'])
        g = (dv>inputs['keplerian_fit']['valid_rv_range'][0]) & (dv<inputs['keplerian_fit']['valid_rv_range'][1])
        tbl = tbl[g]

    if 'valid_jd_range' in inputs['keplerian_fit'].keys():
        print('\tRemoving outliers outside of valid_jd_range = {} '.format(inputs['keplerian_fit'][
                                                                                  'valid_jd_range']))
        jd = tbl['jd']
        g = (jd>inputs['keplerian_fit']['valid_jd_range'][0]) & (jd<inputs['keplerian_fit']['valid_jd_range'][1])
        tbl = tbl[g]


    # we remove the median of the velocities
    tbl['vrad'] -= np.nanmedian(tbl['vrad'])

    return tbl

default_yaml_folder = '/Users/eartigau/apero-utils/PRV-tools/TARGETS/'
yaml_file = 'proxima_2planets'

def kegp(yaml_file,cornerplot = False):

    # load yaml file with inputs
    if not os.path.exists(yaml_file):
        yaml_file = os.path.join(default_yaml_folder,yaml_file)
    if not yaml_file.endswith('.yaml'):
        yaml_file = '{}.yaml'.format(yaml_file)

    # we load yaml file
    print('\tLoading {}'.format(yaml_file))
    inputs = load_yaml(yaml_file)
    # we load the rdb file
    print('\tLoading {}'.format(inputs['rdbfile']))
    tbl = load_table(inputs)

    if  inputs['activity_correction']['type'] == 'Detrend':
        key_detrend = inputs['activity_correction']['quantity']

        print('\tDetrending against {}, derivative = {}'.format(key_detrend, inputs['activity_correction']['derivative']))
        quantity = tbl[key_detrend]
        err_quantity = tbl['s'+key_detrend]



        # we meed to update the amplitude of the GP, it may be completely wrong
        inputs['activity_correction']['gp_params']['amp'][0] = np.nanstd(quantity)
        if inputs['activity_correction']['gp_params']['amp'][1] == 'Gaussian':
            inputs['activity_correction']['gp_params']['amp'][2] = inputs['activity_correction']['gp_params']['amp'][0]/3.0
        if inputs['activity_correction']['gp_params']['amp'][1] == 'Jeffreys':
            inputs['activity_correction']['gp_params']['amp'][2] = [inputs['activity_correction']['gp_params']['amp'][0]/100,inputs['activity_correction']['gp_params']['amp'][0]*100]



        inputs2 = dict(inputs)
        if 'planet_params' in inputs2.keys():
            inputs2.pop('planet_params')

        if 'keplerian_fit' in inputs2.keys():
            inputs2['keplerian_fit']['jit'] = inputs2['activity_correction']['gp_params']['amp'][0]/10, inputs2['keplerian_fit']['jit'][1]

        for ite in range(2):
            mod, like, post = initialize_model(inputs2, tbl,key_quantity = key_detrend)
            res = optimize.minimize(
                post.neglogprob_array,  # objective function is negative log likelihood
                post.get_vary_params(),  # initial variable parameters
                method='Powell',  # Powell also works
            )


            ti = get_plot_time(tbl['jd'], timestep=inputs['keplerian_fit']['timestep'])
            gp_val, gp_var = like.predict(ti)
            gp_val2, gp_var2 = like.predict(tbl['jd'])
            offset = np.nanmedian(quantity-gp_val2)
            quantity -= offset
            print('offset : {}'.format(offset))

        if inputs['activity_correction']['derivative']:
            print('We decorrelate against the time derivative of {}'.format(key_detrend))
            nplots =3
        else:
            nplots = 2

        fig, ax = plt.subplots(nplots, 1, figsize=(12, 4*nplots),sharex=True)
        ax[0].errorbar(tbl['jd'], quantity, yerr=err_quantity, fmt='g.',label='data')
        ax[0].set(title=inputs['name'])
        ax[0].plot(ti, gp_val, 'r-', label='GP')
        ax[0].fill_between(ti,gp_val-(gp_var),gp_val+(gp_var),color='r',alpha=0.2)
        ax[0].set(ylabel = '{} [{}]'.format(inputs['activity_correction']['quantity'],unit(inputs['activity_correction']['quantity'])))
        ax[1].errorbar(tbl['jd'], quantity - gp_val2, yerr=np.sqrt(err_quantity**2+gp_var2), fmt='g.')
        ax[1].set(ylabel = 'residual [{}]'.format(unit(inputs['activity_correction']['quantity'])),xlabel='rjd')
        ax[0].legend()
        ax[0].grid(color='lightgrey', linestyle='--', linewidth=1)
        ax[1].grid(color='lightgrey', linestyle='--', linewidth=1)

        span = np.max(quantity+err_quantity) - np.min(quantity-err_quantity)
        ylim = [np.min(quantity)-0.05*span,np.max(quantity)+0.05*span]

        ax[0].set(xlim = [np.min(ti),np.max(ti)],ylim=ylim)

        if inputs['activity_correction']['derivative']:
            grad = np.gradient(gp_val,ti)
            ax[2].plot(ti,grad,'r-')
            ax[2].set(ylabel = 'GP derivative [{}/day]'.format(inputs['activity_correction']['quantity']),xlabel='jd')
            ax[2].grid(color='lightgrey', linestyle='--', linewidth=1)
        plt.tight_layout()
        plt.show()

        # add to the table the decorrelation quantity
        if inputs['activity_correction']['derivative']:
            decorrelation_spline =  ius(ti, grad, ext=3)
        else:
            decorrelation_spline =  ius(ti, gp_val, ext=3)

        print('for the sake of consistency, we cannot fit the GP on the RV data now')
        print('We set gpfit = "False"')
        inputs['keplerian_fit']['gpfit'] = False

        print('Best fit for detrending GP')
        for gp_key in post.params.keys():
            if gp_key.startswith('gp'):
                best_val = post.params[gp_key].value
                print('\t\t {} : {}'.format(gp_key,best_val))

    if inputs['activity_correction']['type'] == 'Train':
        key_detrend = inputs['activity_correction']['quantity']

        print('\tTraining on {}'.format(key_detrend))
        quantity = tbl[key_detrend]
        err_quantity = tbl['s' + key_detrend]
        quantity -= np.average(quantity,weights=1/err_quantity**2)


        #quantity -= np.nanmedian(quantity)

        #sig = sigma(quantity)
        #quantity /= sig
        #err_quantity /= sig

        inputs2 = dict(inputs)
        inputs2['activity_correction']['gp_params']['amp'][0] = sigma(quantity)

        if 'keplerian_fit' in inputs2.keys():
            inputs2['keplerian_fit']['jit'] = inputs2['activity_correction']['gp_params']['amp'][0]/10, inputs2['keplerian_fit']['jit'][1]

        # we meed to update the amplitude of the GP, it may be completely wrong
        if inputs2['activity_correction']['gp_params']['amp'][1] == 'Gaussian':
            inputs2['activity_correction']['gp_params']['amp'][2] = inputs2['activity_correction']['gp_params']['amp'][0]/3.0
        if inputs2['activity_correction']['gp_params']['amp'][1] == 'Jeffreys':
            inputs2['activity_correction']['gp_params']['amp'][2] = [inputs2['activity_correction']['gp_params'][
                                                                         'amp'][0]/1e3,inputs2['activity_correction'][
                'gp_params']['amp'][0]*1e3]

        if 'planet_params' in inputs2.keys():
            inputs2.pop('planet_params')

        for ite in range(2):
            mod, like, post = initialize_model(inputs2, tbl,key_quantity = key_detrend)
            res = optimize.minimize(
                post.neglogprob_array,  # objective function is negative log likelihood
                post.get_vary_params(),  # initial variable parameters
                method='Powell',  # Powell also works
            )


            ti = get_plot_time(tbl['jd'], timestep=inputs['keplerian_fit']['timestep'])
            gp_val, gp_var = like.predict(ti)
            gp_val2, gp_var2 = like.predict(tbl['jd'])
            offset = np.nanmedian(quantity-gp_val2)
            quantity -= offset
            print('offset : {}'.format(offset))


        gp_val_train = gp_val2
        gp_var_train = gp_var2

        fig, ax = plt.subplots(2, 1, figsize=(12, 8),sharex=True)
        ax[0].errorbar(tbl['jd'], quantity, yerr=err_quantity, fmt='g.',label='data')
        ax[0].set(title=inputs['name'])
        ax[0].plot(ti, gp_val, 'r-', label='Rotation GP')
        ax[0].fill_between(ti,gp_val-(gp_var),gp_val+(gp_var),color='r',alpha=0.2)
        ax[0].set(ylabel = '{} [{}]'.format(inputs['activity_correction']['quantity'],unit(inputs['activity_correction']['quantity'])))
        ax[1].errorbar(tbl['jd'], quantity - gp_val2, yerr=err_quantity, fmt='g.')
        ax[1].set(ylabel = 'residual [{}]'.format(unit(inputs['activity_correction']['quantity'])),xlabel='rjd')
        ax[0].legend()
        ax[0].grid(color='lightgrey', linestyle='--', linewidth=1)
        ax[1].grid(color='lightgrey', linestyle='--', linewidth=1)
        span = np.max(quantity +err_quantity) - np.min(quantity - err_quantity)
        ylim = [np.min(quantity - err_quantity)-0.05*span,np.max(quantity +err_quantity)+0.05*span]
        ax[0].set(xlim = [np.min(ti),np.max(ti)],ylim = ylim)

        plt.tight_layout()
        fig_file = inputs['rdbfile'].replace('.rdb','_GPfit.pdf')
        plt.savefig(fig_file)
        plt.show()

        if cornerplot:
            df = radvel.mcmc(post, nwalkers= inputs['mcmc_params']['nwalkers'], nrun=inputs['mcmc_params']['nrun'], savename=inputs['mcmc_params']['savename'])
            # chains = radvel.mcmc(post,nwalkers=nwalkers,
            #                     nrun=nrun,
            #                     ensembles=3,
            #                     savename=savename)
            print(post.bic())
            fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
            _ = corner.corner(df[post.name_vary_params()], labels=post.name_vary_params(), label_kwargs={"fontsize": 8},
                              plot_datapoints=False, bins=20, quantiles=[0.16, 0.5, 0.84],
                              show_titles=True, title_kwargs={"fontsize": 8}, smooth=True
                              )

            plt.savefig('{}_{}_corner.pdf'.format(inputs['name'],inputs['activity_correction']['quantity']), bbox_inches='tight')
            plt.show()

        for gp_key in post.params.keys():
            if gp_key.startswith('gp'):
                # one cannot change the amplitude of the GP
                if gp_key.split('gp_')[1] == 'amp':
                    continue

                new_val = post.params[gp_key].value
                print('\t\t{} : {}'.format(gp_key,new_val))
                inputs['activity_correction']['gp_params'][gp_key.split('gp_')[1]][0] = new_val
                # we fix the GP parameters as they were trained to the data on the other keyword
                inputs['activity_correction']['gp_params'][gp_key.split('gp_')[1]][-1] = False

        mod, like, post = initialize_model(inputs, tbl,key_quantity = 'vrad')
        res = optimize.minimize(
            post.neglogprob_array,  # objective function is negative log likelihood
            post.get_vary_params(),  # initial variable parameters
            method='Powell',  # Powell also works
        )

    if inputs['activity_correction']['type'] == 'BERV':
        key_x_value = 'BERV'

        quantity = tbl['vrad']
        err_quantity = tbl['svrad']

        inputs2 = dict(inputs)

        if 'planet_params' in inputs2.keys():
            inputs2.pop('planet_params')

        for ite in range(2):
            mod, like, post = initialize_model(inputs2, tbl,key_x_value = key_x_value)
            res = optimize.minimize(
                post.neglogprob_array,  # objective function is negative log likelihood
                post.get_vary_params(),  # initial variable parameters
                method='Powell',  # Powell also works
            )

            ti = np.arange(np.min(tbl['BERV'])-0.1,np.max(tbl['BERV'])+0.1,0.1)
            gp_val, gp_var = like.predict(ti)
            gp_val2, gp_var2 = like.predict(tbl['BERV'])
            offset = np.nanmedian(quantity-gp_val2)
            quantity -= offset
            print('offset : {}'.format(offset))

        fig, ax = plt.subplots(2, 1, figsize=(12, 8),sharex=True)
        ax[0].errorbar(tbl['BERV'], quantity, yerr=err_quantity, fmt='g.',label='data')
        ax[0].set(title=inputs['name'])
        ax[0].plot(ti, gp_val, 'r-', label='GP')
        #ax[0].fill_between(ti,gp_val-np.sqrt(gp_var),gp_val+np.sqrt(gp_var),color='r',alpha=0.2)
        ax[0].set(ylabel = '{} [{}]'.format(inputs['activity_correction']['quantity'],unit(inputs['activity_correction']['quantity'])))
        ax[1].errorbar(tbl['BERV'], quantity - gp_val2, yerr=np.sqrt(err_quantity**2+gp_var2), fmt='g.')
        ax[1].set(ylabel = 'residual [{}]'.format(unit(inputs['activity_correction']['quantity'])),xlabel='rjd')
        ax[0].legend()
        ax[0].grid(color='lightgrey', linestyle='--', linewidth=1)
        ax[1].grid(color='lightgrey', linestyle='--', linewidth=1)
        ax[0].set(xlim = [np.min(ti),np.max(ti)])

        plt.tight_layout()
        fig_file = inputs['rdbfile'].replace('.rdb','_RVfit.pdf')
        plt.savefig(fig_file)
        plt.show()

        spl_berv = ius(ti, gp_val)

    t, vel, errvel = tbl['jd'],  tbl['vrad'], tbl['svrad']
    # add a padding of 5% of the time span to the time grid

    ti = get_plot_time(t, timestep=inputs['keplerian_fit']['timestep'])
    ti_plot = get_plot_time(t, timestep=inputs['keplerian_fit']['timestep'])

    if inputs['activity_correction']['type'] == 'Detrend':
        decorrelation_fit = [0,0]
        nite = 3
    else:
        nite = 1

    # this table is decorrelated if we pass a decorrelation spline
    tbl_decorrelated = Table(tbl)

    if inputs['activity_correction']['type'] == 'BERV':
        tbl_decorrelated['vrad'] -= gp_val2

    for ite in range(nite):
        mod, like, post = initialize_model(inputs, tbl_decorrelated)

        print('BIC before first fitting : ',post.bic())
        print('\tFitting the model')
        res = optimize.minimize(
            post.neglogprob_array,     # objective function is negative log likelihood
            post.get_vary_params(),    # initial variable parameters
            method='Powell',           # Powell also works
            )

        if inputs['activity_correction']['type'] == 'Detrend':
            fit = np.polyfit(decorrelation_spline(tbl_decorrelated['jd']),post.residuals(),1,w=1/tbl['svrad']**2)
            tbl_decorrelated['vrad'] -= np.polyval(fit,decorrelation_spline(tbl_decorrelated['jd']))
            decorrelation_fit+=fit
            print(fit)



    print('BIC after first fitting : ',post.bic())

    zp = tbl['jd'][0]-tbl['plot_date'][0]

    print('rms : {:.3f} m/s, MAD : {:.3f} m/s'.format(np.nanstd(post.residuals()), np.nanpercentile(np.abs(
        post.residuals()), 68)))



    fig, ax = plt.subplots(nrows = 2, ncols = 2, sharex = False,sharey = False,figsize=(16,8),
                           gridspec_kw={'width_ratios':(8,1)})
    ax[0,0].plot_date(tbl['plot_date'],tbl['vrad'], 'k.',alpha = 0.01)
    ax[0,0].errorbar(tbl['plot_date'], tbl['vrad'], yerr=tbl['svrad'], fmt='o',color = 'r', ecolor='k', alpha=0.5)

    gp_val = None

    if (inputs['activity_correction']['type'] == 'RV') or (inputs['activity_correction']['type'] == 'Train'):
        ti_plot_gp = ti_plot - zp

        gp_val, gp_var = like.predict(ti_plot)
        gp_val2, gp_var2 = like.predict(tbl['jd'])
        offset = np.nanmedian(tbl['vrad']-gp_val2)
        gp_val += offset
        print('offset : {}'.format(offset))
        tbl_decorrelated['vrad'] -= gp_val2

        ax[0,0].plot(ti_plot_gp,gp_val,label = 'Rotation signal',alpha=0.9)
        #ax[0,0].fill_between(ti_plot_gp,gp_val-np.sqrt(gp_var),gp_val+np.sqrt(gp_var),color='r',alpha=0.1)

    if inputs['activity_correction']['type'] == 'BERV':
        gp_val = spl_berv(approx_berv(ti_plot-zp,tbl,timekey='plot_date'))
        ax[0,0].plot(ti_plot - zp,gp_val,label = 'BERV GP',alpha=0.9)


    kep_val = like.model(ti_plot)
    kep_val -= np.nanmean(kep_val)

    if type(gp_val) != type(None):
        ax[0,0].plot(ti_plot - zp, gp_val + kep_val, label='keplerian + rotation signal', alpha=0.5)
    else:
        if inputs['activity_correction']['type'] == 'Detrend':
            ax[0,0].plot(ti_plot - zp, np.polyval(decorrelation_fit,decorrelation_spline(ti_plot)), 'r-',
                                                label='decorrelation of {}'.format(inputs['activity_correction']['quantity']))

            ax[0,0].plot(ti_plot - zp, kep_val+np.polyval(decorrelation_fit,decorrelation_spline(ti_plot)), 'g-',
                                                label='total')

        else:
            ax[0,0].plot(ti_plot-zp,kep_val,label = 'keplerian',alpha=0.5)

    ax[1,0].plot_date(tbl['plot_date'],post.residuals(),'k.',alpha=0.1)
    ax[1,0].errorbar(tbl['plot_date'],post.residuals(),yerr = tbl['svrad'],fmt='o',color = 'r', ecolor='k', alpha=0.5)
    ax[1,0].set(xlabel = 'Date')
    ax[0,0].grid(color = 'black',alpha = 0.2)
    ax[1,0].grid(color = 'black',alpha = 0.2)
    ax[0,0].set(ylabel ='Velocity [m/s]')
    ax[1,0].set(ylabel ='Velocity residuals [m/s]')
    print(post)
    ax[0,0].legend()
    ylim = [1.1*np.min(tbl['vrad']-tbl['svrad']),1.1*np.max(tbl['vrad']+tbl['svrad'])]

    ax[0,0].set(xlim=(ti[0]-zp,ti[-1]-zp),ylim = ylim)
    ax[1,0].set(xlim=(ti[0]-zp,ti[-1]-zp),ylim = ylim)
    ax[0,1].set(ylim = ylim)
    ax[1,1].set(ylim = ylim)

    #xh, yh = plt.hist(post.residuals(), plot=False,xlim=ylim)
    title = '\t$\sigma$ : {:.2f} m/s'.format(np.nanstd(tbl['vrad']))
    ax[0,1].hist(tbl['vrad'], orientation='horizontal',bins=25,range=ylim,alpha = 0.6, color = 'orange')
    ax[0,1].text(0,ylim[0]+0.1*(ylim[1]-ylim[0]),title)
    #ax[0,1].set(title=title)
    title = '\t$\sigma$ : {:.2f} m/s'.format(np.nanstd(post.residuals()))
    ax[1,1].hist(post.residuals(), orientation='horizontal',bins=25,range=ylim,alpha = 0.6, color = 'orange')
    ax[1,1].text(0,ylim[0]+0.1*(ylim[1]-ylim[0]),title)
    ax[0,0].set(title = '{}'.format(inputs['name']))
    #ax[1,0].set(title = 'Residuals')
    ax[0,0].get_xaxis().set_ticks([])
    ax[0,1].get_xaxis().set_ticks([])
    ax[0,1].get_yaxis().set_ticks([])
    ax[1,1].get_xaxis().set_ticks([])
    ax[1,1].get_yaxis().set_ticks([])
    ax[0,1].grid(color = 'black',alpha = 0.2)
    ax[1,1].grid(color = 'black',alpha = 0.2)
    plt.tight_layout()
    fig_file = inputs['rdbfile'].replace('.rdb', '_RVfit_2panels.pdf')
    plt.savefig(fig_file)

    plt.show()


    tbl_output_file = inputs['rdbfile'].replace('.rdb','_output_GPs.rdb')
    tbl_output = Table()
    tbl_output['jd'] = tbl['jd']
    tbl_output['vrad'] = tbl['vrad']
    tbl_output['svrad'] = tbl['svrad']
    tbl_output['residuals'] = post.residuals()
    if inputs['activity_correction']['type'] != 'None':
        tbl_output['gp_val'] = gp_val2
        tbl_output['gp_err'] = np.sqrt(gp_var2)
        if inputs['activity_correction']['type'] == 'Train':
            tbl_output['gp_val_train_{}'.format(inputs['activity_correction']['quantity'])] = gp_val_train
            tbl_output['gp_err_train_{}'.format(inputs['activity_correction']['quantity'])] = np.sqrt(gp_var_train)
    tbl_output.write(tbl_output_file,format='ascii.rdb',overwrite=True)
    print('Output file: {}'.format(tbl_output_file))

    print('Sigma {:.2f} m/s, MAD {:.2f}'.format(np.nanstd(post.residuals()),np.nanpercentile(np.abs(post.residuals()),68)))
    #plt.hist(post.residuals())
    #plt.xlabel('Residuals [m/s]')
    #plt.ylabel('N')
    #plt.show()


    if (inputs['activity_correction']['type'] == 'RV') or (inputs['activity_correction']['type'] == 'Train'):
        RVPlot = orbit_plots.GPMultipanelPlot(
            post,
            subtract_gp_mean_model=False,
            plot_likelihoods_separately=False,
            subtract_orbit_model=False
        )
    else:
        RVPlot = orbit_plots.MultipanelPlot(post, legend=True)


    RVPlot.plot_multipanel()

    nwalkers = inputs['mcmc_params']['nwalkers']
    nrun = inputs['mcmc_params']['nrun']
    savename = inputs['mcmc_params']['savename']

    print(post.bic())
    plt.show()

    if 'DTEMP' in tbl_decorrelated.colnames:
        periodogram(tbl_decorrelated)

    if not cornerplot:
        return

    df = radvel.mcmc(post,nwalkers=nwalkers,nrun=nrun,savename=savename)
    #chains = radvel.mcmc(post,nwalkers=nwalkers,
    #                     nrun=nrun,
    #                     ensembles=3,
    #                     savename=savename)
    print(post.bic())
    fig, ax = plt.subplots(nrows = 1, ncols =1, figsize=(8,8))
    _ = corner.corner(df[post.name_vary_params()], labels=post.name_vary_params(), label_kwargs={"fontsize": 8},
                plot_datapoints=False, bins=20, quantiles=[0.16, 0.5, 0.84],
                show_titles=True, title_kwargs={"fontsize": 8}, smooth=True
            )

    plt.savefig('{}_corner.pdf'.format(inputs['name']), bbox_inches='tight')
    plt.show()

    quants = df.quantile([0.159, 0.5, 0.841]) # median & 1sigma limits of posterior distributions

    for par in post.params.keys():
        if post.params[par].vary:
            med = quants[par][0.5]
            high = quants[par][0.841] - med
            low = med - quants[par][0.159]
            err = np.mean([high,low])
            err = radvel.utils.round_sig(err)
            med, err, errhigh = radvel.utils.sigfig(med, err)
            print('{} : {} +/- {}'.format(par, med, err))
            post.model.params[par].value = med
        print('\t{}, {}'.format(par,post.params[par].value))

    print(post.bic())
