#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:34:45 2023

@author: caroline

Stellar model fitting utility functions
"""

#%%------------------- Import modules ------------------#
import os
# os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
# os.environ['CRDS_PATH'] = "/Users/caroline/crds_cache"
# os.environ['PYSYN_CDBS'] = "/Users/caroline/trds"

os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
os.environ['CRDS_PATH'] = "/home/caroline/crds_cache"
os.environ['PYSYN_CDBS'] = "/home/caroline/trds"
import numpy as np
import matplotlib.pyplot as plt
import os

os.environ["OMP_NUM_THREADS"] = "1"

import stctm.pystellspec as psspec
import stctm.exotune_utilities as xtu
import astropy.constants as const
import pdb
import emcee
import astropy.table as table
import h5py
from copy import deepcopy
import sys
import astropy.io as aio
import shutil

import stctm.stellar_retrieval_utilities as sru 
import matplotlib.gridspec as gridspec
import pandas as pd
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors

import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter

from astropy.convolution import convolve, Box1DKernel,  Gaussian1DKernel, Trapezoid1DKernel
import corner
from multiprocessing import Pool
# from emcee.utils import MPIPool

#%% Stats

def BIC(chi2,nDataPoints,nPara):
    '''
    Adapted from auxbenneke/utilities.py
    Liddle, A.R., 2007. Information criteria for astrophysical model selection. Monthly Notices of the Royal Astronomical Society: Letters 377, L74–L78. https://doi.org/10.1111/j.1745-3933.2007.00306.x
    The BIC assumes that the data points are independent and identically distributed, which may or may not be valid depending on the data set under consideration
    '''
    return chi2 + nPara*np.log(nDataPoints)

#%% Utility functions

def get_median_spectrum_and_unc(wave_um, spec_ts, flux_calibrated=False):

    # calculate median stellar spectrum
    print("Calculating median spectrum...")
    median_spectrum = np.nanmedian(spec_ts,axis=0)

    # calculate uncertainties on median stellar spectrum from ptp scatter in each wavelength band
    print("Calculating uncertainties in each wavelength bin from point to point scatter...")
    err_median_spectrum = np.nanstd(spec_ts,axis=0)
    if flux_calibrated:
        median_spectrum = median_spectrum * 1e-6 * 1e-23 * const.c.value * 1e+6  * 1e-4 / wave_um**2 * 1e4
        err_median_spectrum = err_median_spectrum * 1e-6 * 1e-23 * const.c.value * 1e+6  * 1e-4 / wave_um**2 * 1e4
    
    return median_spectrum, err_median_spectrum


def plot_median_spectrum_and_unc(wave_arr, median_spectrum, err_median_spectrum,label="Median from entire time series"):

    fig, ax = plt.subplots(1,1,figsize=(10,4))
    ax.plot(wave_arr, median_spectrum, color="k",zorder=0, label=label)
    ax.fill_between(wave_arr, median_spectrum-err_median_spectrum, median_spectrum+err_median_spectrum,color="gray", zorder=-9)
    
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel("Abs. calibrated flux [mJy]")
    ax.legend(loc=1)
    
    return fig, ax

def init_default_and_fitted_param(Tphot, met, logg_phot, fitspot=True, fitfac=True, 
                                  fitThet=False, fitTphot=False, fitFscale=True,
                                  fitlogg_phot=False,fitdlogg_het=False,
                                  fitLogfSpotFac=[0,0],
                                  fiterrInfl=False,
                                  Fscale_guess = 1, dlogg_het_guess = 0.0):
    """
    
    Parameters
    ----------
    Tphot : float
        photosphere temperature in K.
    met : float
        log stellar metallicity.
    logg_phot : float
        log surface gravity of the star (photosphere)
    fitspot : bool
        True to include spots in the model.
    fitfac : bool
        True to include faculae in the mode.
    fitThet : bool
        True to fit temperature(s) of heterogeneity(ies).
    fitTphot : bool
        True to fit stellar photosphere temperature
    fitFscale : bool
        True to fit for the scaling factor Fscale. Default is True.
    fitlogg_phot : bool
        True to fit for logg_phot. Default is False.
    fitdlogg_het : bool
        True to fit for dlogg_het. Default is False.
    Fscale_guess : float
        Initial guess for Fscale
    dlogg_het_guess : float
        Initial guess for dlogg_het. Default is 0.
    Returns
    -------
    dict of default param values and list of fitted parameters.

    """
    param = dict()
    fitparanames = []
    param["Tphot"] = Tphot
    param["met"] = met
    param['loggphot'] = logg_phot
    param["dlogghet"] =0.
    if dlogg_het_guess is None:
        param['logghet'] = logg_phot
    else:
        param['logghet'] = logg_phot + dlogg_het_guess
    param["deltaTspot"] = -150.
    param["deltaTfac"] = 150.
    param["Tspot"] = param["Tphot"] - 150.
    param["Tfac"] = param["Tphot"] + 150.
    param["Fscale"] = Fscale_guess
    param["logFscale"] = np.log10(Fscale_guess)
    param["logErrInfl"] = 0
    param["errInfl"] = 1    

    
    if fiterrInfl:
        param["logErrInfl"] = 1
        param["errInfl"] = 10
        fitparanames.append("logErrInfl")        
    if fitspot:       
        param["fspot"] = 0.02
        if fitLogfSpotFac[0]:
            param["log_fspot"] = np.log10(param["fspot"])
            fitparanames.append("log_fspot")
        else:
            fitparanames.append("fspot")
    else:
        param["fspot"] = 0.    
    
    if fitfac:
        param["ffac"] = 0.05
        if fitLogfSpotFac[1]:
            param["log_ffac"] = np.log10(param["ffac"])
            fitparanames.append("log_ffac") 
        else:
            fitparanames.append("ffac")
    else:
        param["ffac"] = 0.
    
    if fitTphot:
        fitparanames.append("Tphot")
        
    if fitThet:
        if fitspot:
            fitparanames.append("deltaTspot")
        if fitfac:
            fitparanames.append("deltaTfac")

    if fitFscale:
        fitparanames.append("logFscale")
    
    if fitlogg_phot:
        fitparanames.append("loggphot")
    if fitdlogg_het:
        fitparanames.append("dlogghet")
        
    return param, fitparanames

def get_derived_param(param):
    """

    Parameters
    ----------
    param : dict
        dictionary of default param values .

    Returns
    -------
    param with Tspot and Tfac updated
    """
    param["Tspot"] = param["Tphot"] + param["deltaTspot"]
    param["Tfac"] = param["Tphot"] + param["deltaTfac"]
    param["logghet"] = param["loggphot"] + param["dlogghet"]
    param["Fscale"] = 10**(param["logFscale"])
    param["errInfl"] = 10**(param["logErrInfl"])
    if "log_fspot" in param.keys():
        param["fspot"] = 10**param["log_fspot"] 
    if "log_ffac" in param.keys():
        param["ffac"] = 10**param["log_ffac"] 
    return param


# integrate model in wavelength range
def integ_model(lam_min,lam_max,wave,F):
    F_int = []
    dwave=wave[1:]-wave[:-1]

    for i in range(lam_min.size):
        ind = np.where((wave>=lam_min[i])*(wave<lam_max[i]))
        numerator = np.trapz(F[ind] * dwave[ind],x=wave[ind])
        denominator = np.trapz(dwave[ind],x=wave[ind])
        F_int.append(numerator/denominator)
#        F_int.append(np.trapz(F[ind],x=wave[ind]))
    return np.array(F_int)


def calc_spectrum_model_and_integrate(Fscale, fspot, ffac, spec, model_wv, model_phot_fixR, model_spot_fixR, model_fac_fixR):
    
    waveMin, waveMax, yval, yerrLow = spec
    # calculate stellar contamination model
    
    phot    = (1-ffac-fspot)*model_phot_fixR
    full_model = phot
    
    if fspot>0:
        spot    = fspot*model_spot_fixR
        full_model= full_model + spot
    if ffac>0:
        fac     = ffac*model_fac_fixR
        full_model = full_model + fac
    
    # integrate in bandpass
    model_int  = integ_model(waveMin,waveMax,model_wv, full_model)
    # pdb.set_trace()

    return Fscale*model_int

def Fscale_obs(Fscale,fspot,ffac,model_phot_fixR, model_spot_fixR, model_fac_fixR):
    phot    = (1-ffac-fspot)*model_phot_fixR
    full_model = phot
    
    if fspot>0:
        spot    = fspot*model_spot_fixR
        full_model= full_model + spot
    if ffac>0:
        fac     = ffac*model_fac_fixR
        full_model = full_model + fac
    return Fscale*full_model

def get_scaling_guess(param, T_grid, logg_grid,  model_wv, models_grid, spec,
                      overwrite_param=True, wave_min_match_um = 1.0, wave_max_match_um = 2.0):
    
    waveMin = np.array(spec.waveMin)
    # pdb.set_trace()
    waveMax = np.array(spec.waveMax)
    yval = np.array(spec.yval)
    yerrUpp = np.array(spec.yerrUpp)
    spec_pickleformat = (waveMin, waveMax, yval, yerrUpp ) 
    
    model_phot_fixR, _, _ = get_closest_models_from_grid(param, models_grid, T_grid, logg_grid)
    
    # identify wavelengths in common
    ind_obs = np.where((spec.wave>wave_min_match_um)*(spec.wave<wave_max_match_um))
    
    model_int = calc_spectrum_model_and_integrate(1,0,0, 
                                           spec_pickleformat, model_wv ,model_phot_fixR, model_phot_fixR, model_phot_fixR)


    median_model = np.nanmedian(model_int[ind_obs])
    median_obs = np.nanmedian(spec.yval[ind_obs])
    Fscale_guess = median_obs/median_model
                      
    print("\nGuessing Fscale=", Fscale_guess)
    if overwrite_param:
        param["logFscale"] = np.log10(Fscale_guess)
        param["Fscale"] = Fscale_guess
    return Fscale_guess, model_int

def get_closest_models_from_grid(param, models_grid_fixedR, Teffs_grid, loggs_grid):

    # Spot model
    ind_Tspot = np.argmin(np.abs(Teffs_grid-param["Tspot"]))
    ind_logg_het = np.argmin(np.abs(loggs_grid-param["logghet"]))
    model_spot_fixR  = models_grid_fixedR[ind_Tspot, ind_logg_het]
    
    # Faculae model
    ind_Tfac = np.argmin(np.abs(Teffs_grid-param["Tfac"]))
    ind_logg_het = np.argmin(np.abs(loggs_grid-param["logghet"]))
    model_fac_fixR  = models_grid_fixedR[ind_Tfac, ind_logg_het]
    
    # Photosphere model
    ind_Tphot = np.argmin(np.abs(Teffs_grid-param["Tphot"]))
    ind_logg_phot = np.argmin(np.abs(loggs_grid-param["loggphot"]))
    model_phot_fixR  = models_grid_fixedR[ind_Tphot,ind_logg_phot]
    
    return model_phot_fixR, model_spot_fixR, model_fac_fixR


def lnlike(theta, param, fitparanames, spec, models_baseline, T_grid, logg_grid,  model_wv, models_grid):
    """
    log-likelihood function

    Parameters
    ----------
    theta : array
        array of fitted parameter values.
    param : dict
        dictionary of values for all the model parameters.
    # param: Tphot, met, logg_phot, logg_spot, Tspot, fspot, Tfac, ffac
    fitparanames : list of str
        list of the fitted parameters.
    spec : pyStellSpec object ## TODO fix here 
        (stellar spectrum)
    models_baseline : list of arrays
        baseline models for the photosphere, spots, and faculae at fixed R.
    T_grid: numpy 1d array
        grid of temperatures for the grid of models (models_grid)
    logg_grid: numpy 1d array
        grid of loggs for the grid of models (models_grid)
    model_wv: numpy 1d array
        wavelength array of fixed R models
    models_grid: 3d array
        grid of models with varying temperatures and loggs, at fixed R

    Returns
    -------
    int
        ln likelihood value for parameters theta.

    """
    waveMin, waveMax, yval, yerrLow = spec

    # expand models from the baseline models tuple
    model_phot_baseline, model_spot_baseline, model_fac_baseline = models_baseline 
    
    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta[i]
    param = get_derived_param(param)
    
    if "deltaTspot" in fitparanames:
        ind_Tspot = np.argmin(np.abs(T_grid-param["Tspot"]))
        ind_logg_het = np.argmin(np.abs(logg_grid-param["logghet"]))
        model_spot_fixR  = models_grid[ind_Tspot, ind_logg_het]
    else:
        model_spot_fixR = model_spot_baseline
    
    if "deltaTfac" in fitparanames:
        ind_Tfac = np.argmin(np.abs(T_grid-param["Tfac"]))
        ind_logg_het = np.argmin(np.abs(logg_grid-param["logghet"]))
        model_fac_fixR  = models_grid[ind_Tfac, ind_logg_het]

    else:
        model_fac_fixR = model_fac_baseline

    if "Tphot" in fitparanames:
        ind_Tphot = np.argmin(np.abs(T_grid-param["Tphot"]))
        ind_logg_phot = np.argmin(np.abs(logg_grid-param["loggphot"]))
        model_phot_fixR  = models_grid[ind_Tphot,ind_logg_phot]

    else:
        model_phot_fixR = model_phot_baseline
    
    if "logFscale" in fitparanames:
        Fscale = param["Fscale"]
    else:
        Fscale=1
        

        
    model_int = calc_spectrum_model_and_integrate(Fscale, param['fspot'],param['ffac'], spec, 
                                                  model_wv, model_phot_fixR, model_spot_fixR, 
                                                  model_fac_fixR)


    yerrLow = yerrLow*param["errInfl"]
    lnlk = -0.5*(np.sum((yval-model_int )**2./yerrLow**2.))
    
    return lnlk, model_int
def get_param_priors(param, gaussparanames,hyperp_gausspriors=[], 
                     fitLogfSpotFac=[False,False], hyperp_logpriors=[],T_grid=[2300, 10000], 
                     logg_grid=[2.5,5.5], Fscale_guess=1.):
    """

    Parameters
    ----------
    param : dict
        dictionary of default param values .

    Returns
    -------
    dictionary of parameter priors in the form of [low, high] for uniform priors.

    """
    
    
    
    defaultpriors = dict()

    defaultpriors['ffac'] = [0., 0.9]
    defaultpriors['fspot'] = [0., 0.9]
    defaultpriors['deltaTfac'] = [100, T_grid[-1]-param["Tphot"]]
    defaultpriors['deltaTspot'] = [T_grid[0]-param["Tphot"], -100.]
    defaultpriors["Tphot"] = [T_grid[0], T_grid[-1]]
    defaultpriors["logFscale"] = [np.log10(Fscale_guess)-5, np.log10(Fscale_guess)+5]
    defaultpriors["logErrInfl"] = [0,np.log10(50)]

    defaultpriors["loggphot"] = [np.max([logg_grid[0],2.5]), np.min([5.5, logg_grid[-1]])]
    defaultpriors["dlogghet"] = [logg_grid[0]-param["logg_phot"], 0.]
    
    parampriors = dict()
    
    # check parameters for log priors
    allparanames = ['ffac',"fspot","deltaTfac","deltaTspot", "Tphot", "logFscale","loggphot","dlogghet","logErrInfl"]
    for par in allparanames:
        if par not in ["fspot", "ffac"]:
            parampriors[par] = defaultpriors[par]
        else:
            if par == "fspot":
                if fitLogfSpotFac[0]:
                    lowlim = hyperp_logpriors[0]
                    upplim = hyperp_logpriors[1]
                    parampriors["log_"+par] = [lowlim,upplim]
                else:
                    parampriors[par] = defaultpriors[par]
            elif par == "ffac":
                if fitLogfSpotFac[1]:
                    lowlim = hyperp_logpriors[0]
                    upplim = hyperp_logpriors[1]
                    parampriors["log_"+par] = [lowlim,upplim]
                else:
                    parampriors[par] = defaultpriors[par]
            else:
                parampriors[par] = defaultpriors[par]
                
                
                
                
    for par in param:
        if par in gaussparanames:
            ind_gaussparam = np.where(gaussparanames == par)[0][0]
            mean_gaussparam = hyperp_gausspriors[ind_gaussparam][0]
            std_gaussparam = hyperp_gausspriors[ind_gaussparam][1]
            new_prior = [np.max([mean_gaussparam - 5.*std_gaussparam, parampriors[par][0]]), np.min([mean_gaussparam + 5.*std_gaussparam, parampriors[par][1]])]
            # pdb.set_trace()
            parampriors[par] = new_prior
    # pdb.set_trace()
    return parampriors

def lnprior(theta, param, fitparanames, gaussparanames, hyperp_gausspriors, 
            fitLogfSpotFac, hyperp_logpriors, T_grid, logg_grid, Fscale_guess):
    """
    log-prior function

    Parameters
    ----------
    theta : array
        array of fitted parameter values.
    fitparanames : list of str
        list of the fitted parameters.
    gaussparanames : list of str
        list of the fitted parameters with Gaussian priors.
    mean_Gauss_para : numpy array
        array of mean parameters for all the params with Gaussian priors
    std_Gauss_para : numpy array
        array of std for all the params with Gaussian priors
    T_grid: numpy 1d array
        grid of temperatures for the grid of models (models_grid)
    logg_grid: numpy 1d array
        grid of loggs for the grid of models (models_grid)
    Fscale_guess : float
        Initial guess for Fscale
    Returns
    -------
    int
        ln prior value for parameters theta.

    """
    lp = 0.
    
    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta[i]
    param = get_derived_param(param)
    
    parampriors = get_param_priors(param, gaussparanames, hyperp_gausspriors=hyperp_gausspriors, 
                                   fitLogfSpotFac=fitLogfSpotFac, hyperp_logpriors=hyperp_logpriors,
                                   T_grid=T_grid, logg_grid=logg_grid, Fscale_guess=Fscale_guess)
    for i, paraname in enumerate(fitparanames):
        # print("Prior on", paraname)
        # print("lp:", lp)
        # pdb.set_trace()
        if parampriors[paraname][0] < theta[i] <parampriors[paraname][1]:
            # print("Parampriors for this para:", parampriors[paraname])
            if paraname in gaussparanames:
                ind_gaussparam = np.where(gaussparanames == paraname)[0][0]
                mean_gaussparam = hyperp_gausspriors[ind_gaussparam][0]
                std_gaussparam = hyperp_gausspriors[ind_gaussparam][1]
                # ind_gaussparam = np.where(list(gaussparanames) == paraname)
                lp = lp - 0.5*((theta[i]-mean_gaussparam)**2./(std_gaussparam**2.))

            else:
                lp = lp + 0.0
        else:
            lp = - np.inf
    
    if param["fspot"] + param["ffac"]>1:
        lp = - np.inf            

    return lp

def lnprob(theta, spec, param, fitparanames, gaussparanames, hyperp_gausspriors, 
           fitLogfSpotFac,hyperp_logpriors,models_baseline_fixR, T_grid, logg_grid, 
           model_wv,models_grid_fixR, Fscale_guess):
    
    lp = lnprior(theta, param, fitparanames, gaussparanames, hyperp_gausspriors, 
                 fitLogfSpotFac,hyperp_logpriors,T_grid, logg_grid, Fscale_guess)
    
    if not np.isfinite(lp):
        # return -np.inf, np.zeros_like(np.array(spec.yval)) * np.nan
        return -np.inf
    else:
        # pdb.set_trace()

        lnlk, model = lnlike(theta, param, fitparanames, spec, models_baseline_fixR, T_grid, 
                             logg_grid, model_wv, models_grid_fixR)

        lnprb = lp + lnlk
        return lnprb, model
        # return lnprb

    
#%%
## Plotting

def xspeclog(ax,xlim=None,level=1,fmt="%2.1f"):
    """ 
    adapted from auxbenneke/utilities.py
    """
    
    if level==0.5:
        majorticks = np.arange(0.2,6,0.2)
        minorticks = np.arange(0.2,6,0.1)
    if level==1:
        majorticks = [1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,8.0]
        minorticks = np.arange(0.3,6,0.1)
    if level==2:
        majorticks = np.r_[1.0,1.5,np.arange(2,9)]
        minorticks = None
    if level==3:
        majorticks = np.r_[1.0,1.5,2,3,4,6,8]
        minorticks = np.r_[5,7,9.0]
    if level==4:
        majorticks = np.r_[1,2,3,4,6,8]
        minorticks = np.r_[5,7,9]
    
    ax.set_xscale('log')
    ax.minorticks_on()

    if majorticks is not None:
        ax.xaxis.set_major_locator(ticker.LogLocator(subs=majorticks))
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter(fmt))
    if minorticks is not None:
        ax.xaxis.set_minor_locator(ticker.LogLocator(subs=minorticks))
        ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter(''))

    if xlim is not None:
        ax.set_xlim(xlim)
        
def setAxesFontSize(ax,fontsize):
    """ 
    adapted from auxbenneke/utilities.py
    """
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
        
def chainplot(samples,labels=None,nwalkers=None,fontsize=None,lw=1,stepRange=None):
    """ 
    adapted from auxbenneke/utilities.py
    """
    if samples.ndim==2:
        nsamp = samples.shape[0]     
        npara = samples.shape[1]
    elif samples.ndim==3:
        npara = samples.shape[2]
        
    if nwalkers is not None:     
        samples=samples.reshape([nwalkers,int(nsamp/nwalkers),npara])        
        
    nx=int(np.sqrt(npara))+1
    ny=int(npara*1.0/nx)+1
    
    fig, axes = plt.subplots(nx, ny, sharex=True, figsize=(20, 12))
    if axes.ndim==1:
        axes=np.array([axes])
        
    if stepRange is not None:
        steps = np.linspace(stepRange[0],stepRange[1],samples.shape[1])
    
    for i in range(npara):
        setAxesFontSize(axes[int(i/ny),i%ny],fontsize)
        if samples.ndim==2:
            axes[int(i/ny),i%ny].plot(samples[:, i].T, color="k", alpha=0.4, lw=lw,rasterized=False)
        elif samples.ndim==3:
            if stepRange is not None:
                axes[int(i/ny),i%ny].plot(steps,samples[:, :, i].T, color="k", alpha=0.4, lw=lw,rasterized=False)
            else:
                axes[int(i/ny),i%ny].plot(samples[:, :, i].T, color="k", alpha=0.4, lw=lw,rasterized=False)
        
        if labels is not None:
            axes[int(i/ny),i%ny].set_ylabel(labels[i])
        
    return fig,axes    

#%% Saving

def save_mcmc_to_pandas(results_folder, runname, sampler, burnin, ndim, fitparanames, save_fit):
    # write csv file and astropy table with samples outside of burn-in
    samples = pd.DataFrame(sampler.chain[:, burnin:, :].reshape((-1, ndim)),
                           columns=[x for x in fitparanames])
    lnprobability = pd.Series(sampler.lnprobability[:, burnin:].reshape(-1),
                              name='lnprobability')
    lnlikesave = pd.Series(sampler.lnprobability[:, burnin:].reshape(-1),
                       name='lnlike')
    panda = pd.concat([lnprobability, lnlikesave, samples], axis=1)
    t_res = table.Table.from_pandas(panda)
    if save_fit:
       aio.ascii.write(t_res, results_folder+"exotune_pandas_"+runname+'.csv', format='csv', overwrite=True)
       
    # save best fit parameters and quantiles
    bestfit, ind_bestfit, ind_maxprob = samp2bestfit(panda)
    print(bestfit)
    if save_fit:
        bestfit.to_csv(results_folder+'exotune_bestfit_'+runname+'.csv')
        
    # make corner plot
    print("\nMax Likelihood:\n")
    print(bestfit["MaxLike"])
    
    parabestfit = np.array(bestfit["MaxLike"][2:])
    return bestfit, ind_bestfit, ind_maxprob, parabestfit, samples, t_res

def samp2bestfit(samp):
    """
    Make pandas dataframe containing statistical info on the posterior
    """
    
    x50 = samp.quantile(0.50)
    x50.name = '50'
    x16 = samp.quantile(0.16)
    x16.name = '16'
    x84 = samp.quantile(0.84)
    x84.name = '84'
    x025 = samp.quantile(0.025)
    x025.name = '2.5'
    x975 = samp.quantile(0.975)
    x975.name = '97.5'

    ind = np.asarray(samp.lnprobability).argmax()
    MaxProb = samp.iloc[ind]
    MaxProb.name = 'MaxProb'

    ind_bestfit = np.asarray(samp.lnlike).argmax()
    MaxLike = samp.iloc[ind_bestfit]
    MaxLike.name = 'MaxLike'

    bestfit = pd.concat([x50, x16, x84, x025, x975, MaxLike, MaxProb], axis=1)
    bestfit['-1sigma'] = bestfit['50'] - bestfit['16']
    bestfit['+1sigma'] = bestfit['84'] - bestfit['50']
    bestfit['-2sigma'] = bestfit['50'] - bestfit['2.5']
    bestfit['+2sigma'] = bestfit['97.5'] - bestfit['50']

    return bestfit, ind_bestfit, ind

def save_bestfit_stats(spec, ind_bestfit, fitparanames, flat_oot_spec_models, results_folder, runname, save_fit=True):
    # create a dictionary that collates all the best-fit information
    print("Saving stats on the best fit...")
    bestfit_stats = dict()
    best_model = flat_oot_spec_models[ind_bestfit]
    nPara = len(fitparanames)
    nDataPoints = len(spec["yval"])
    n_dof = nDataPoints - nPara
    bestfit_stats["ind_postburnin"] = ind_bestfit
    bestfit_stats["chi2"] = np.sum((spec['yval']-best_model)**2./spec["yerrLow"]**2.)
    bestfit_stats["redchi2"] = bestfit_stats["chi2"]/n_dof
    bestfit_stats["BIC"] = BIC(bestfit_stats["chi2"] ,nDataPoints,nPara)
    t_bestfit_stats = table.Table([bestfit_stats])
    
    if save_fit:
        print("Writing to file...")
        aio.ascii.write(t_bestfit_stats, results_folder+"exotune_bestfit_stats_"+runname+'.csv', format='csv', overwrite=True)


def get_exotune_blobs(exotune_models_blobs, percentiles=[0.2, 2.3, 15.9, 84.1, 97.7, 99.8]):
 
    
    # for each epoch of each planet, get the median, 1 and 2sigma pred time
    exotune_percentiles = np.nanpercentile(exotune_models_blobs, percentiles, axis=0)
    return exotune_percentiles

def plot_exotune_samples_res(spec, param, fitparanames,
                          ind_bestfit, post_burnin_samples, T_grid, logg_grid,
                          modelgrid_wave, 
                          grid_models_fixedresP,
                          sample_spectra=None, modelgrid_resP=10000,
                          target_resP=100,N_samp=1000, ax=None,
                          bestfit_color = 'k', color="b",plot3sig=True,
                          plot2sig=True, plot1sig=True, plotmedian=True,
                          plotbestfit=True, legend_loc=1, save_csv=True, 
                          results_folder="", runname=""):
    

    
    # get the target wavelength array
    kernel=Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)    # *2 because FWHM = 2 standard deviation
    l = int(kernel.shape[0]/2)
    wv_array = modelgrid_wave[l:-l]
    # ind = np.where((wv_array>spec["waveMin"][0])*(wv_array<spec["waveMax"][-1]))
    
    if sample_spectra is None: # unless the grid of sample spectra is provided
        # draw random samples from the post-burnin samples
        if len(post_burnin_samples)<N_samp:
            N_samp = len(post_burnin_samples)
            print("Too few samples, using N_samp=", N_samp)
        
        ind_samples = np.random.choice(np.arange(len(post_burnin_samples)), size=N_samp, replace=False)
        
        
        
        # calculate array of sample spectra
        sample_spectra = np.zeros((N_samp, wv_array.size))
       
        for i in range(N_samp):
            if i%100 == 0:
                print("i=", i+1, "/", N_samp)
            theta = post_burnin_samples.iloc[ind_samples[i]]
            for j, paraname in enumerate(fitparanames):
                param[paraname] = theta[j]
            param = get_derived_param(param)
            
            ind_Tspot = np.argmin(np.abs(T_grid-param["Tspot"]))
            ind_Tfac = np.argmin(np.abs(T_grid-param["Tfac"]))
            ind_Tphot = np.argmin(np.abs(T_grid-param["Tphot"]))
            ind_logg_phot = np.argmin(np.abs(logg_grid-param["loggphot"]))
            ind_logg_het = np.argmin(np.abs(logg_grid-param["logghet"]))
            
            model_spot = grid_models_fixedresP[ind_Tspot, ind_logg_het]
            model_phot = grid_models_fixedresP[ind_Tphot, ind_logg_phot]
            model_fac = grid_models_fixedresP[ind_Tfac, ind_logg_het]
            
            model_oot = Fscale_obs(param["Fscale"],param['fspot'],param['ffac'],
                                   model_phot,model_spot,model_fac)
                        
            # convolve the model to the target resolving power
            kernel=Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)    # *2 because FWHM = 2 standard deviation
            l = int(kernel.shape[0]/2)
            # x = wave[l:-l]
            model_oot = convolve(model_oot, kernel)[l:-l]
            # pdb.set_trace()
            sample_spectra[i] = deepcopy(model_oot)
    
    # calculate spectrum for best-fit parameters
    theta = post_burnin_samples.iloc[ind_bestfit]
    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta[i]
    param = get_derived_param(param)
    
    ind_Tspot = np.argmin(np.abs(T_grid-param["Tspot"]))
    ind_Tfac = np.argmin(np.abs(T_grid-param["Tfac"]))
    ind_Tphot = np.argmin(np.abs(T_grid-param["Tphot"]))
    ind_logg_phot = np.argmin(np.abs(logg_grid-param["loggphot"]))
    ind_logg_het = np.argmin(np.abs(logg_grid-param["logghet"]))
    
    model_spot = grid_models_fixedresP[ind_Tspot, ind_logg_het]
    model_phot = grid_models_fixedresP[ind_Tphot, ind_logg_phot]
    model_fac = grid_models_fixedresP[ind_Tfac, ind_logg_het]
    
    model_oot = Fscale_obs(param["Fscale"],param['fspot'],param['ffac'],
                           model_phot,model_spot,model_fac)
    
    # convolve the model to the target resolving power
    kernel=Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)    # *2 because FWHM = 2 standard deviation
    l = int(kernel.shape[0]/2)
    # x = wave[l:-l]
    oot_spec_best = convolve(model_oot, kernel)[l:-l]
    
    # calculate percentiles from sample spectra
    perc_spectra = get_exotune_blobs(sample_spectra, percentiles=[0.2, 2.3, 15.9, 50., 84.1, 97.7, 99.8])
    
    lowest_z = 1000
    
    
    if ax is None:
        fig, ax = spec.plot(marker=".",markersize=1,ls="",color="gray",zorder=lowest_z+3, label="Observed spectrum")

    else:
        fig = None
        spec.plot(ax=ax,marker=".",markersize=1,color="gray",zorder=lowest_z+3, label="Observed spectrum")



    if plot1sig:
        ax.fill_between(wv_array,perc_spectra[2],perc_spectra[4],color=color, alpha=0.5,zorder=lowest_z+1,label=r'1$\sigma$')


    if plot2sig:
        ax.fill_between(wv_array,perc_spectra[1],perc_spectra[5],color=color, alpha=0.3,zorder=lowest_z,label=r'2$\sigma$')

    if plot3sig:
        ax.fill_between(wv_array,perc_spectra[0],perc_spectra[6],color=color, alpha=0.2,zorder=lowest_z,label=r'3$\sigma$')

    if plotmedian:
        ax.plot(wv_array,perc_spectra[3],color=color,zorder=lowest_z+2,label=r'Median')
    
    if plotbestfit: 
        ax.plot(wv_array,oot_spec_best,color=bestfit_color,zorder=lowest_z+2,label=r'Max. likelihood')

    
    ax.legend(loc=legend_loc)
    
    if save_csv:
        dct = dict()
        dct["wave"] = deepcopy(wv_array)
        dct["bestfit"] = oot_spec_best
        dct["median"] = perc_spectra[3]
        dct["+1 sigma"] = perc_spectra[4]
        dct["-1 sigma"] = perc_spectra[2]
        dct["+2 sigma"] = perc_spectra[5]
        dct["-2 sigma"] = perc_spectra[1]    
        dct["+3 sigma"] = perc_spectra[6]
        dct["-3 sigma"] = perc_spectra[0] 
        
        t = table.Table(data=dct)
        t.write(results_folder+"fixedR_1_2_3_sigma"+runname+".csv", overwrite=True)


    return fig, ax, sample_spectra

def get_labels_from_fitparanames(fitparanames):
    labels = []
    for p in fitparanames:
        if p == "fspot":
            labels.append(r"$f_\mathrm{spot}$")
        elif p == "ffac":
            labels.append(r"$f_\mathrm{fac}$")  
        elif p == "log_fspot":
            labels.append(r"$\log_{10} f_\mathrm{spot}$")
        elif p == "log_ffac":
            labels.append(r"$\log_{10} f_\mathrm{fac}$") 
        elif p == "deltaTfac":
            labels.append(r"$\Delta T_\mathrm{fac} [K]$")           
        elif p == "deltaTspot":
            labels.append(r"$\Delta T_\mathrm{spot} [K]$")    
        elif p == "Tphot":
            labels.append(r"$T_\mathrm{phot}$ [K]") 
        elif p == "logErrInfl":
            labels.append(r"log$_{10}$ err. fac.") 
        elif p == "logFscale":
            labels.append(r"log$_{10}$ $F$") 
        elif p == "logghet":
            labels.append(r"log $g_\mathrm{het}$") 
        elif p == "dlogghet":
            labels.append(r"$\Delta$ log $g_\mathrm{het}$") 
        elif p == "loggphot":
            labels.append(r"log $g_\mathrm{phot}$") 
    return labels

def plot_corner(samples, plotparams, plot_datapoints=False, smooth=1.5,
                        quantiles=[0.16, 0.5, 0.84], title_kwargs={'fontsize':14},
                        hist_kwargs={"linewidth":3}, rg=None, 
                        show_titles=True, levels=(0.393,0.865,0.989), **kwargs):
    """
    Corner plot for an emcee fit of the water mass fraction that matches
    the observed planet params
    
    samples: generated by emcee sampler
    plotparams: plot params
    other args: args for the corner function
    
    Returns the figure with the corner plot 
    """
    hist_kwargs["color"] = plotparams["hist_color"]
    color = plotparams["hist_color"]
    fig = corner.corner(samples, labels=plotparams["labels"], 
                        plot_datapoints=plot_datapoints, smooth=smooth,
                        show_titles=show_titles, quantiles=quantiles,
                        title_kwargs=title_kwargs, color=color,
                        hist_kwargs=hist_kwargs, range=rg, 
                        levels=levels, **kwargs)
    return fig

def plot_custom_corner(samples, fitparanames, parabestfit, param):
    # reorder samples
    # if "logErrInfl" in fitparanames:
        # ordered_samples = np.zeros((samples.shape[0],samples.shape[1]-1))
    # else:
    ordered_samples = np.zeros((samples.shape[0],samples.shape[1]))
    # ordered_fitparanames_all = ["fspot", "deltaTspot", "ffac", "deltaTfac", "logg_het", "Tphot", "logFscale"]
    ordered_fitparanames_all = ["fspot", "deltaTspot", "ffac", "deltaTfac", "dlogghet", "Tphot", "loggphot", "logFscale","logErrInfl"]
    # ordered_fitparanames_all = ["fspot", "deltaTspot", "ffac", "deltaTfac", "dlogghet", "Tphot", "loggphot", "logFscale"]
    ordered_fitparanames = []
    ind = np.where(np.array(fitparanames)=="logFscale")
    Fscale_guess =  10**parabestfit[ind[0][0]]
    # for nice plotting
    # param = dict()
    # ind = np.where(np.array(fitparanames)=="Tphot")
    # param["Tphot"] = parabestfit[ind[0][0]]
    # ind = np.where(np.array(fitparanames)=="loggphot")
    # param["loggphot"] = parabestfit[ind[0][0]]
    parampriors = get_param_priors(param,[],Fscale_guess=Fscale_guess)
    
    ordered_rgs = [] # ranges for each parameter
    for p in ordered_fitparanames_all:
        if p in fitparanames:
            ordered_fitparanames.append(p)
            
    ordered_parabestfit = np.zeros_like(parabestfit)
    for i, p in enumerate(ordered_fitparanames):
        ind = np.where(np.array(fitparanames)==p)
        # print(ind)
        ordered_samples[:,i] = np.array(samples)[:, ind[0][0]]
        ordered_parabestfit[i] = parabestfit[ind[0][0]]
        # std_samples = np.nanstd(ordered_samples[:,i])
        # med_samples = np.nanmedian(ordered_samples[:,i])
        min_samples = np.percentile(ordered_samples[:,i], 0.05)
        max_samples = np.percentile(ordered_samples[:,i], 99.95)
        low_range = np.max([min_samples, parampriors[p][0]])
        upp_range = np.min([max_samples, parampriors[p][1]])
        ordered_rgs.append([low_range,upp_range])
    
    print("ordered_fitparanames: ", ordered_fitparanames)
    print("ranges:", ordered_rgs)
    # ind = np.where(get_labels_from_fitparanames(ordered_fitparanames)!="logErrInfl")[0]

    # pdb.set_trace()
    plotparams = dict()
    plotparams["hist_color"] = "C0"
    plotparams["labels"] = get_labels_from_fitparanames(ordered_fitparanames)

    fig = plot_corner(ordered_samples, plotparams, smooth=0.8, fill_contours=True, 
                      truths = ordered_parabestfit, truth_color="k")
    # fig = plot_corner(ordered_samples[:,ind], plotparams, smooth=0.8, fill_contours=True, 
    #                   truths = ordered_parabestfit[ind], truth_color="k")    
                                    

    fig.set_dpi(50)  
    fig.set_figheight(11)
    fig.set_figwidth(15)   
    
    return fig

def format_param_str(param, fitparanames):
    param = get_derived_param(param)

    str1 = "Fitted params: "+str(fitparanames)+"\n"
    str1 = str1 + "Stellar params: Tphot="+str(int(param["Tphot"]*10.)/10.)+" met="+str(param["met"]) + "\n"
    str1 = str1 + "logg_phot="+str(param["loggphot"]) + " logg_het="+str(param["logghet"])+ "\n"
    str1 = str1 + "Tspot="+str(int(param["Tspot"]*10.)/10.)+" Tfac="+str(int(param["Tfac"]*10.)/10.)+"\n"
    str1 = str1 + "fspot=" + str(int(param["fspot"]*10000.)/10000.) + " ffac="+str(int(param["ffac"]*10000)/10000.)
    return str1

def plot_maxlike_and_maxprob(spec, param, parabestfit, ind_maxprob, ind_bestfit, fitparanames, flat_st_ctm_models, pad=0.25):
    # Plot best fit stellar contamination model
    fig, ax = spec.plot()
    param_bestfit = deepcopy(param)
    for i, p in enumerate(fitparanames):
        param_bestfit[p] = parabestfit[i]
    ax.scatter(spec.wave, flat_st_ctm_models[ind_maxprob], label="Max. Probability", color="slateblue", alpha=1.)
    ax.scatter(spec.wave, flat_st_ctm_models[ind_bestfit], label="Max. Likelihood", color="r", alpha=0.5, marker=".", s=50)

    ax.text(np.min(spec.waveMin)-pad/4, 1.08*np.median(spec.yval), format_param_str(param_bestfit, fitparanames), fontsize=10, fontweight="bold", color="k")

    ax.set_xlim(np.min(spec.waveMin)-pad/2, np.max(spec.waveMax)+pad)
    # ax.set_ylim(0.8*np.median(spec['yval']), 1.15*np.median(spec['yval']))
    
    ax.set_title("Best-fit model")
    ax.legend(loc=4)
    return fig, ax  