#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:34:45 2023

@author: caroline

Stellar model fitting utility functions (parallel runs)
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
import pdb
os.environ["OMP_NUM_THREADS"] = "1"
import stctm.exotune_parallel_setup as xps

import stctm.pystellspec as psspec
import stctm.exotune_utilities as xtu
import stctm.stellar_retrieval_utilities as sru
import astropy.constants as const
import emcee
import astropy.table as table
import h5py
import sys
import astropy.io as aio
import shutil

from multiprocessing import Pool
import scipy.signal as sig

#%%
def mcmc_run_serial(nwalkers, ndim, dtype, pos, nsteps, ncpu, other_args):
    print("Create sampler...")
    
    print("Running serial version of the MCMC fit!")
    sampler = emcee.EnsembleSampler(nwalkers, ndim, xtu.lnprob, 
                                    args=other_args,
                                    blobs_dtype=dtype)
    
    print("Running MCMC...")
    sampler.run_mcmc(pos, nsteps, progress=True, store=True)
    print("MCMC complete!!")
    
    return sampler

def mcmc_run_parallel(nwalkers, ndim, dtype, pos, nsteps, ncpu, other_args):
    spec_pickleformat, param, fitparanames, gaussparanames, mean_Gauss_para, std_Gauss_para,  fl_phot_spot_fac, Teffs_grid, loggs_grid,  wv_template_thisR, models_grid_fixedR,Fscale_guess = other_args 

    # for lnlike
    param0 = param
    fitparanames0 = fitparanames
    spec0 = spec_pickleformat
    models_baseline0 = fl_phot_spot_fac
    T_grid0 = Teffs_grid
    logg_grid0 = loggs_grid
    model_wv0 = wv_template_thisR
    models_grid0 = models_grid_fixedR

    # For get_param_priors
    gaussparanames0 = gaussparanames
    mean_Gauss_para0 = mean_Gauss_para
    std_Gauss_para0 = std_Gauss_para 
    
    def lnlike_parallel(theta):
        
        param =  param0
        fitparanames =  fitparanames0
        spec =  spec0
        models_baseline =  models_baseline0
        T_grid =  T_grid0 
        logg_grid =  logg_grid0
        model_wv =  model_wv0
        models_grid =  models_grid0
        
        waveMin, waveMax, yval, yerrLow = spec
    
        # expand models from the baseline models tuple
        model_phot_baseline, model_spot_baseline, model_fac_baseline = models_baseline 
        
        for i, paraname in enumerate(fitparanames):
            param[paraname] = theta[i]
        param = xtu.get_derived_param(param)
        
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
        
        if "Fscale" in fitparanames:
            Fscale = param["Fscale"]
        else:
            Fscale=1
            
    
            
        model_int = xtu.calc_spectrum_model_and_integrate(Fscale, param['fspot'],param['ffac'], spec, 
                                                      model_wv, model_phot_fixR, model_spot_fixR, 
                                                      model_fac_fixR)
    
    
        
        lnlk = -0.5*(np.sum((yval-model_int )**2./yerrLow**2.))
        
        return lnlk, model_int
    
    
    
    def get_param_priors_parallel():
        """
    
        Parameters
        ----------
        param : dict
            dictionary of default param values .
    
        Returns
        -------
        dictionary of parameter priors in the form of [low, high] for uniform priors.
    
        """
        T_grid =  T_grid0
        logg_grid =  logg_grid0
        param =  param0
        gaussparanames =  gaussparanames0
        mean_Gauss_para =  mean_Gauss_para0
        std_Gauss_para =  std_Gauss_para0
        
        parampriors = dict()
        parampriors['ffac'] = [0., 0.5]
        parampriors['fspot'] = [0., 0.5]
        parampriors['deltaTfac'] = [100, T_grid[-1]-param["Tphot"]]
        parampriors['deltaTspot'] = [T_grid[0]-param["Tphot"], -100.]
        # parampriors["Fscale"] = [0.8*Fscale_guess, 1.2*Fscale_guess]
        parampriors["Fscale"] = [0.01, 50]
        parampriors["loggphot"] = [np.max([logg_grid[0],2.5]), np.min([5.5, logg_grid[-1]])]
        parampriors["dlogghet"] = [logg_grid[0]-param["loggphot"], 0.]
        parampriors["Tphot"] = [T_grid[0], T_grid[-1]]
        
        for par in param:
            if par in gaussparanames:
                ind_gaussparam = np.where(gaussparanames == par)[0][0]
                new_prior = [np.max([mean_Gauss_para[ind_gaussparam] - 5.*std_Gauss_para[ind_gaussparam], parampriors[par][0]]), np.min([mean_Gauss_para[ind_gaussparam] + 5.*std_Gauss_para[ind_gaussparam], parampriors[par][1]])]
                # pdb.set_trace()
                parampriors[par] = new_prior
        # pdb.set_trace()
        return parampriors
    
    def lnprior_parallel(theta):
    
        param =  param0
        fitparanames =  fitparanames0
        gaussparanames =  gaussparanames0
        mean_Gauss_para =  mean_Gauss_para0
        std_Gauss_para =  std_Gauss_para0
    
        
        lp = 0.
        
        for i, paraname in enumerate(fitparanames):
            param[paraname] = theta[i]
        param = xtu.get_derived_param(param)
        
        parampriors = get_param_priors_parallel()
        for i, paraname in enumerate(fitparanames):
            # print("Prior on", paraname)
            # print("lp:", lp)
            # pdb.set_trace()
            if parampriors[paraname][0] < theta[i] <parampriors[paraname][1]:
                # print("Parampriors for this para:", parampriors[paraname])
                if paraname in gaussparanames:
                    ind_gaussparam = np.where(gaussparanames == paraname)[0][0]
                    # ind_gaussparam = np.where(list(gaussparanames) == paraname)
                    lp = lp - 0.5*((theta[i]-mean_Gauss_para[ind_gaussparam])**2./(std_Gauss_para[ind_gaussparam]**2.))
                else:
                    lp = lp + 0.0
            else:
                lp = - np.inf
        
        if param["fspot"] + param["ffac"]>1:
            lp = - np.inf            
    
        return lp
    
    def lnprob_parallel(theta):
        waveMin, waveMax, yval, yerrLow = spec0
        
        lp = lnprior_parallel(theta)
        
        if not np.isfinite(lp):
            return -np.inf, np.zeros_like(np.array(yval)) * np.nan
            # return -np.inf
        else:
            # pdb.set_trace()
    
            lnlk, model = lnlike_parallel(theta)
    
            lnprb = lp + lnlk
            return lnprb, model
    

    print("Create sampler...")
    

    print("Running parallel version of the MCMC fit on", ncpu, "CPUs!")
    pool = Pool(ncpu) # on 10 CPU

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_parallel,
                                    pool=pool,
                                    blobs_dtype=dtype)

    print("Running MCMC...")

    sampler.run_mcmc(pos, nsteps, progress=True, store=True)
    
    print("MCMC complete!")
    
    # Close the thread pool
    pool.close()
    pool.join()
    
    return sampler

