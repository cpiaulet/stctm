#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:34:45 2023

@author: caroline

Stellar retrieval utility functions
"""

import os
os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
os.environ['CRDS_PATH'] = "/Users/caroline/crds_cache"
os.environ['PYSYN_CDBS'] = "/Users/caroline/trds"
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
from astropy.io import fits
from scipy.optimize import minimize
import pdb
import pkg_resources
import emcee
import corner
import pysynphot as psp
import astropy.table as table
import h5py
import random as rdm
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors

import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
from copy import deepcopy
import sys
import pandas as pd
import astropy.io as aio
import time
import matplotlib.gridspec as gridspec
from astropy.convolution import convolve, Box1DKernel,  Gaussian1DKernel, Trapezoid1DKernel

#%% --- Constants --- #

c = const.c.value
h = const.h.value

#%% --- Helper functions --- #

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

# model for observed transit depth (Rackham+2018)

def Dlam_obs(Dlam,fspot,ffac,Fphot_int,Fspot_int,Ffac_int):
    spot    = fspot*(1.-(Fspot_int/Fphot_int))
    fac     = ffac*(1.-(Ffac_int/Fphot_int))
    epsilon = 1./(1.-spot-fac)
    return Dlam*epsilon

def calc_stctm_model_and_integrate(Dlam, fspot, ffac, spec, model_wv, model_phot_fixR, model_spot_fixR, model_fac_fixR):
    # calculate stellar contamination model
    spot    = fspot*(1.-(model_spot_fixR/model_phot_fixR))
    fac     = ffac*(1.-(model_fac_fixR/model_phot_fixR))
    epsilon = 1./(1.-spot-fac)
    
    # integrate in bandpass
    # pdb.set_trace()

    epsilon_int  = integ_model(spec['waveMin'],spec['waveMax'],model_wv, epsilon)

    return Dlam*epsilon_int, Dlam*epsilon
def make_wavelength_array(wv_min_um=0.2, wv_max_um=6.0, resPower=1000, use_pymsg=True):
    """
    

    Parameters
    ----------
    wv_min_um : float, optional
        min wavelength of the array, in um. The default is 0.2.
    wv_max_um : float, optional
        max wavelength of the array, in um. The default is 6.0.
    resPower : float, optional
        target resPower. The default is 1000.
    use_pymsg : bool, optional
        whether or not the wave edges should be returned for pymsg

    Returns
    -------
    wavelength array at the specified resolving power

    """
    wv_list = [wv_min_um]
    while wv_list[-1]<wv_max_um:
        dlam = wv_list[-1]/resPower
        wv_list.append(wv_list[-1]+dlam)
    wv_array = np.array(wv_list)
    if use_pymsg == False:
        return wv_array
    else:
        wv_edges = 0.5 * (wv_array[1:] + wv_array[:-1])
        return wv_array[1:-1], wv_edges

def load_phoenix_model(Teff, M, logG, wv_target=None, wv_edges=None, resPower_target=None, 
                       wv_min_um=0.2, wv_max_um=6.0, use_pymsg = True, pymsg_specgrid=None):
    """
    

    Parameters
    ----------
    Teff : float
        Star effective temperature.
    M : float
        Star metallicity.
    logG : float
        Star log g.
    wv_target : array, optional
        Target wavelength array. The default is None.
    resPower_target : float, optional
        Target wavelength array resolving power
    wv_min_um : float, optional
        min wavelength of the array, in um. The default is 0.2.
    wv_max_um : float, optional
        max wavelength of the array, in um. The default is 6.0.
    use_pymsg : bool, optional
        whether or not to use pymsg to load the model
    Returns
    -------
    wavePhoenix : array
        wavelength array.
    fStarSurfPhoenix : array
        flux array.

    """
    if use_pymsg:
        x = {'Teff': Teff, 'log(g)': logG, "[Fe/H]": M}
        
        if wv_edges is None:
            if resPower_target is None:
                print("Error, resPower_target and wv_edges are both None")
                return np.nan
            else:
                wv_target, wv_edges = make_wavelength_array(wv_min_um=wv_min_um, wv_max_um=wv_max_um, resPower=resPower_target, use_pymsg=use_pymsg)
                f_pymsg = pymsg_specgrid.flux(x, wv_edges*1e4)
                fStarSurfPymsg = convertIntensity(f_pymsg, wv_target, InputUnit='(ergs/s)/(cm**2*A)', WavelengthUnit='um', OutputUnit='W/(m**2*um)')
                return wv_target, fStarSurfPymsg
        else:    
            f_pymsg = pymsg_specgrid.flux(x, wv_edges*1e4)
                    
            fStarSurfPymsg = convertIntensity(f_pymsg, wv_target, InputUnit='(ergs/s)/(cm**2*A)', WavelengthUnit='um', OutputUnit='W/(m**2*um)')
        
        return wv_target, fStarSurfPymsg
    else:
        #load Phoenix model
        sp = psp.Icat('phoenix', Teff, M, logG)
        # get wave and flux in expected units
        wavePhoenix = convertWave(sp.wave, 'A', 'um')
        fStarSurfPhoenix = convertIntensity(sp.flux, wavePhoenix, InputUnit='(ergs/s)/(cm**2*A)', WavelengthUnit='um', OutputUnit='W/(m**2*um)')
        
        if wv_target is None:
            if resPower_target is None:
                return wavePhoenix, fStarSurfPhoenix
            else:
                wv_target = make_wavelength_array(wv_min_um=wv_min_um, wv_max_um=wv_max_um, resPower=resPower_target)
                fStarSurf = np.zeros_like(wv_target)
                fStarSurf = np.interp(wv_target, wavePhoenix, fStarSurfPhoenix)
                return wv_target, fStarSurf
    
        else:
            
            fStarSurf = np.zeros_like(wv_target)
            fStarSurf = np.interp(wv_target, wavePhoenix, fStarSurfPhoenix)
            return wv_target, fStarSurf


def convertWaveToSi(inp,unit1):
    """ 
    adapted from auxbenneke/radutils.py
    """

    # WavelengthUnit
    if unit1=='Hz':     #Frequency
        f=inp
        wave=c/f
        wavenumber=1/wave
    elif unit1=='m':    #Wavelength
        wave=inp
        f=c/wave
        wavenumber=1/wave
    elif unit1=='um':
        wave=inp*1e-6
        f=c/wave
        wavenumber=1/wave
    elif unit1=='nm':
        wave=inp*1e-9
        f=c/wave
        wavenumber=1/wave
    elif unit1=='A':
        wave=inp*1e-10
        f=c/wave
        wavenumber=1/wave
    elif unit1=='m**-1': #Wavenumber
        wavenumber=inp
        f=wavenumber*c
        wave=1/wavenumber
    elif unit1=='cm**-1':
        wavenumber=inp*1e2
        f=wavenumber*c
        wave=1/wavenumber
    else:
        print('Input Unit Error!!!')

    ##Now wave is in meters
    ##Now f is in 1/s
    ##Now wavenumber is in meters**-1

    return wave, f, wavenumber


def convertWave(inp,unit1,unit2):
    """ 
    adapted from auxbenneke/radutils.py
    """
    wave, f, wavenumber = convertWaveToSi(inp,unit1)

    ##Now wave is in meters
    ##Now f is in 1/s
    ##Now wavenumber is in meters**-1

    # WavelengthUnit
    if unit2=='Hz':     #Frequency
        output = f
    elif unit2=='m':    #Wavelength
        output = wave
    elif unit2=='um':
        output=wave*1e6
    elif unit2=='nm':
        output=wave*1e9
    elif unit2=='A':
        output=wave*1e10
    elif unit2=='m**-1': #Wavenumber
        output=wavenumber
    elif unit2=='cm**-1':
        output=wavenumber*1e-2
    else:
        print('Output Unit Error!!!')

    return output

def convertIntensity(Iin,LambdaInput,InputUnit='W/(m**2*Hz)',WavelengthUnit='um',OutputUnit='W/(m**2*um)'):
    """ 
    adapted from auxbenneke/radutils.py
    """
    ##Unit converter for flux and intensity
    ##Ex: Iout = ConvertIntensityUnits(Iin,'W/(m**2*Hz)','W/(m**2*um)',wave,'m'
    ##    plot(wave*1e6,ConvertIntensityUnits(IncidentFlux,'W/(m**2*Hz)','W/(m**2*um)',wave,'m'
    ##    trapz(wave*1e6,ConvertIntensityUnits(IncidentFlux,'W/(m**2*Hz)','W/(m**2*um)',wave,'m'


    wave, f, wavenumber = convertWaveToSi(LambdaInput,WavelengthUnit)


    ##Now wave is in meters
    ##Now f is in 1/s
    ##Now wavenumber is in meters**-1


    #Convert Iin into 'W/(m**2*Hz)'
    if InputUnit=='W/(m**2*Hz)':                                  #Frequency Bin
        I=Iin
    elif InputUnit=='Jy':                                  #Frequency Bin
        I=Iin  /  (1e26)
    elif InputUnit=='W/(m**2*m)':                               #Wavelength Bin
        I=Iin  /  (c/(wave**2))
    elif InputUnit=='W/(m**2*um)':
        I=Iin  /  (c/(wave**2)  *1e-6)    #checked!!!
    elif InputUnit=='W/(m**2*nm)':
        I=Iin  /  (c/(wave**2)  *1e-9)    #checked!!!
    elif InputUnit=='W/(m**2*A)':
        I=Iin  /  (c/(wave**2)  *1e-10)    #checked!!!

    elif InputUnit=='(ergs/s)/(m**2*A)':
        I=Iin  /  (c/(wave**2)  *1e-10  * 1e7)

    elif InputUnit=='(ergs/s)/(cm**2*A)':
        I=Iin  /  (c/(wave**2) *1e-10  * 1e7  *1e-4)

    elif InputUnit=='(ergs/s)/(cm**2*cm)':
        I=Iin  /  (c/(wave**2) *1e-10  * 1e7  *1e-4  *10*1e3*1e3*10  )

    elif InputUnit=='W/(m**2*m**-1)':   #= W*m**-1 = W/m          #Wavenumber Bin
        I=Iin  /  (c)
    elif InputUnit=='W/(cm**2*cm**-1)':   #= W*cm**-1 = W/cm
        I=Iin  /  (c  *  1e-2)   #checked
    elif InputUnit=='W/(m**2*cm**-1)':
        I=Iin  /  (c  *  1e-2  *  1e4)
    elif InputUnit=='(photons/s)/(m**2*Hz)':     #-----Photon count------
        I=Iin  /  (1/(h*f))
    elif InputUnit=='(photons/s)/(m**2*um)':     #-----Photon count------
        I=Iin  /  (c/(wave**2)  *1e-6  /(h*f))
    elif InputUnit=='(photons/s)/(m**2*cm**-1)':     #-----Photon count------
        I=Iin  /  (c  *  1e-2  *  1e4    /(h*f))
    elif InputUnit=='(photons/s)/(cm**2*A)':     #-----Photon count------  #checked
        I=Iin  /  ((c/(wave**2)  *1e-6  /(h*f)) *1e-4 *1e-4)
    elif InputUnit=='(photons/s)/(m**2*cm**-1)':     #-----Photon count------ #checked
        I=Iin  /  (c  *  1e-2  *  1e4    /(h*f))
    else:
        print('Input Unit Error!!!')

    # Now:   I is in 'W/(m**2*Hz)'


    #Convert I into selected Output Unit
    if OutputUnit=='W/(m**2*Hz)':                                  #Frequency Bin
        Iout=I
    elif OutputUnit=='Jy':                                  #Frequency Bin
        Iout=I  *  (1e26)
    elif OutputUnit=='W/(m**2*m)':                               #Wavelength Bin
        Iout=I  *  (c/(wave**2))
    elif OutputUnit=='W/(m**2*um)':
        Iout=I  *  (c/(wave**2)  *1e-6)    #checked!!!
    elif OutputUnit=='W/(m**2*nm)':
        Iout=I  *  (c/(wave**2)  *1e-9)    #checked!!!
    elif OutputUnit=='W/(m**2*A)':
        Iout=I  *  (c/(wave**2)  *1e-10)    #checked!!!
    # elif OutputUnit=='W/(m**2*DeltaLogLambda800)'                               #Wavelength Bin
    #     R=800
    #     Iout=I  *  (c/(wave)/log(R+1))


    elif OutputUnit=='(ergs/s)/(m**2*A)':
        Iout=I  *  (c/(wave**2)  *1e-10  * 1e7)
    elif OutputUnit=='(ergs/s)/(cm**2*A)':         #checked!!!
        Iout=I  *  (c/(wave**2) *1e-10  * 1e7  *1e-4)

    elif OutputUnit=='W/(m**2*m**-1)':   #= W*m**-1 = W/m          #Wavenumber Bin
        Iout=I  *  (c)
    elif OutputUnit=='W/(cm**2*cm**-1)':   #= W*cm**-1 = W/cm
        Iout=I  *  (c  *  1e-2)   #checked
    elif OutputUnit=='W/(m**2*cm**-1)':
        Iout=I  *  (c  *  1e-2  *  1e4)
    elif OutputUnit=='(photons/s)/(m**2*Hz)':     #-----Photon count------
        Iout=I  *  (1/(h*f))
    elif OutputUnit=='(photons/s)/(m**2*um)':     #-----Photon count------  #checked
        Iout=I  *  (c/(wave**2)  *1e-6  /(h*f))
    elif OutputUnit=='(photons/s)/(m**2*A)':     #-----Photon count------  #checked
        Iout=I  *  (c/(wave**2)  *1e-6  /(h*f)) *1e-4
    elif OutputUnit=='(photons/s)/(cm**2*A)':     #-----Photon count------  #checked
        Iout=I  *  (c/(wave**2)  *1e-6  /(h*f)) *1e-4 *1e-4
    elif OutputUnit=='(photons/s)/(m**2*cm**-1)':     #-----Photon count------ #checked
        Iout=I  *  (c  *  1e-2  *  1e4    /(h*f))
    elif OutputUnit=='(photons/s)/(cm**2*cm**-1)':     #-----Photon count------ #checked
        Iout=I  *  (c  *  1e-2  *  1e4    /(h*f)) *1e-4
    else:
        print('Output Unit Error!!!')


    return Iout   

#%% --- Using emcee to find the optimal ffac and fspot--- #

   

def lnlike(theta, param, fitparanames, spec, models_baseline, T_grid, logg_grid,  model_wv, models_grid):
    """
    log-likelihood function

    Parameters
    ----------
    theta : array
        array of fitted parameter values.
    param : dict
        dictionary of values for all the model parameters.
    # param: Tphot, met, logg, Tspot, fspot, Tfac, ffac
    fitparanames : list of str
        list of the fitted parameters.
    spec : pyTransSpec object 
        (planet atmosphere observations in transmission)
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
    
    # expand models from the baseline models tuple
    model_phot_baseline, model_spot_baseline, model_fac_baseline = models_baseline 
    
    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta[i]
    param = get_derived_param(param)
    
    if "deltaTspot" in fitparanames:
        ind_Tspot = np.argmin(np.abs(T_grid-param["Tspot"]))
        ind_logg_het = np.argmin(np.abs(logg_grid-param["logg_het"]))
        model_spot_fixR  = models_grid[ind_Tspot, ind_logg_het]
    else:
        model_spot_fixR = model_spot_baseline
    
    if "deltaTfac" in fitparanames:
        ind_Tfac = np.argmin(np.abs(T_grid-param["Tfac"]))
        ind_logg_het = np.argmin(np.abs(logg_grid-param["logg_het"]))
        model_fac_fixR  = models_grid[ind_Tfac, ind_logg_het]

    else:
        model_fac_fixR = model_fac_baseline

    if "Tphot" in fitparanames:
        ind_Tphot = np.argmin(np.abs(T_grid-param["Tphot"]))
        ind_logg_phot = np.argmin(np.abs(logg_grid-param["logg_phot"]))
        model_phot_fixR  = models_grid[ind_Tphot,ind_logg_phot]

    else:
        model_phot_fixR = model_phot_baseline
    
    if "Dscale" in fitparanames:
        Dscale = param["Dscale"]
    else:
        Dscale = np.median(spec['yval'])
        

    # pdb.set_trace()
    model_int, _ = calc_stctm_model_and_integrate(Dscale,param['fspot'],param['ffac'], 
                                           spec, model_wv ,model_phot_fixR, model_spot_fixR, model_fac_fixR)
    
    lnlk = -0.5*(np.sum((spec['yval']-model_int )**2./spec["yerrLow"]**2.))
    
    return lnlk, model_int

def lnprior(theta, param, fitparanames, gaussparanames, hyperp_gausspriors, 
                               fitLogfhet, hyperp_logpriors,
            T_grid, logg_grid, Dscale_guess):
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
                                   fitLogfhet=fitLogfhet, hyperp_logpriors=hyperp_logpriors,
                                   T_grid=T_grid,logg_grid=logg_grid, Dscale_guess=Dscale_guess)
    for i, paraname in enumerate(fitparanames):
        if parampriors[paraname][0] < theta[i] <parampriors[paraname][1]:
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
           fitLogfhet, hyperp_logpriors, models_baseline_fixR, T_grid, logg_grid, 
           model_wv,models_grid_fixR):
    
    lp = lnprior(theta, param, fitparanames, gaussparanames, hyperp_gausspriors, 
                 fitLogfhet,hyperp_logpriors,
                 T_grid, logg_grid, Dscale_guess=np.median(spec["yval"]))
    
    if not np.isfinite(lp):
        return -np.inf, np.zeros_like(np.array(spec['yval'])) * np.nan
    else:
        lnlk, model = lnlike(theta, param, fitparanames, spec, models_baseline_fixR, T_grid, 
                             logg_grid, model_wv, models_grid_fixR)
        lnprb = lp + lnlk
        return lnprb, model

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

def BIC(chi2,nDataPoints,nPara):
    '''
    Adapted from auxbenneke/utilities.py
    Liddle, A.R., 2007. Information criteria for astrophysical model selection. Monthly Notices of the Royal Astronomical Society: Letters 377, L74–L78. https://doi.org/10.1111/j.1745-3933.2007.00306.x
    The BIC assumes that the data points are independent and identically distributed, which may or may not be valid depending on the data set under consideration
    '''
    return chi2 + nPara*np.log(nDataPoints)



#%% define default param, fitparanames, parampriors

def init_default_and_fitted_param(Tphot, met, logg_phot,
                                  fitspot=True, fitfac=True, 
                                  fitThet=False, fitTphot=False, fitDscale=False,
                                  fitlogg_phot=False,
                                  fitlogg_het=False,
                                  fitLogfhet=False,
                                  Dscale_guess = 7000, logg_het_guess = None):
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
    fitDscale : bool
        True to fit for the scaling factor D. Default is False.
    fitlogg_het : bool
        True to fit for logg_het. Default is False.
    Dscale_guess : float
        Initial guess for Dscale 
    logg_het_guess : float
        Initial guess for logg_het. Default is None.
    Returns
    -------
    dict of default param values and list of fitted parameters.

    """
    param = dict()
    fitparanames = []
    param["Tphot"] = Tphot
    param["met"] = met
    param['logg_phot'] = logg_phot
    if logg_het_guess is None:
        param["dlogg_het"] =0.
        param['logg_het'] = logg_phot
    else:
        param['logg_het'] = logg_het_guess
        param["dlogg_het"] =logg_het_guess-logg_phot
    param["deltaTspot"] = -150.
    param["deltaTfac"] = 150.
    param["Tspot"] = param["Tphot"] - 150.
    param["Tfac"] = param["Tphot"] + 150.
    param["Dscale"] = Dscale_guess
    

        
    if fitspot:       
        param["fspot"] = 0.02
        if fitLogfhet is False:
            fitparanames.append("fspot")
        else:
            param["log_fspot"] = np.log10(param["fspot"])
            fitparanames.append("log_fspot")

    else:
        param["fspot"] = 0.    
    
    if fitfac:
        param["ffac"] = 0.05
        if fitLogfhet is False:
            fitparanames.append("ffac")
        else:
            param["log_ffac"] = np.log10(param["ffac"])
            fitparanames.append("log_ffac") 
            
    
    else:
        param["ffac"] = 0.
    
    if fitTphot:
        fitparanames.append("Tphot")
        
    if fitThet:
        if fitspot:
            fitparanames.append("deltaTspot")
        if fitfac:
            fitparanames.append("deltaTfac")

    if fitDscale:
        fitparanames.append("Dscale")
        
    
    if fitlogg_phot:
        fitparanames.append("logg_phot")
    if fitlogg_het:
        fitparanames.append("dlogg_het")
        
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
    param["logg_het"] = param["logg_phot"] + param["dlogg_het"]
    if "log_fspot" in param.keys():
        param["fspot"] = 10**param["log_fspot"] 
    if "log_ffac" in param.keys():
        param["ffac"] = 10**param["log_ffac"] 
    return param
    

# def get_prior(allpriors, paraname, default, user_setting):
#     if 
def get_param_priors(param, gaussparanames,hyperp_gausspriors=[], 
                     fitLogfhet=False, hyperp_logpriors=[],
                     T_grid=[2300, 10000], logg_grid=[], Dscale_guess=7000.):
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

    defaultpriors['ffac'] = [0., 0.5]
    defaultpriors['fspot'] = [0., 0.5]
    defaultpriors['deltaTfac'] = [100, T_grid[-1]-param["Tphot"]]
    defaultpriors['deltaTspot'] = [T_grid[0]-param["Tphot"], -100.]
    defaultpriors["Tphot"] = [T_grid[0], T_grid[-1]]
    defaultpriors["Dscale"] = [0.8*Dscale_guess, 1.2*Dscale_guess]
    defaultpriors["logg_phot"] = [np.max([logg_grid[0],2.5]), np.min([5.5, logg_grid[-1]])]
    defaultpriors["dlogg_het"] = [logg_grid[0]-param["logg_phot"], 0.]
    
    parampriors = dict()
    
    # check parameters for log priors
    allparanames = ['ffac',"fspot","deltaTfac","deltaTspot", "Tphot", "Dscale","logg_phot","dlogg_het"]
    for par in allparanames:
        if par not in ["fspot", "ffac"]:
            parampriors[par] = defaultpriors[par]
        else:
            if fitLogfhet:
                lowlim = hyperp_logpriors[0]
                upplim = hyperp_logpriors[1]
                parampriors["log_"+par] = [lowlim,upplim]
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

def format_param_str(param, fitparanames):
    param = get_derived_param(param)

    str1 = "Fitted params: "+str(fitparanames)+"\n"
    str1 = str1 + "Stellar params: Tphot="+str(int(param["Tphot"]*10.)/10.)+" met="+str(param["met"]) + "\n"
    str1 = str1 + "logg_phot="+str(param["logg_phot"]) + " logg_het="+str(param["logg_het"])+ "\n"
    str1 = str1 + "Tspot="+str(int(param["Tspot"]*10.)/10.)+" Tfac="+str(int(param["Tfac"]*10.)/10.)+"\n"
    str1 = str1 + "fspot=" + str(int(param["fspot"]*10000.)/10000.) + " ffac="+str(int(param["ffac"]*10000)/10000.)
    return str1

#%% Make any stellar contamination model

def load_phoenix_models_from_param(param, resPower_target=10000, use_pymsg=True, pymsg_specgrid=None):
    param = get_derived_param(param)

    model_wv, model_spot_fl = load_phoenix_model(param['Tspot'], param["met"], param["logg_het"], resPower_target=resPower_target, use_pymsg=use_pymsg, pymsg_specgrid=pymsg_specgrid)
    _, model_phot_fl = load_phoenix_model(param['Tphot'], param["met"], param["logg_phot"], resPower_target=resPower_target, use_pymsg=use_pymsg, pymsg_specgrid=pymsg_specgrid)
    _, model_fac_fl = load_phoenix_model(param['Tfac'], param["met"], param["logg_het"], resPower_target=resPower_target, use_pymsg=use_pymsg, pymsg_specgrid=pymsg_specgrid)
    
    model_fl_phot_spot_fac = (model_spot_fl,model_phot_fl,model_fac_fl)
    return model_wv, model_fl_phot_spot_fac

    

def generate_PHOENIX_grid_fixedR(wv_target, wv_edges, feh=0.0, Teff_range=[], logg_range=[], 
                          Teffstep=20, loggstep=0.1, fname_save="", use_pymsg = True, pymsg_specgrid=None):
    
    nwave = len(wv_target)
    Teffs_grid = np.arange(Teff_range[0], Teff_range[1], Teffstep)
    loggs_grid = np.arange(logg_range[0], logg_range[1], loggstep)
    
    models_grid = np.zeros((Teffs_grid.size, loggs_grid.size, nwave))
    for i, Teff in enumerate(Teffs_grid):
        if i%10==0: 
            print(str(i)+"/"+str(len(Teffs_grid)))
        for j, logg in enumerate(loggs_grid):
            wvl, flx = load_phoenix_model(Teffs_grid[i], feh, loggs_grid[j], wv_target=wv_target, wv_edges=wv_edges, use_pymsg=use_pymsg, pymsg_specgrid=pymsg_specgrid)
            models_grid[i,j,:] = flx
    h5f = h5py.File(fname_save, 'w')

    h5f.create_dataset('model_grid', data=models_grid)
    h5f.close()
    return Teffs_grid, loggs_grid, models_grid

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
       aio.ascii.write(t_res, results_folder+"stctm_pandas_"+runname+'.csv', format='csv', overwrite=True)
       
    # save best fit parameters and quantiles
    bestfit, ind_bestfit, ind_maxprob = samp2bestfit(panda)
    print(bestfit)
    if save_fit:
        bestfit.to_csv(results_folder+'stctm_bestfit_'+runname+'.csv')
        
    # make corner plot
    print("\nMax Likelihood:\n")
    print(bestfit["MaxLike"])
    
    parabestfit = np.array(bestfit["MaxLike"][2:])
    return bestfit, ind_bestfit, ind_maxprob, parabestfit, samples, t_res

def save_bestfit_stats(spec, ind_bestfit, fitparanames, flat_st_ctm_models, results_folder, runname, save_fit=True):
    # create a dictionary that collates all the best-fit information
    print("Saving stats on the best fit...")
    bestfit_stats = dict()
    best_model = flat_st_ctm_models[ind_bestfit]
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
        aio.ascii.write(t_bestfit_stats, results_folder+"stctm_bestfit_stats_"+runname+'.csv', format='csv', overwrite=True)

#%% Plotting

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
        
def plot_stctm_blobs(spec, stctm_models_blobs,
                          ind_bestfit,ax=None,
                          bestfit_color = 'k', color="b",
                          plot3sig=True, plot2sig=True, plot1sig=True, plotmedian=True,
                          plotbestfit=True, legend_loc=1, save_csv=False, 
                          results_folder="", runname=""):


    # percentiles = [2.5, 16., 50., 84., 97.5]
    percentiles = [0.2, 2.3, 15.9, 50., 84.1, 97.7, 99.8]

    # for each epoch of each planet, get the median, 1 and 2sigma pred time
    stctm_percentiles = np.nanpercentile(stctm_models_blobs, percentiles, axis=0)
    

    # calc for best fit (need ind_bestfit)
    stctm_best = stctm_models_blobs[ind_bestfit, :] 
    if ax is None:
        fig, ax = spec.plot()
    else:
        fig = None
        spec.plot(ax=ax)


    if plot1sig:
        ax.fill_between(spec["wave"],stctm_percentiles[2],stctm_percentiles[4],color=color, alpha=0.5,zorder=-9,label=r'1$\sigma$')


    if plot2sig:
        ax.fill_between(spec["wave"],stctm_percentiles[1],stctm_percentiles[5],color=color, alpha=0.3,zorder=-10,label=r'2$\sigma$')

    if plot3sig:
        ax.fill_between(spec["wave"],stctm_percentiles[0],stctm_percentiles[6],color=color, alpha=0.2,zorder=-10,label=r'3$\sigma$')

    if plotmedian:
        ax.plot(spec["wave"],stctm_percentiles[3],color=color,zorder=-8,label=r'Median')
    
    if plotbestfit: 
        ax.plot(spec["wave"],stctm_best,color=bestfit_color,zorder=-8,label=r'Max. likelihood')

    
    ax.legend(loc=legend_loc)
    
    if save_csv:
        dct = dict()
        dct["wave"] = deepcopy(np.array(spec["wave"]))
        dct["bestfit"] = stctm_best
        dct["median"] = stctm_percentiles[3]
        dct["+1 sigma"] = stctm_percentiles[4]
        dct["-1 sigma"] = stctm_percentiles[2]
        dct["+2 sigma"] = stctm_percentiles[5]
        dct["-2 sigma"] = stctm_percentiles[1]    
        dct["+3 sigma"] = stctm_percentiles[6]
        dct["-3 sigma"] = stctm_percentiles[0] 
        
        t = table.Table(data=dct)
        t.write(results_folder+"blobs_1_2_3_sigma"+runname+".csv", overwrite=True)
    return fig, ax

def get_para_for_postprocess(param, fitparanames, bestfit_params_and_quantiles_table, which):
    """
    which has to be in: ['50', '16', '84', '2.5', '97.5', 'MaxLike']
    returns a dict of temperatures and het. fractions
    """
    para_stcontmodel = dict()
    allpara = ["Tspot" , "Tphot", "Tfac", "fspot", "ffac", "log_fspot", "log_ffac"]
    for p in allpara:
        if p in fitparanames:
            if p[0] == "d":
                para_stcontmodel[p[5:]] =  bestfit_params_and_quantiles_table[which][p] + para_stcontmodel["Tphot"]
            para_stcontmodel[p] =  bestfit_params_and_quantiles_table[which][p]
        elif "delta"+p in fitparanames:
            allpara.append("delta"+p)
        else:
            para_stcontmodel[p] = param[p]

    return para_stcontmodel



def plot_stctm_samples_res(spec, param, fitparanames,
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
            ind_logg_phot = np.argmin(np.abs(logg_grid-param["logg_phot"]))
            ind_logg_het = np.argmin(np.abs(logg_grid-param["logg_het"]))
            
            model_spot = grid_models_fixedresP[ind_Tspot, ind_logg_het]
            model_phot = grid_models_fixedresP[ind_Tphot, ind_logg_phot]
            model_fac = grid_models_fixedresP[ind_Tfac, ind_logg_het]
            
            model_stctm = Dlam_obs(param["Dscale"],param['fspot'],param['ffac'],
                                   model_phot,model_spot,model_fac)
                        
            # convolve the model to the target resolving power
            kernel=Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)    # *2 because FWHM = 2 standard deviation
            l = int(kernel.shape[0]/2)
            # x = wave[l:-l]
            model_stctm = convolve(model_stctm, kernel)[l:-l]
            # pdb.set_trace()
            sample_spectra[i] = deepcopy(model_stctm)
    
    # calculate spectrum for best-fit parameters
    theta = post_burnin_samples.iloc[ind_bestfit]
    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta[i]
    param = get_derived_param(param)
    
    ind_Tspot = np.argmin(np.abs(T_grid-param["Tspot"]))
    ind_Tfac = np.argmin(np.abs(T_grid-param["Tfac"]))
    ind_Tphot = np.argmin(np.abs(T_grid-param["Tphot"]))
    ind_logg_phot = np.argmin(np.abs(logg_grid-param["logg_phot"]))
    ind_logg_het = np.argmin(np.abs(logg_grid-param["logg_het"]))
    
    model_spot = grid_models_fixedresP[ind_Tspot, ind_logg_het]
    model_phot = grid_models_fixedresP[ind_Tphot, ind_logg_phot]
    model_fac = grid_models_fixedresP[ind_Tfac, ind_logg_het]
    
    model_stctm = Dlam_obs(param["Dscale"],param['fspot'],param['ffac'],
                           model_phot,model_spot,model_fac)
    
    # convolve the model to the target resolving power
    kernel=Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)    # *2 because FWHM = 2 standard deviation
    l = int(kernel.shape[0]/2)
    # x = wave[l:-l]
    stctm_best = convolve(model_stctm, kernel)[l:-l]
    
    # calculate percentiles from sample spectra
    perc_spectra = get_stctm_blobs(sample_spectra, percentiles=[0.2, 2.3, 15.9, 50., 84.1, 97.7, 99.8])
    
    
    
    if ax is None:
        fig, ax = spec.plot()
    else:
        fig = None
        spec.plot(ax=ax)

    lowest_z = 1000

    if plot1sig:
        ax.fill_between(wv_array,perc_spectra[2],perc_spectra[4],color=color, alpha=0.5,zorder=lowest_z+1,label=r'1$\sigma$')


    if plot2sig:
        ax.fill_between(wv_array,perc_spectra[1],perc_spectra[5],color=color, alpha=0.3,zorder=lowest_z,label=r'2$\sigma$')

    if plot3sig:
        ax.fill_between(wv_array,perc_spectra[0],perc_spectra[6],color=color, alpha=0.2,zorder=lowest_z,label=r'3$\sigma$')

    if plotmedian:
        ax.plot(wv_array,perc_spectra[3],color=color,zorder=lowest_z+2,label=r'Median')
    
    if plotbestfit: 
        ax.plot(wv_array,stctm_best,color=bestfit_color,zorder=lowest_z+2,label=r'Max. likelihood')

    
    ax.legend(loc=legend_loc)
    
    if save_csv:
        dct = dict()
        dct["wave"] = deepcopy(wv_array)
        dct["bestfit"] = stctm_best
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

    
    # return fig, ax, sample_spectra

def get_closest_models_from_grid(param, grid_models_fixedresP, T_grid, logg_grid):
    ind_Tspot = np.argmin(np.abs(T_grid-param["Tspot"]))
    ind_Tfac = np.argmin(np.abs(T_grid-param["Tfac"]))
    ind_Tphot = np.argmin(np.abs(T_grid-param["Tphot"]))
    ind_logg_phot = np.argmin(np.abs(logg_grid-param["logg_phot"]))
    ind_logg_het = np.argmin(np.abs(logg_grid-param["logg_het"]))
    
    model_spot = grid_models_fixedresP[ind_Tspot, ind_logg_het]
    model_phot = grid_models_fixedresP[ind_Tphot, ind_logg_phot]
    model_fac = grid_models_fixedresP[ind_Tfac, ind_logg_het]
    fl_phot_spot_fac = (model_phot, model_spot, model_fac)
    return fl_phot_spot_fac
    
def get_stctm_blobs(stctm_models_blobs, percentiles=[0.2, 2.3, 15.9, 84.1, 97.7, 99.8]):
 
    
    # for each epoch of each planet, get the median, 1 and 2sigma pred time
    stctm_percentiles = np.nanpercentile(stctm_models_blobs, percentiles, axis=0)
    return stctm_percentiles

def get_stctm_amplitude(spec, stctm_percentiles):
    med = np.median(spec['yval'])
    stctm_percentiles = stctm_percentiles - med
    ampl_3sig = np.max(np.array([np.abs(stctm_percentiles[0]), np.abs(stctm_percentiles[-1])]), axis=0)
    ampl_2sig = np.max(np.array([np.abs(stctm_percentiles[1]), np.abs(stctm_percentiles[-2])]), axis=0)
    ampl_1sig = np.max(np.array([np.abs(stctm_percentiles[2]), np.abs(stctm_percentiles[-3])]), axis=0)
    return ampl_1sig, ampl_2sig, ampl_3sig

def plot_stctm_amplitude(spec, stctm_models_blobs,
                          ax=None,color="b",
                          plot3sig=True, plot2sig=True, plot1sig=True):
    
    if ax is None:
        fig, ax = plt.subplots(1,1)
    else:
        fig = None
    stctm_percentiles = get_stctm_blobs(stctm_models_blobs)
    ampl_1sig, ampl_2sig, ampl_3sig = get_stctm_amplitude(spec, stctm_percentiles)
    zeros = np.zeros_like(ampl_1sig)

    if plot1sig:
        ax.fill_between(spec["wave"],ampl_1sig,zeros,color=color, alpha=0.5,zorder=-9,label=r'1$\sigma$')

    if plot2sig:
        ax.fill_between(spec["wave"],ampl_2sig,zeros,color=color, alpha=0.3,zorder=-10,label=r'2$\sigma$')

    if plot3sig:
        ax.fill_between(spec["wave"],ampl_3sig,zeros,color=color, alpha=0.2,zorder=-10,label=r'3$\sigma$')

    ax.legend(loc=1)
    ax.set_xlabel(r'Wavelength [$\mu m$]')
    ax.set_ylabel(r'|$A_\mathrm{contam}|$ [ppm]')

    return fig, ax

# corner plot function
def plot_corner(samples, plotparams, plot_datapoints=False, smooth=1.,
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
    try:
        fig = corner.corner(samples, labels=plotparams["labels"], 
                            plot_datapoints=plot_datapoints, smooth=smooth,
                            show_titles=show_titles, quantiles=quantiles,
                            title_kwargs=title_kwargs, color=color,
                            hist_kwargs=hist_kwargs, range=rg, 
                            levels=levels, **kwargs)
    except:
        print("Error in corner plot!")
        pdb.set_trace()
    return fig

def plot_custom_corner(samples, fitparanames, parabestfit):
    # reorder samples
    ordered_samples = np.zeros_like(samples)
    ordered_fitparanames_all = ["fspot", "log_fspot","deltaTspot", "ffac", "log_ffac","deltaTfac", "dlogg_het", "Tphot", "logg_phot","Dscale"]
    ordered_fitparanames = []
    for p in ordered_fitparanames_all:
        if p in fitparanames:
            ordered_fitparanames.append(p)
    ordered_parabestfit = np.zeros_like(parabestfit)
    for i, p in enumerate(ordered_fitparanames):
        ind = np.where(np.array(fitparanames)==p)
        # print(ind)
        ordered_samples[:,i] = np.array(samples)[:, ind[0][0]]
        ordered_parabestfit[i] = parabestfit[ind[0][0]]
        

    plotparams = dict()
    plotparams["hist_color"] = "coral"
    plotparams["labels"] = get_labels_from_fitparanames(ordered_fitparanames)
    # pdb.set_trace()

    fig = plot_corner(ordered_samples, plotparams, smooth=0.2, fill_contours=True, 
                      truths = ordered_parabestfit, truth_color="k")
    
                                    

    fig.set_dpi(50)  
    fig.set_figheight(11)
    fig.set_figwidth(15)   
    
    return fig

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
            labels.append(r"$\Delta T_\mathrm{fac}$")           
        elif p == "deltaTspot":
            labels.append(r"$\Delta T_\mathrm{spot}$")    
        elif p == "Tphot":
            labels.append(r"$T_\mathrm{phot}$") 
        elif p == "Dscale":
            labels.append(r"$D$") 
        elif p == "logg_het":
            labels.append(r"log $g_\mathrm{het}$") 
        elif p == "dlogg_het":
            labels.append(r"$\Delta$ log $g_\mathrm{het}$") 
        elif p == "logg_phot":
            labels.append(r"log $g_\mathrm{phot}$") 
    return labels
             
def plot_maxlike_and_maxprob(spec, param, parabestfit, ind_maxprob, ind_bestfit, fitparanames, flat_st_ctm_models, pad=0.25):
    # Plot best fit stellar contamination model
    fig, ax = spec.plot()
    param_bestfit = deepcopy(param)
    for i, p in enumerate(fitparanames):
        param_bestfit[p] = parabestfit[i]
    ax.scatter(spec['wave'], flat_st_ctm_models[ind_maxprob], label="Max. Probability", color="slateblue", alpha=1.)
    ax.scatter(spec['wave'], flat_st_ctm_models[ind_bestfit], label="Max. Likelihood", color="r", alpha=0.5, marker=".", s=50)

    ax.text(np.min(spec["waveMin"])-pad/4, 1.08*np.median(spec['yval']), format_param_str(param_bestfit, fitparanames), fontsize=10, fontweight="bold", color="k")

    ax.set_xlim(np.min(spec["waveMin"])-pad/2, np.max(spec["waveMax"])+pad)
    ax.set_ylim(0.8*np.median(spec['yval']), 1.15*np.median(spec['yval']))
    
    ax.set_title("Best-fit model")
    ax.legend(loc=4)
    return fig, ax   

def get_and_plot_residual_spec(spec, flat_st_ctm_models, ind_bestfit):
    best_fit_dppm = flat_st_ctm_models[ind_bestfit]
    
    residual_spec = deepcopy(spec)
    for i in range(len(residual_spec["yval"])):
        residual_spec[i]['yval'] = residual_spec[i]['yval'] - best_fit_dppm[i]
    
    waveunit = "um"
    residual_spec['wave'].unit=waveunit
    residual_spec['waveMin'].unit=waveunit
    residual_spec['waveMax'].unit=waveunit
    #Meta data
    residual_spec.meta['waveunit']=waveunit
    
    # Plot best fit stellar contamination model
    fig, ax = residual_spec.plot(pretty=True)
    return fig, ax, residual_spec

def get_and_plot_bestfitcorr_spec(spec, flat_st_ctm_models, bestfit, ind_bestfit, pad=0.25):
    # get the best-fitting stellar contamination spectrum
    best_fit_stctm_spec = flat_st_ctm_models[ind_bestfit]
    
    
    # divide this best-fitting stellar contamination spectrum by the best-fitting Dscale to get e_lambda
    if "Dscale" in list(bestfit.index):
        ind_Dscale = np.where(np.array(list(bestfit.index)) == "Dscale")[0][0]
        best_fit_Dscale = bestfit["MaxLike"][ind_Dscale]
        best_fit_elambda = best_fit_stctm_spec/best_fit_Dscale
    else:
        best_fit_elambda = best_fit_stctm_spec
    
    # divide the observed spectrum by e_lambda
    corrected_spec = deepcopy(spec)
    for i in range(len(corrected_spec["yval"])):
        corrected_spec[i]['yval'] = corrected_spec[i]['yval']/best_fit_elambda[i]
    
    waveunit = "um"
    corrected_spec['wave'].unit=waveunit
    corrected_spec['waveMin'].unit=waveunit
    corrected_spec['waveMax'].unit=waveunit
    #Meta data
    corrected_spec.meta['waveunit']=waveunit
    corrected_spec.meta['color'] = "k"
    
    
    # Plot best fit stellar contamination model
    fig, ax = corrected_spec.plot(pretty=True)
    ax.set_xlim(np.min(spec["waveMin"])-pad/2, np.max(spec["waveMax"])+pad)
    ax.set_ylim(0.8*np.median(corrected_spec['yval']), 1.15*np.median(corrected_spec['yval']))
    ax.set_title("Corrected by the max. likelihood model")
    return fig, ax, corrected_spec

def get_and_plot_stctmcorr_spec(spec, flat_st_ctm_models, samples, fitparanames, nsamp=1000, pad=0.25):
    """
    get a spec object that contains the spectrum corrected 
    by the full posterior of stellar contamination spectra
    and plot the result
    
    pad: for wave range
    """
    
    spec_stctmcorr = deepcopy(spec)
    
    random_indices = np.random.randint(0,flat_st_ctm_models.shape[0],nsamp)
    
    random_corr_spectra = np.zeros((nsamp, flat_st_ctm_models[0].shape[0]))
    
    yvals = np.array([spec_stctmcorr[i]["yval"] for i in range(len(spec_stctmcorr["yval"]))])
    yerrUpps = np.array([spec_stctmcorr[i]["yerrUpp"] for i in range(len(spec_stctmcorr["yval"]))])
    yerrLows = np.array([spec_stctmcorr[i]["yerrLow"] for i in range(len(spec_stctmcorr["yval"]))])
    
    for i in range(nsamp): # iterate to get distribution
        # get st ctm model
        this_model = flat_st_ctm_models[random_indices[i]]
        this_Dlam = samples["Dscale"][random_indices[i]]
        
        # transform to e_lambda
        this_elambda = this_model/this_Dlam
        
        # correct yval with this
        # record this sample yval
        random_corr_spectra[i,:] = yvals/this_elambda
        
    
    # get the median and standard deviation in each wavelength bin
    median_corr_yval = np.nanmedian(random_corr_spectra, axis=0)
    std_corr_yval = np.nanstd(random_corr_spectra, axis=0)
    
    # add standard deviation in quadrature with original yErr to get new error bars
    full_yerrUpp = np.sqrt(std_corr_yval**2. + yerrUpps**2.)
    full_yerrLow = np.sqrt(std_corr_yval**2. + yerrLows**2.)
    
    
    for i in range(len(spec_stctmcorr["yval"])):
        spec_stctmcorr[i]['yval'] = median_corr_yval[i]
        spec_stctmcorr[i]['yerrUpp'] = full_yerrUpp[i]
        spec_stctmcorr[i]['yerrLow'] = full_yerrLow[i]
    
    
    # Plot best fit stellar contamination model
    fig, ax = spec_stctmcorr.plot(pretty=True, label="Corrected spectrum")
    spec.plot(ax=ax, color="blue", zorder=-1000, alpha=0.2, label="Initial spectrum")
    ax.set_xlim(np.min(spec["waveMin"])-pad/2, np.max(spec["waveMax"])+pad)
    ax.set_ylim(0.9*np.median(spec_stctmcorr['yval']), 1.1*np.median(spec_stctmcorr['yval']))
    ax.set_title("Corrected by the distribution of models")
    ax.legend()
    
    return fig, ax, spec_stctmcorr


