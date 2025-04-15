# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:40:05 2019

@author: caroline

Fit a stellar spectrum using combination of stellar models
v3: implement log priors
v2: implement error inflation
v1: implement parallelization
CLEAN_v0: move functions to utilities scripts

"""
## ------------------- Import modules ------------------#
import os

os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
os.environ['CRDS_PATH'] = "/home/caroline/crds_cache"
os.environ['PYSYN_CDBS'] = "/home/caroline/trds"

import numpy as np
import matplotlib.pyplot as plt
import pdb
os.environ["OMP_NUM_THREADS"] = "1"

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
from copy import deepcopy
from multiprocessing import Pool
import scipy.signal as sig
from matplotlib.gridspec import GridSpec  



## User inputs: Spectrum and labels
print("Reading in user inputs...")

# label for the observation
label_obs = "test_visit"

start_from_timeseries = True # whether to start from a time series of spectra, or from a spectrum file
save_median_spectrum = False # whether to save the median spectrum that was created from the time series
path_save_median_spectrum = "../../observations/planetname/planet_outoftransit_spectrum.csv" # path to save the median spectrum to

if start_from_timeseries:
    if label_obs=="test_visit":
        #** IF starting from time series of spectra **#
        # path to stellar model time series (where on your installation to find the data. Example with a Eureka! output format)
        folder_stellar_model_ts = "/home/caroline/GitHub/eureka_results/Results20230822_TRAPPIST_1_d__grp0-5__refpix11_TRAPPIST_1_d_NIRSPEC_PRISM_20221109/Stage3/S3_2024-11-09_std_absfluxcal_run1/ap7_bg8/"
        path_to_stellar_model_ts = folder_stellar_model_ts + "S3_std_absfluxcal_ap7_bg8_SpecData.h5"
else:
    if label_obs=="test_visit":
        path_to_spec = "/home/caroline/GitHub/observations/GJ_3090_b/GJ3090b_spectrum.csv"
        spec_format = "MR_csv" # create a new one in pystellspec.py if needed. Current options: ["MR_csv", "basic"]]

if start_from_timeseries:
    # define what is masked in the time and wavelength domains
    obsmaskpattern="avoidflares"


    if obsmaskpattern == "avoidflares":
        # Masking of time intervals to construct median spectrum
        kern_size = 19 # kernel size for plotted median-filtered light curve. Set to None for no smoothing.
        jd_range_mask = [59892+np.array([0.745,0.77]),59892+np.array([0.795,0.86])] # JD range masked to construct spectrum

        # Mask custom wavelength range (user-defined, e.g. saturation)
        wave_range_mask = []
        # wave_range_mask = [[1.01,1.8]] # example for a case where we want to mask a range
    else:
        print("This obsmaskpattern:", obsmaskpattern, "is not defined!")
        pdb.set_trace()

## User inputs: information on the star

which_star = "TRAPPIST-1"

if which_star == "TRAPPIST-1": #** copy this block of code for another star if needed!
    print("The star is TRAPPIST-1 !") 
    # most are from Agol+2020 (Table 7)
    Mstar       = 0.0898*const.M_sun.value # Mann et al. (2019)
    Rstar       = 0.1192*const.R_sun.value # Agol+20
    Teffstar    = 2566. #K, Agol+20
    Teffstar_err = 26. #K, Agol+20
    feh         = 0.040 #[Fe/H], Gillon+2017
    loggstar    = 5.2396 # cgs, Agol+20
else:
    print("The star:", which_star, "does not have defined parameters!")
    pdb.set_trace()

## User inputs: stellar models grid

label_grid = "PHOENIX_TRAPPIST_1"
use_phoenix = True

if label_grid == "PHOENIX_TRAPPIST_1":
    logg_range = [2.5,5.5] 
    Teff_range = [np.min([2300.-Teffstar, -100.])+Teffstar, Teffstar+1000.] 
    loggstep = 0.1 #cgs
    Teffstep = 20. #K
    resPower_target = 10000
    wv_min_um = 0.2
    wv_max_um = 5.4
    stmodfile = "../../R"+str(resPower_target)+"_model_grids/TRAPPIST_1_pymsg.h5"

else:
    print("The stellar model grid:", label_grid, "is not defined!")
    pdb.set_trace()

## User inputs: MCMC fitting params
optimize_param = False # whether to stop after initial processing instead of running MCMC
parallel = True # whether to run the MCMC in parallel
nsteps= 3000 # number of steps in the MCMC
frac_burnin = 0.6 # fraction of steps to discard as burn-in
burnin = int(nsteps*frac_burnin)

if parallel:
    ncpu = 30 # number of cpus to run the fit on
else:
    ncpu = 1 # for serial runs

# options to fit stellar surface and heterogeneity properties
fitspot = False
fitfac = True 
fitThet = True
fitTphot = True
fitlogg_phot = True
fitdlogg_het = True

# data-driven options
fitFscale = True # whether to fit a scaling factor to the model flux
fiterrInfl = True # whether to fit an error inflation factor
save_fit = True # if True, postprocessing will result in the creation of output files

# priors

gaussparanames = np.array([])
hyperp_gausspriors = []

# gaussparanames = np.array(["Tphot"])
# hyperp_gausspriors = [[Teffstar,70]]

fitLogfSpotFac = [0,1] #[spot,fac]
hyperp_logpriors = [-5,0] #[]

# Where to get the photosphere log g value from: user-specified value or pre-set loggstar
logg_phot_source = "value" # options: "value", "loggstar"
logg_phot_value = 2.5 # user-specified value

res_suffix = "_test_for_GitHub"
if parallel:
    res_suffix = res_suffix + "_parallelRun"


##  Start run

# create the label for the run
label_run = ""
label_run = label_run + "Obs_"+label_obs + "_Mod_" + label_grid
if start_from_timeseries:
    label_run = label_run + "_Mask" + obsmaskpattern

print("\nStarting Run:",label_run)

##  Reading in the observed spectrum
if start_from_timeseries:
    print("\nReading in time-series of spectra: ", path_to_stellar_model_ts)

    s3output = h5py.File(path_to_stellar_model_ts, 'r')
    waveUnique = np.unique(np.array(s3output["wave_1d"]))
    timeUnique = np.unique(np.array(s3output["time"]))
    
    print("\nReading in optimal extraction results...")
    optspec = np.array(s3output["optspec"])
    optspec2D = optspec.reshape(timeUnique.size, waveUnique.size) 
    print("\nShape (time, wave) of spectra time series:", optspec2D.shape)


## Remove time intervals from 2D array to calculate median
if start_from_timeseries:
    # create median-normalized light curve
    lc_mednorm = np.nansum(optspec2D/np.nanmedian(optspec2D,axis=0)[None,:],axis=1)
    lc_mednorm = lc_mednorm/np.nanmedian(lc_mednorm)

    if kern_size is not None:
        lc_medfilt = sig.medfilt(lc_mednorm, kernel_size=kern_size)
        t_medfilt = sig.medfilt(timeUnique, kernel_size=kern_size)
    
    fig1, ax = plt.subplots(1,1)
    if kern_size is not None:
        bjd0=int(t_medfilt[0])
        ax.plot(t_medfilt-bjd0,lc_medfilt,color="k")
    else:
        bjd0=int(timeUnique[0])
        ax.plot(timeUnique-bjd0,lc_mednorm,color="k")
    
    
    if len(jd_range_mask[0]):
        for i, mask_range in enumerate(jd_range_mask):
            ind = np.where((timeUnique<mask_range[0])+(timeUnique>mask_range[1]))
            ax.axvspan(mask_range[0]-bjd0,mask_range[1]-bjd0,zorder=-1000,color="gray",alpha=0.5, label=mask_range)
            timeUnique = timeUnique[ind]
            optspec2D = np.squeeze(optspec2D[ind,:])
            print("Masking time range ",mask_range, "new shape (time, wave):",optspec2D.shape)
    
    ax.legend()
    ax.set_xlabel("Time [BJD-"+str(bjd0)+"]")
    ax.set_ylabel("Flux (median-filt.)")
    ax.set_title("Masked regions (TIME)")

##  Calculate median spectrum (or read in spectrum to fit) and save to spec object...
print("\nCalculate or Read in median out-of-transit spectrum and use to initialize pyStellSpec object...")

if start_from_timeseries:
    print("/nCalculate median spectrum and use to initialize StellSpec object...")
    median_spectrum, err_median_spectrum = xtu.get_median_spectrum_and_unc(waveUnique, optspec2D, flux_calibrated=True)

    specdict = dict()
    specdict["wave"] = waveUnique
    specdict["yval"] = median_spectrum*1e10
    specdict["yerrLow"] = err_median_spectrum*1e10
    specdict["yerrUpp"] = err_median_spectrum*1e10
    t = table.Table(data=specdict)
    if save_median_spectrum:
        t.write(path_save_median_spectrum,format="ascii.ecsv")
    spec = psspec.StellSpec(t, inputtype="basic")

else:
    t = aio.ascii.read(path_to_spec) # read in the spectrum that was provided

    spec = psspec.StellSpec(t, inputtype=spec_format)

##  Apply operations to the StellSpec object

print("\nRemoving any wavelengths longer than red end of PHOENIX model range... ")
# make sure we do not go over the limits of the PHOENIX model wavelength range
ind_torem_matchPHOENIX = np.where(np.array(spec.wave)>=5.39)
spec.remDataByIndex(ind_torem_matchPHOENIX[0])

print("\nRemoving any places where the observed flux contains NaNs... ")
# remove bins where the flux is NaN
all_yval= np.array(spec.yval)
ind_NaNs = np.argwhere(np.isnan(all_yval)).T
spec.remDataByIndex(ind_NaNs[0])

print("\nPlotting the spectrum with any discarded wavelength regions ")
# plot the spectrum before wavelength regions masking with regions shown on top
fig2, ax = plt.subplots(1,1)
spec.plot(ax=ax)
for i, mask_wv_range in enumerate(wave_range_mask):
    ind_to_keep = np.where((spec.waveMax<mask_wv_range[0])+(spec.waveMin>mask_wv_range[1]))
    ind_to_mask = np.where((spec.waveMax>=mask_wv_range[0])*(spec.waveMin<=mask_wv_range[1]))
    spec.remDataByIndex(ind_to_mask[0])
    ax.axvspan(mask_wv_range[0],mask_wv_range[1],zorder=-1000,color="gray",alpha=0.5)
    print("Masking wave range ",mask_wv_range, "in median spectrum prior to fit!")
    print("New wave shape:", spec.wave.shape)

ax.set_xlabel(r"Wavelength [$\mu$m]")
ax.set_ylabel(r"Flux (median in time)")
ax.set_title("Masked regions (WAVELENGTH)")

print("\nObtain values from spec object to save as picklable object for MCMC...")
# get the values we need from the spec object while keeping it picklable (only numpy arrays)
waveMin = np.array(spec.waveMin)
waveMax = np.array(spec.waveMax)
yval = np.array(spec.yval)
yerrUpp = np.array(spec.yerrUpp)
spec_pickleformat = (waveMin, waveMax, yval, yerrUpp )


## 
##  Read in the stellar models grid

print("\nReading in grid of stellar models...")
print("Stellar models grid parameter ranges:")
print("logg:", logg_range)
print("Teff:", Teff_range)
# sys.exit()
wv_template_thisR, wv_template_edges_thisR = sru.make_wavelength_array(wv_min_um=wv_min_um, wv_max_um=wv_max_um, 
                    resPower=resPower_target, use_pymsg=True)


# check if the grid file exists
if os.path.exists(stmodfile):
    print("\nStellar model grid file existed ! Reading it in...")
    Teffs_grid = np.arange(Teff_range[0], Teff_range[1], Teffstep)
    loggs_grid = np.arange(logg_range[0], logg_range[1], loggstep)
    
    h5f = h5py.File(stmodfile, 'r')

    models_grid_fixedR = h5f['model_grid'][:]
    h5f.close()    
    print("\nGrid loaded, file closed.")
# if not, generate it
else:
    pdb.set_trace()
    print("\nThe stellar models grid does not exist !! It needs to first be \
          created before running this file with create_fixedR_grid_pymsg.py !!")

    
print("\nFixed resolution model grid shape:", models_grid_fixedR.shape)


##  Visualize nominal model in grid on top of observations (scaled)

# obtain all parameters from default settings
# logg initial guess
if logg_phot_source == "loggstar":
    logg_phot = deepcopy(loggstar)
elif logg_phot_source == "value":
    logg_phot = deepcopy(logg_phot_value)

param, fitparanames = xtu.init_default_and_fitted_param(Teffstar, feh, logg_phot, fitLogfSpotFac=fitLogfSpotFac)
param = xtu.get_derived_param(param)

print("\nGuessing Fscale value from the data...")
Fscale_guess, model_int = xtu.get_scaling_guess(param, Teffs_grid, loggs_grid,  wv_template_thisR, models_grid_fixedR, spec, wave_min_match_um = 2.5, wave_max_match_um = 3.0)

# save the fixed R spectra as baseline specs to do the integration after the fact
fl_phot_spot_fac = xtu.get_closest_models_from_grid(param, models_grid_fixedR, Teffs_grid, loggs_grid)

print("\nPlotting default model with inferred initial guess Fscale applied...")
# Plot model closest to the data
fig3, ax = plt.subplots(1,1,figsize=(10,4))
ax.plot(spec.wave, spec.yval, color="k",zorder=0) # plot the observed spectrum
ax.plot(wv_template_thisR,fl_phot_spot_fac[0]*Fscale_guess,color="r",label="Model",zorder=-1,alpha=0.2)
ax.plot(spec.wave,model_int*Fscale_guess,color="r",label="Model in data bandpass",zorder=-1)
ax.set_xlabel(r"Wavelength [$\mu$m]")
ax.set_ylabel(r'Stellar flux [$\times$ 10$^{-10}$ erg/s/cm$^2$/$\mu$m]')

ax.legend(loc=1)# get nominal model from default parameters

##  Stop here if you just wanted to update your data selection
if optimize_param:
    this_dir = os.getcwd()+"/"
    print("\nSaving previously-created figures...")
    if start_from_timeseries:
        fig1.savefig(this_dir + "exotune_select_time_"+res_suffix+'_preprocessOnly.png')
    fig2.savefig(this_dir + "exotune_select_wave_"+res_suffix+'_preprocessOnly.png')
    plt.close(fig1)
    plt.close(fig2)
    fig3.savefig(this_dir + "exotune_get_fscale_"+res_suffix+'_preprocessOnly.png')
    plt.close(fig3)

    print("\nStopping here, parameters were optimized!")

    sys.exit()

##  Initialize for MCMC
if 1:
    print("\nInitializing params for MCMC...")
    param, fitparanames = xtu.init_default_and_fitted_param(Teffstar, feh, loggstar,
                                                        fitspot=fitspot, fitfac=fitfac,
                                                        fitThet=fitThet, fitTphot=fitTphot,
                                                        fitFscale=fitFscale, fitlogg_phot=fitlogg_phot,
                                                        fitdlogg_het=fitdlogg_het,fiterrInfl=fiterrInfl,
                                                        fitLogfSpotFac=fitLogfSpotFac,
                                                        Fscale_guess = Fscale_guess, 
                                                        dlogg_het_guess = 0.0)
    print("\n** Fitparanames:", fitparanames)
    print("\n** Param: ", param)

    runname = "fit"+label_run
    for p in fitparanames:
        runname = runname + "_"+p
    runname = runname + res_suffix
    res_dir = "../../exotune_results/"+runname+"/"
    print("\nResults folder:", res_dir)
    
    # Setup for post-processing
    if not os.path.isdir(res_dir):
       os.makedirs(res_dir)

    # Save previously-created figures
    print("\nSaving previously-created figures...")
    if start_from_timeseries:
        fig1.savefig(res_dir + "exotune_select_time_"+runname+'.png')
    fig2.savefig(res_dir + "exotune_select_wave_"+runname+'.png')
    plt.close(fig1)
    plt.close(fig2)
    fig3.savefig(res_dir + "exotune_get_fscale_"+runname+'.png')
    plt.close(fig3)

    #** Get .py files used to run this case
    this_dir = os.getcwd()+"/"
    this_script = __file__.split(os.sep)[-1]
    utils_script = str(xtu.__file__)

    script_name = "exotune_runfile_thisrun.py"
    utils_script_name = "exotune_utilities_thisrun.py"
    
    print("\nSaving files...")
    print("\n--Run-analysis file...")
    print("\n**This file:", this_dir+this_script)
    print("Saved to file:", res_dir+script_name)
    shutil.copy(this_dir+this_script, res_dir+script_name)
    print("... Done.**")

    print("\n--Utilities file...")

    print("\n**This file:", utils_script)
    print("Saved to file:", res_dir+utils_script_name)

    shutil.copy(utils_script, res_dir+utils_script_name)
    print("... Done.**")
    
##  Plot fitted spectrum
if 1:
    print("\nPlotting fitted spectrum...")

    # plot spectrum
    fig, ax = spec.plot()    
    fig.savefig(res_dir + "exotune_fitted_spectrum.pdf")
    plt.close(fig)

##  define parameters for emcee run

if 1:
    print("\nSetting up MCMC...")
    ndim, nwalkers = len(fitparanames), len(fitparanames) * 20
    print("\nUsing N walkers = 20x N para...")

    defaultpara=np.zeros(len(fitparanames))

    for i, p in enumerate(fitparanames):
        defaultpara[i] = param[p]
    
        
    pos = [defaultpara+ 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
    
    # dtype and names of blobs
    dtype = [("oot_spec_model", object)]
    
    # define sampler
    # take fixed-R baseline spectra as inputs
    

    other_args = (spec_pickleformat, param, fitparanames, gaussparanames, hyperp_gausspriors, 
                  fitLogfSpotFac, hyperp_logpriors,
                  fl_phot_spot_fac, Teffs_grid, loggs_grid, 
                  wv_template_thisR, models_grid_fixedR,Fscale_guess)




##  Variables and functions required for parallel runs

#** a little awkward to have in this script, but necessary to avoid huge overheads
#** when doing the multiprocessing 
#** (following this guide: https://emcee.readthedocs.io/en/stable/tutorials/parallel/)

if ncpu>1:
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
    hyperp_gausspriors0 = hyperp_gausspriors
    fitLogfSpotFac0 = fitLogfSpotFac
    hyperp_logpriors0 = hyperp_logpriors
    
    
    def lnlike_parallel(theta):
        """
        Log-likelihood function for parallel runs.
        The only explicit parameter is theta (parameter array) but the others
        are set as variables outside of the function.
        """
        
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
        
        if "logFscale" in fitparanames:
            Fscale = param["Fscale"]
        else:
            Fscale=1
            
    
            
        model_int = xtu.calc_spectrum_model_and_integrate(Fscale, param['fspot'],param['ffac'], spec, 
                                                      model_wv, model_phot_fixR, model_spot_fixR, 
                                                      model_fac_fixR)
    
    
        yerrLow = yerrLow*param["errInfl"]

        lnlk = -0.5*(np.sum((yval-model_int )**2./yerrLow**2.))
        
        return lnlk, model_int
    
    
    
    def get_param_priors_parallel():
        """
        Get parameter priors for the MCMC run.
        Returns a dictionary of parameter priors in the form of [low, high] for uniform priors.

        """
        T_grid =  T_grid0
        logg_grid =  logg_grid0
        param =  param0
        gaussparanames =  gaussparanames0
        hyperp_gausspriors = hyperp_gausspriors0
        fitLogfSpotFac = fitLogfSpotFac0
        hyperp_logpriors = hyperp_logpriors0
        
        defaultpriors = dict()

        defaultpriors['ffac'] = [0., 0.9]
        defaultpriors['fspot'] = [0., 0.9]
        defaultpriors['deltaTfac'] = [100, T_grid[-1]-param["Tphot"]]
        defaultpriors['deltaTspot'] = [T_grid[0]-param["Tphot"], -100.]
        defaultpriors["Tphot"] = [T_grid[0], T_grid[-1]]
        defaultpriors["logFscale"] = [np.log10(Fscale_guess)-5, np.log10(Fscale_guess)+5]
        defaultpriors["logErrInfl"] = [0,np.log10(50)]
        defaultpriors["loggphot"] = [np.max([logg_grid[0],2.5]), np.min([5.5, logg_grid[-1]])]
        defaultpriors["dlogghet"] = [logg_grid[0]-param["loggphot"], 0.]
        
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
                    

        # bounds for parameters with Gaussian priors
        for par in param:
            if par in gaussparanames:
                ind_gaussparam = np.where(gaussparanames == par)[0][0]
                mean_gaussparam = hyperp_gausspriors[ind_gaussparam][0]
                std_gaussparam = hyperp_gausspriors[ind_gaussparam][1]
                new_prior = [np.max([mean_gaussparam - 5.*std_gaussparam, parampriors[par][0]]), np.min([mean_gaussparam + 5.*std_gaussparam, parampriors[par][1]])]
                parampriors[par] = new_prior
                
        return parampriors
    
    

    
    
    def lnprior_parallel(theta):
        """
        Log-prior function for parallel runs.
        The only explicit parameter is theta (parameter array) but the others
        are set as variables outside of the function.
        """
    
        param =  param0
        fitparanames =  fitparanames0
        gaussparanames =  gaussparanames0
        hyperp_gausspriors = hyperp_gausspriors0


        lp = 0.
        
        for i, paraname in enumerate(fitparanames):
            param[paraname] = theta[i]
        param = xtu.get_derived_param(param)
        
        parampriors = get_param_priors_parallel()

        for i, paraname in enumerate(fitparanames):
            if parampriors[paraname][0] < theta[i] <parampriors[paraname][1]:
                if paraname in gaussparanames:
                    ind_gaussparam = np.where(gaussparanames == paraname)[0][0]
                    mean_gaussparam = hyperp_gausspriors[ind_gaussparam][0]
                    std_gaussparam = hyperp_gausspriors[ind_gaussparam][1]
                    lp = lp - 0.5*((theta[i]-mean_gaussparam)**2./(std_gaussparam**2.))
    
                else:
                    lp = lp + 0.0
            else:
                lp = - np.inf
        
        if param["fspot"] + param["ffac"]>1:
            lp = - np.inf             
    
        return lp
    
    def lnprob_parallel(theta):
        """
        Log-probability function for parallel runs.
        The only explicit parameter is theta (parameter array) but the others
        are set as variables outside of the function.
        """
        waveMin, waveMax, yval, yerrLow = spec0
        
        lp = lnprior_parallel(theta)
        
        if not np.isfinite(lp):
            return -np.inf, np.zeros_like(np.array(yval)) * np.nan
            # return -np.inf
        else:
            lnlk, model = lnlike_parallel(theta)
    
            lnprb = lp + lnlk
            return lnprb, model
    
## 
if 1:
    print("\nCreating sampler...")
    
    if ncpu==1:
        print("Running serial version of the MCMC fit!")
        sampler = emcee.EnsembleSampler(nwalkers, ndim, xtu.lnprob, 
                                        args=other_args,
                                        blobs_dtype=dtype)
    else:
        print("Running parallel version of the MCMC fit on", ncpu, "CPUs!")
        pool = Pool(ncpu) # on 10 CPU

        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_parallel,
                                        pool=pool,
                                        blobs_dtype=dtype)
        
##  run emcee
if 1:
    print("\nRunning MCMC...")

    sampler.run_mcmc(pos, nsteps, progress=True, store=True)


    print("\nRecall -- MCMC setup:")
    print("\n** Fitparanames:", fitparanames)
    print("\n** Param: ", param)
    
##  Post-process + Get blobs

if 0:
    sys.exit()
##  Set up MCMC

##  ----- Post-processing ----
if 1: # chainplot

    print("\nGenerating chain plot...")

    fig, axes = xtu.chainplot(sampler.chain[:, burnin:, :], labels=fitparanames)
    if save_fit:
        fig.savefig(res_dir + "exotune_chainplot_noburnin_"+runname+'.png')
        plt.close(fig)

    fig, axes = xtu.chainplot(sampler.chain[:, :, :], labels=fitparanames)
    if save_fit:
        fig.savefig(res_dir + "exotune_chainplot_"+runname+'.png')
        plt.close(fig)

## 
if 1: #saving
    print("\nSaving results to pandas DataFrame...")

    rs = xtu.save_mcmc_to_pandas(res_dir, runname, sampler, burnin, ndim, 
                                 fitparanames, save_fit)
    bestfit, ind_bestfit, ind_maxprob, parabestfit, samples, t_res = rs

## 
##  Setup for plotting

pad = 0.25 # in um
target_resP = 300 # plot stellar spectra at which resolution
sample_spectra = None

##  plot 1,2,3 sigma from fixed resolution contamination spectra

print("\nPlotting 1,2,3 sigma percentiles for the models superimposed with the data...")

if sample_spectra is None:
    # spec.meta["color"] = "k"

    fig, ax, sample_spectra = xtu.plot_exotune_samples_res(spec, param, fitparanames,
                              ind_bestfit, samples, Teffs_grid, loggs_grid,
                              wv_template_thisR,
                              models_grid_fixedR,
                              sample_spectra=None, modelgrid_resP=10000,
                              target_resP=target_resP,N_samp=1000, ax=None,
                              bestfit_color = 'k', color="coral",plot3sig=True,
                              plot2sig=True, plot1sig=True, plotmedian=True,
                              plotbestfit=True, legend_loc=4, save_csv=True, 
                              results_folder=res_dir, runname=runname)


    ax.set_xlim(np.min(spec.waveMin)-pad/2, np.max(spec.waveMax)+pad)
    # ax.set_ylim(0.8*np.median(spec.yval), 1.15*np.median(spec.yval))
    
    if save_fit:
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.png')
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.pdf')
    xtu.xspeclog(ax,level=1)
    if save_fit:
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.png')
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.pdf')
        plt.close(fig)


else:
    # when sample spectra were already calculated:
    # spec.meta["color"] = "k"
    fig, ax, _ = xtu.plot_exotune_samples_res(spec, param, fitparanames,
                              ind_bestfit, samples, Teffs_grid, loggs_grid,
                              wv_template_thisR,
                              models_grid_fixedR,
                              sample_spectra=sample_spectra, modelgrid_resP=10000,
                              target_resP=target_resP,N_samp=1000, ax=None,
                              bestfit_color = 'k', color="coral",plot3sig=True,
                              plot2sig=True, plot1sig=True, plotmedian=True,
                              plotbestfit=True, legend_loc=4, save_csv=True) 
    ax.set_xlim(np.min(spec.waveMin)-pad/2, np.max(spec.waveMax)+pad)

    if save_fit:
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.png')
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.pdf')
    xtu.xspeclog(ax,level=1)
    if save_fit:
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.png')
        fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.pdf')
        plt.close(fig)

##  Save blobs + plot best fit
if 1:
    print("\nSaving blobs...")
    # Blobs
    blobs = sampler.get_blobs()
    if ncpu>1:
        # Close the thread pool
        pool.close()
        pool.join()
    flat_oot_spec_models = blobs.T[:, burnin:]["oot_spec_model"].reshape((-1))
    if save_fit:
        np.save(res_dir+"oot_spec_model_blobs_"+runname+".npy", flat_oot_spec_models)

    print("\nPlotting best fit model with observations...")
    fig, ax = xtu.plot_maxlike_and_maxprob(spec, param, parabestfit, ind_maxprob, 
                                           ind_bestfit, fitparanames, flat_oot_spec_models, pad=pad)
    
    if save_fit:
        fig.savefig(res_dir+"exotune_bestfit_model_with_obs_"+runname+".png")
        plt.close(fig)

    print("\nComputing statistics (chi-squared, BIC) on the run results...")
    xtu.save_bestfit_stats(spec, ind_bestfit, fitparanames, flat_oot_spec_models,
                           res_dir, runname, save_fit=save_fit)
    
            
    print("\nSaving default parameters to file...")
    t_defaultparam = table.Table([param])
    if save_fit:
        aio.ascii.write(t_defaultparam, res_dir+"exotune_defaultparams_"+runname+'.csv', format='csv', overwrite=True)

    oot_spec_models = np.array([flat_oot_spec_models[i] for i in range(flat_oot_spec_models.size)])

## plot 1,2,3 sigma with blobs

if 1:
    fig, ax = xtu.plot_exotune_blobs(spec, oot_spec_models,
                                   ind_bestfit,
                                   bestfit_color='k', color="coral",
                                   plot3sig=True, plot2sig=True, plot1sig=True, plotmedian=True,
                                   plotbestfit=True, legend_loc=4, save_csv=True,
                                   results_folder=results_folder, runname=runname)

    ax.set_xlim(np.min(spec["waveMin"]) - pad / 2, np.max(spec["waveMax"]) + pad)
    ax.set_ylim(0.8 * np.median(spec['yval']), 1.15 * np.median(spec['yval']))

    if save_fit:
        fig.savefig(results_folder + "exotune_1_2_3_sigma_" + runname + ".pdf")

    if save_fit:
        fig.savefig(results_folder + "exotune_1_2_3_sigma_" + runname + ".png")

# sys.exit()
##  joint plot
if 1:
    print("\nCreating combo plot with spectra and parameter distributions...")
    nparaplot = len(fitparanames)
    if "logErrInfl" in fitparanames:
        nparaplot = nparaplot -1
    if "logFscale" in fitparanames:
        nparaplot = nparaplot -1
    fig = plt.figure(figsize=(10,8))
    gs = GridSpec(2, nparaplot,left=0.1, hspace=0.25,right=0.95,bottom=0.1,top=0.95,height_ratios=[5,1.5])
    axspec = fig.add_subplot(gs[0,:])

    # plot the spectrum
    xtu.plot_exotune_samples_res(spec, param, fitparanames,
                              ind_bestfit, samples, Teffs_grid, loggs_grid,
                              wv_template_thisR,
                              models_grid_fixedR,
                              sample_spectra=sample_spectra, ax=axspec, modelgrid_resP=10000,
                              target_resP=target_resP,N_samp=1000, 
                              bestfit_color = 'k', color="C1",plot3sig=True,
                              plot2sig=True, plot1sig=True, plotmedian=True,
                              plotbestfit=True, legend_loc=4, save_csv=True) 
    
    axspec.set_xlabel("Wavelength [$\mu$m]")
    axspec.set_ylabel(r'Stellar flux [$\times$ 10$^{-10}$ erg/s/cm$^2$/$\mu$m]')

    
    axspec.set_xlim(np.min(spec.waveMin)-pad/2, np.max(spec.waveMax)+pad)
    plotlabels = xtu.get_labels_from_fitparanames(fitparanames)
    
    
    iplot = 0
    # get parameter priors to use as bounds for the distributions
    if ncpu>1:
        parampriors = get_param_priors_parallel()
    else:
        parampriors = xtu.get_param_priors(param,[],Fscale_guess=Fscale_guess)

    # bottom panels: distributions on the parameters
    for i in range(len(fitparanames)):
        if fitparanames[i]!="logErrInfl":
            if fitparanames[i]!="logFscale":
                axi = fig.add_subplot(gs[1,iplot])
                
                axi.hist(samples[fitparanames[i]], ec="k", color="C1", alpha=0.7, histtype="stepfilled",bins=20)
                # axi.set_xticklabels(axi.get_xticklabels(),fontsize=10)
                axi.tick_params(axis='x', labelsize=10)
                axi.set_yticks([])
                # axi.set_yticklabels([])
                axi.set_xlabel(plotlabels[i])
                
                min_samples = np.percentile(samples[fitparanames[i]], 0.05)
                max_samples = np.percentile(samples[fitparanames[i]], 99.95)
                low_range = np.max([min_samples, parampriors[fitparanames[i]][0]])
                upp_range = np.min([max_samples, parampriors[fitparanames[i]][1]])
                perc_16 = np.percentile(samples[fitparanames[i]], 16)
                perc_50 = np.percentile(samples[fitparanames[i]], 50)
                perc_84 = np.percentile(samples[fitparanames[i]], 84)
                m1sig = perc_50 - perc_16  
                p1sig = perc_84 - perc_50
                
                perc_50_str = str(int(perc_50*1000)/1000)
                m1sig_str = str(int(m1sig*1000)/1000)
                p1sig_str = str(int(p1sig*1000)/1000)
                str_label = plotlabels[i] +"= "+perc_50_str+"$^{+"+p1sig_str+"}_{-"+m1sig_str+"}$"
                axi.set_title(str_label,fontsize=10)
                
                axi.set_xlim([low_range, upp_range])
                iplot += 1
    # ax.set_ylim(0.8*np.median(spec.yval), 1.15*np.median(spec.yval))
    
    if save_fit:
        fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.png')
        fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.pdf')
    xtu.xspeclog(axspec,level=1)
    if save_fit:
        fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.png')
        fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.pdf')
        plt.close(fig)

##  Corner plot
# try: #corner plot
print("\nCreating corner plot...")
fig = xtu.plot_custom_corner(samples, fitparanames, parabestfit, param,
                   gaussparanames,hyperp_gausspriors,fitLogfSpotFac,hyperp_logpriors,Teffs_grid,loggs_grid)


suffix = "_custom"
if save_fit:
    fig.savefig(res_dir+"exotune_corner_bestfit_"+runname+suffix+".pdf")
    plt.close(fig)
print("\nPost-processing done!")
# except:
    # print("Corner plot failed -- probably dynamic range issue")


