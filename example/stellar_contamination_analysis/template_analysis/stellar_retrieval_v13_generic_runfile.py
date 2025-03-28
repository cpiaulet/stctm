# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 12:40:05 2025

@author: caroline

fit a model assuming a flat planet spectrum and stellar 
contribution (M2V) as from Rackham et al. 2018 with
spots and faculae.

plot the observed spectrum with the best fit + the mean
prediction from Rackham+ (2018)
v11: add option to fit photosphere log g
v12: add possibility to fit fhet in log space
v13: change for ffac, fspot
"""

##------------------- Import modules ------------------#
import os
import numpy as np
import matplotlib.pyplot as plt

# make sure the environment variables are set up
os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
os.environ['CRDS_PATH'] = "/home/caroline/crds_cache"
os.environ['PYSYN_CDBS'] = "/home/caroline/trds"

import stctm.pytransmspec as ptspec
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


##------------------- User inputs ------------------
## --- User inputs: Spectrum and labels --- #
# input required every time

#**#** Information on the fitted spectrum (input data)
which_planet = "TRAPPIST-1 c"
which_star = "TRAPPIST-1"
which_visit = "1"
instrument = "NIRISS SOSS"
# label for plots
label = instrument +" "+which_planet +" visit "+which_visit

# path to spectrum files (replace with the path to your files)
if which_visit == "1":
    path_to_spec = "/home/caroline/GitHub/observations/TRAPPIST_1_c/MR_TRAPPIST_1_c_NIRISS_SOSS_visit1_revResp.spec"
elif which_visit == "2":
    path_to_spec = "/home/caroline/GitHub/observations/TRAPPIST_1_c/MR_TRAPPIST_1_c_NIRISS_SOSS_visit2_revResp.spec"

spec_format = "wavemicrons" # create a new one in pytransspec.py if needed. Current options: ["basic", "wavemicrons"]]


## --- User inputs: Stellar parameters --- #
if which_star == "TRAPPIST-1":  ## copy this block of code for another star if needed!
    print("The star is TRAPPIST-1 !")

    # most are from Agol+2020 (Table 7)
    Mstar = 0.0898 * const.M_sun.value  # Mann et al. (2019)
    Rstar = 0.1192 * const.R_sun.value  # Agol+20
    Teffstar = 2566.  # K, Agol+20
    Teffstar_err = 26.  # K, Agol+20
    feh = 0.040  # [Fe/H], Gillon+2017
    loggstar = 5.2396  # cgs, Agol+20

    print('logg =' + str(loggstar))
else:
    print("Stellar parameters are not defined for the star " + which_star + " !!!")
    pdb.set_trace()

## --- User inputs: MCMC fitting params

nsteps=3000
frac_burnin = 0.6
burnin = int(nsteps*frac_burnin)

# which parameters to fit
fitspot = True
fitfac = True 
fitThet = True
fitTphot = True
fitDscale = True
fitlogg_phot = True
fitlogg_het = True

# whether to save the fit
save_fit=True

# Setting up Gaussian priors, if any

# --- Block to use if there is no Gaussian prior ---

# gaussparanames = np.array([])
# hyperp_gausspriors = []

# --- Block to use if you want to use a Gaussian prior e.g. on the stellar temperature ---
gaussparanames = np.array(["Tphot"])
hyperp_gausspriors = [[Teffstar,70]] # mean and std of the Gaussian prior on Tphot

if len(gaussparanames):
    print("** Using Gaussian prior on ",list(gaussparanames))

# specify whether we want to fit fspot/ffac with prior uniform in log or lin space
fitLogfSpotFac = [0,0] # this means fitting both in linear space
hyperp_logpriors = [-5,0] #[]


# Where to get the photosphere log g value from: user-specified value or pre-set loggstar
logg_phot_source = "loggstar" # options: "value", "loggstar"
logg_phot_value = 5
# where to get the heterogeneity log g value from: user-specified value, or pre-set loggstar, or whatever logg_phot is
logg_het_default_source = "logg_phot" # options: "value", "logg_phot"
logg_het_value = 5

# choose suffix for all the files produced by this fit
res_suffix="TRAPPIST_1_c_example"

## --- User inputs: Read in stellar models grid --- #

#** Modify below only if you generated the grid differently

# stellar models grid file path
if which_star == "TRAPPIST-1":
    # range of params for the grid
    logg_range = [2.5,5.5]
    Teff_range = [np.min([2300.-Teffstar, -100.])+Teffstar, Teffstar+1000.]
    loggstep = 0.1 #cgs
    Teffstep = 20. #K
    resPower_target = 10000
    wv_min_um = 0.2
    wv_max_um = 5.4

    # stellar models grid file path
    stmodfile = "../../R"+str(resPower_target)+"_model_grids/TRAPPIST_1"
# Add other "elif" statement for another star

stmodfile = stmodfile + "_pymsg.h5"

## --- End of user inpute --- #


## --- Load the grid --- #
print("Grid parameter ranges:")
print("logg:", logg_range)
print("Teff:", Teff_range)
# sys.exit()
wv_template_thisR, wv_template_edges_thisR = sru.make_wavelength_array(wv_min_um=wv_min_um, wv_max_um=wv_max_um, 
                    resPower=resPower_target, use_pymsg=True)


# check if the grid file exists
if os.path.exists(stmodfile):
    print("Stellar model grid file existed! Reading it in...")
    Teffs_grid = np.arange(Teff_range[0], Teff_range[1], Teffstep)
    loggs_grid = np.arange(logg_range[0], logg_range[1], loggstep)
    
    h5f = h5py.File(stmodfile, 'r')

    models_grid_fixedR = h5f['model_grid'][:]
    h5f.close()    
    print("Grid loaded, file closed.")
# if not, generate it
else:
    pdb.set_trace()
    print("The stellar models grid does not exist !! It needs to first be \
          created before running this file with create_fixedR_grid_pymsg.py !!")

    
print("Fixed resolution model grid shape:", models_grid_fixedR.shape)


## Get default log g for the heterogeneity and photosphere

# logg for MCMC
if logg_phot_source == "loggstar":
    logg_phot = deepcopy(loggstar)
elif logg_phot_source == "value":
    logg_phot = deepcopy(logg_phot_value)

if logg_het_default_source == "logg_phot":
    logg_het_default = deepcopy(logg_phot)
elif logg_het_default_source == "value":
    logg_het_default = deepcopy(logg_het_value)


## --- Import observed spectrum from file --- #
# input required if the required column names have different names in your input file

spec = ptspec.TransSpec(path_to_spec, inputtype=spec_format)

## Initialize for MCMC
if 1:
    print("Initializing params for MCMC...")

    
    param, fitparanames = sru.init_default_and_fitted_param(Teffstar, feh, logg_phot, 
                                                        fitspot=fitspot, fitfac=fitfac, 
                                                        fitThet=fitThet, fitTphot=fitTphot,
                                                        fitDscale=fitDscale,
                                                        fitlogg_phot=fitlogg_phot,
                                                        fitlogg_het=fitlogg_het,
                                                        fitLogfSpotFac=fitLogfSpotFac,
                                                        Dscale_guess=np.median(spec["yval"]),
                                                        logg_het_guess = logg_het_default)
    
    print("\n** Fitparanames:", fitparanames)
    print("\n** Param: ", param)
    
    # save the fixed R spectra as baseline specs to do the integration after the fact
    fl_targetR_phot_spot_fac_baseline = sru.get_closest_models_from_grid(param, models_grid_fixedR, Teffs_grid, loggs_grid)

    runname = "fit"
    for p in fitparanames:
        runname = runname + "_"+p
    runname = runname + res_suffix
    results_folder = "../../stellar_contamination_results/"+runname+"/"   
    print("Results folder:", results_folder)
    
    # Setup for post-processing

    if not os.path.isdir(results_folder):
       os.makedirs(results_folder)
       
    
    #** Get .py files used to run this case
    this_dir = os.getcwd()+"/"
    res_dir = os.sep.join(__file__.split(os.sep)[:-2]) + "/../stellar_contamination_results/"+runname+"/"
    this_script = __file__.split(os.sep)[-1]
    script_name = "runscript.py"
    print("\nSaving files...")
    print("\nThis file:", this_dir+this_script)
    print("Saved to file:", res_dir+script_name)
    shutil.copy(this_dir+this_script, res_dir+script_name)


## --- Plot observed spectrum --- #

if 1:
    print("Plotting fitted spectrum...")

    # plot spectrum
    fig, ax = spec.plot()    
    fig.savefig(results_folder + "stctm_fitted_spectrum.pdf")

    
## define parameters for emcee run

if 1:
    print("Setting up MCMC...")
    ndim, nwalkers = len(fitparanames), len(fitparanames) * 20

    defaultpara=np.zeros(len(fitparanames))

    for i, p in enumerate(fitparanames):
        defaultpara[i] = param[p]
    
        
    pos = [defaultpara+ 1e-2*np.random.randn(ndim) for i in range(nwalkers)]
    
    # dtype and names of blobs
    dtype = [("st_ctm_model", object)]
    
    # define sampler
    # take fixed-R baseline spectra as inputs
    other_args = (spec, param, fitparanames, gaussparanames, hyperp_gausspriors,
                  fitLogfSpotFac, hyperp_logpriors,
                  fl_targetR_phot_spot_fac_baseline, Teffs_grid, loggs_grid, 
                  wv_template_thisR, models_grid_fixedR)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, sru.lnprob, 
                                    args=other_args,
                                    blobs_dtype=dtype)
    

## run emcee
if 1:
    print("Running MCMC...")

    sampler.run_mcmc(pos, nsteps, progress=True, store=True)


    print("\nRecall -- MCMC setup:")
    print("\n** Fitparanames:", fitparanames)
    print("\n** Param: ", param)
    
## Post-process + Get blobs

if 0:
    sys.exit()
## ----- Post-processing ----
if 1: # chainplot

    
    fig, axes = sru.chainplot(sampler.chain[:, burnin:, :], labels=fitparanames)
    if save_fit:
        fig.savefig(results_folder + "stctm_chainplot_noburnin_"+runname+'.png')

    fig, axes = sru.chainplot(sampler.chain[:, :, :], labels=fitparanames)
    if save_fit:
        fig.savefig(results_folder + "stctm_chainplot_"+runname+'.png')

##
if 1: #saving


    rs = sru.save_mcmc_to_pandas(results_folder, runname, sampler, burnin, ndim, 
                                 fitparanames, save_fit)
    bestfit, ind_bestfit, ind_maxprob, parabestfit, samples, t_res = rs
    
## Stopping point
# sys.exit()

## Setup for plotting

pad = 0.25 # in um
target_resP = 100 # plot stellar spectra at which resolution
sample_spectra = None

## plot 1,2,3 sigma from fixed resolution contamination spectra

if sample_spectra is None:
    spec.meta["color"] = "k"

    fig, ax, sample_spectra = sru.plot_stctm_samples_res(spec, param, fitparanames,
                              ind_bestfit, samples, Teffs_grid, loggs_grid,
                              wv_template_thisR,
                              models_grid_fixedR,
                              sample_spectra=None, modelgrid_resP=10000,
                              target_resP=target_resP,N_samp=1000, ax=None,
                              bestfit_color = 'k', color="coral",plot3sig=True,
                              plot2sig=True, plot1sig=True, plotmedian=True,
                              plotbestfit=True, legend_loc=4, save_csv=True, 
                              results_folder=results_folder, runname=runname)


    ax.set_xlim(np.min(spec["waveMin"])-pad/2, np.max(spec["waveMax"])+pad)
    ax.set_ylim(0.8*np.median(spec['yval']), 1.15*np.median(spec['yval']))
    
    if save_fit:
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"1_2_3sigma_"+runname+'.png')
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"1_2_3sigma_"+runname+'.pdf')
    sru.xspeclog(ax,level=1)
    if save_fit:
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.png')
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.pdf')


else:
    # when sample spectra were already calculated:
    spec.meta["color"] = "k"
    fig, ax, _ = sru.plot_stctm_samples_res(spec, param, fitparanames,
                              ind_bestfit, samples, Teffs_grid, loggs_grid,
                              wv_template_thisR,
                              models_grid_fixedR,
                              sample_spectra=sample_spectra, modelgrid_resP=10000,
                              target_resP=target_resP,N_samp=1000, ax=None,
                              bestfit_color = 'k', color="coral",plot3sig=True,
                              plot2sig=True, plot1sig=True, plotmedian=True,
                              plotbestfit=True, legend_loc=4, save_csv=True) 
    ax.set_xlim(np.min(spec["waveMin"])-pad/2, np.max(spec["waveMax"])+pad)

    if save_fit:
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"1_2_3sigma_"+runname+'.png')
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"1_2_3sigma_"+runname+'.pdf')
    sru.xspeclog(ax,level=1)
    if save_fit:
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.png')
        fig.savefig(results_folder + "stctm_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.pdf')



## Corner plot
if 1: #corner plot
    fig = sru.plot_custom_corner(samples, fitparanames, parabestfit)
    
    suffix = "_custom"
    if save_fit:
        fig.savefig(results_folder+"stctm_corner_bestfit_"+runname+suffix+".pdf")
        
## Save blobs + plot best fit
if 1:    
    # Blobs
    blobs = sampler.get_blobs()
    flat_st_ctm_models = blobs.T[:, burnin:]["st_ctm_model"].reshape((-1))
    if save_fit:
        np.save(results_folder+"st_ctm_model_blobs_"+runname+".npy", flat_st_ctm_models)


    fig, ax = sru.plot_maxlike_and_maxprob(spec, param, parabestfit, ind_maxprob, 
                                           ind_bestfit, fitparanames, flat_st_ctm_models, pad=pad)
    
    if save_fit:
        fig.savefig(results_folder+"stctm_bestfit_model_with_obs_"+runname+".png")

    sru.save_bestfit_stats(spec, ind_bestfit, fitparanames, flat_st_ctm_models,
                           results_folder, runname, save_fit=save_fit)

    t_defaultparam = table.Table([param])
    if save_fit:
        aio.ascii.write(t_defaultparam, results_folder + "stctm_defaultparams_" + runname + '.csv', format='csv',
                        overwrite=True)

    st_ctm_models = np.array([flat_st_ctm_models[i] for i in range(flat_st_ctm_models.size)])

## plot 1,2,3 sigma with blobs

    
if 1: 
    fig, ax = sru.plot_stctm_blobs(spec, st_ctm_models,
                              ind_bestfit,
                              bestfit_color = 'k', color="coral",
                              plot2sig=True, plot1sig=True, plotmedian=True,
                              plotbestfit=True, legend_loc=4, save_csv=True,
                              results_folder=results_folder, runname=runname)

    ax.set_xlim(np.min(spec["waveMin"])-pad/2, np.max(spec["waveMax"])+pad)
    ax.set_ylim(0.8*np.median(spec['yval']), 1.15*np.median(spec['yval']))

    if save_fit:
        fig.savefig(results_folder+"stctm_1_2_sigma_noprint"+runname+".pdf")

    if save_fit:
        fig.savefig(results_folder+"stctm_1_2_sigma_"+runname+".png")

sys.exit()  


## Plot amplitude

## plot 1,2,3 sigma with the amplitude of stellar contamination
if 1: # v1: for the paper

    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1)
    ax = fig.add_subplot(gs[:2, :])
    
    sru.plot_stctm_blobs(spec, st_ctm_models,
                              ind_bestfit,ax=ax, 
                              bestfit_color = 'k', color="coral",
                              plot2sig=True, plot1sig=True, plotmedian=True,
                              plotbestfit=True, legend_loc=4)
    ax.set_ylabel(r'Transit Depth [ppm]')
    sru.xspeclog(ax,level=1)
    ax.axhline(np.median(spec["yval"]), ls="--", color="k", label="Flat spectrum")
    ax.legend(ncol=2)

    # ax.set_ylim(7000., 8200.)
    # ax.set_xlim(0.85, 5.1)
    # ax.set_ylim(7330., 7730.)

    ax2 = fig.add_subplot(gs[2, :], sharex=ax)
    sru.plot_stctm_amplitude(spec, st_ctm_models,
                              ax=ax2,color="coral")
    ax2.axhline(100., color="grey", ls=":", label="100 ppm")
    sru.xspeclog(ax2,level=1)

    ax2.legend()
    fig.tight_layout()

    ax2.set_ylim(ymin=0)
    if save_fit:
        fig.savefig(results_folder+"stctm_1_2_sigma_withamplitude_noprint"+runname+".pdf")




































