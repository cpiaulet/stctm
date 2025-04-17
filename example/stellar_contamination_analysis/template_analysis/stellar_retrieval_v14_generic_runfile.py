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
v14: implement ini file for user inputs
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

import stctm.stellar_retrieval_utilities as sru
import matplotlib.gridspec as gridspec


## The main code starts here

def main(argv):
    '''
    Example:
 	python stellar_retrieval_v14_generic_runfile.py template_ini_stctm.ini
    '''

    if len(argv) > 1:
        iniFile = argv[1]
    else:
        iniFile = 'template_ini_stctm.ini'

    if not os.path.exists(iniFile):
        raise FileNotFoundError('USER ERROR: iniFile does not exist.')

    params = sru.parse_all_input_iniFile(iniFile)


    ## Stellar params
    # Get default log g for the heterogeneity and photosphere
    sru.get_default_logg(params)

    # grid of stellar models
    wv_template_thisR, wv_template_edges_thisR, Teffs_grid, loggs_grid, models_grid_fixedR = sru.get_stellar_model_grid(params)

    ## MCMC params

    # where the burn-in starts in the chains
    burnin = int(params["nsteps"] * params["frac_burnin"])

    # gaussian priors
    gaussparanames, hyperp_gausspriors = sru.parse_gaussian_inputs(params) # mean and std of the Gaussian prior on Tphot

    # specify whether we want to fit fspot/ffac with prior uniform in log or lin space
    fitLogfSpotFac = sru.parse_range_string(params["fitLogfSpotFac"], as_type=int)
    hyperp_logpriors = sru.parse_range_string(params["hyperp_logpriors"])


    ## Load observed spectrum
    spec = ptspec.TransSpec(params["path_to_spec"], inputtype=params["spec_format"])
    # pdb.set_trace()
    ## Prepare for MCMC
    if 1:

        print("Initializing params for MCMC...")

        param, fitparanames = sru.init_default_and_fitted_param_fromDict(params,
                                                            Dscale_guess=np.median(spec["yval"]))

        print("\n** Fitparanames:", fitparanames)
        print("\n** Param: ", param)

        # save the fixed R spectra as baseline specs to do the integration after the fact
        fl_targetR_phot_spot_fac_baseline = sru.get_closest_models_from_grid(param, models_grid_fixedR, Teffs_grid, loggs_grid)

        runname = "fit"
        for p in fitparanames:
            runname = runname + "_"+p
        runname = runname + "_"+ params["res_suffix"]
        results_folder = "../../stellar_contamination_results/"+runname+"/"
        print("Results folder:", results_folder)

        # Setup for post-processing

        if not os.path.isdir(results_folder):
           os.makedirs(results_folder)

        this_dir = os.getcwd() + "/"
        res_dir = os.sep.join(__file__.split(os.sep)[:-2]) + "/../stellar_contamination_results/" + runname + "/"
        utils_script = str(sru.__file__)

        this_script = __file__.split(os.sep)[-1]

        sru.save_ref_files(this_dir, this_script, iniFile, utils_script, res_dir)

        # define parameters for emcee run
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



    ## plot observed spectrum + run emcee
    if 1:
        print("Plotting fitted spectrum...")

        # plot spectrum
        fig, ax = spec.plot()
        fig.savefig(results_folder + "stctm_fitted_spectrum.pdf")

        print("Running MCMC...")

        sampler.run_mcmc(pos, params["nsteps"], progress=True, store=True)


        print("\nRecall -- MCMC setup:")
        print("\n** Fitparanames:", fitparanames)
        print("\n** Param: ", param)

    ## Post-process + Get blobs
    ## ----- Post-processing ----
    if 1: # chainplot


        fig, axes = sru.chainplot(sampler.chain[:, burnin:, :], labels=fitparanames)
        if params["save_fit"]:
            fig.savefig(results_folder + "stctm_chainplot_noburnin_"+runname+'.png')

        fig, axes = sru.chainplot(sampler.chain[:, :, :], labels=fitparanames)
        if params["save_fit"]:
            fig.savefig(results_folder + "stctm_chainplot_"+runname+'.png')

    ##
    if 1: #saving

        rs = sru.save_mcmc_to_pandas(results_folder, runname, sampler, burnin, ndim,
                                     fitparanames, params["save_fit"])
        bestfit, ind_bestfit, ind_maxprob, parabestfit, samples, t_res = rs


    ## Setup for plotting

    sample_spectra = None # for the first time running through plot_stctm_samples_res()

    ## plot 1,2,3 sigma from fixed resolution contamination spectra
    if 1:
        if sample_spectra is None:
            spec.meta["color"] = "k"

            fig, ax, sample_spectra = sru.plot_stctm_samples_res(spec, param, fitparanames,
                                      ind_bestfit, samples, Teffs_grid, loggs_grid,
                                      wv_template_thisR,
                                      models_grid_fixedR,
                                      sample_spectra=None, modelgrid_resP=10000,
                                      target_resP=params["target_resP"],N_samp=1000, ax=None,
                                      bestfit_color = 'k', color="coral",plot3sig=True,
                                      plot2sig=True, plot1sig=True, plotmedian=True,
                                      plotbestfit=True, legend_loc=4, save_csv=True,
                                      results_folder=results_folder, runname=runname)


            ax.set_xlim(np.min(spec["waveMin"])-params["pad"]/2, np.max(spec["waveMax"])+params["pad"])
            ax.set_ylim(0.8*np.median(spec['yval']), 1.15*np.median(spec['yval']))

            if params["save_fit"]:
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_1_2_3sigma_"+runname+'.png')
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_1_2_3sigma_"+runname+'.pdf')
            sru.xspeclog(ax,level=1)
            if params["save_fit"]:
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_logwave_1_2_3sigma_"+runname+'.png')
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_logwave_1_2_3sigma_"+runname+'.pdf')


        else:
            # when sample spectra were already calculated:
            spec.meta["color"] = "k"
            fig, ax, _ = sru.plot_stctm_samples_res(spec, param, fitparanames,
                                      ind_bestfit, samples, Teffs_grid, loggs_grid,
                                      wv_template_thisR,
                                      models_grid_fixedR,
                                      sample_spectra=sample_spectra, modelgrid_resP=10000,
                                      target_resP=params["target_resP"],N_samp=1000, ax=None,
                                      bestfit_color = 'k', color="coral",plot3sig=True,
                                      plot2sig=True, plot1sig=True, plotmedian=True,
                                      plotbestfit=True, legend_loc=4, save_csv=True)
            ax.set_xlim(np.min(spec["waveMin"])-params["pad"]/2, np.max(spec["waveMax"])+params["pad"])

            if params["save_fit"]:
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_1_2_3sigma_"+runname+'.png')
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_1_2_3sigma_"+runname+'.pdf')
            sru.xspeclog(ax,level=1)
            if params["save_fit"]:
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_logwave_1_2_3sigma_"+runname+'.png')
                fig.savefig(results_folder + "stctm_resP"+str(params["target_resP"])+"_logwave_1_2_3sigma_"+runname+'.pdf')



    ## Corner plot
    if 1: #corner plot
        fig = sru.plot_custom_corner(samples, fitparanames, parabestfit)

        suffix = "_custom"
        if params["save_fit"]:
            fig.savefig(results_folder+"stctm_corner_bestfit_"+runname+suffix+".pdf")

    ## Save blobs + plot best fit
    if 1:
        # Blobs
        blobs = sampler.get_blobs()
        flat_st_ctm_models = blobs.T[:, burnin:]["st_ctm_model"].reshape((-1))
        if params["save_fit"]:
            np.save(results_folder+"st_ctm_model_blobs_"+runname+".npy", flat_st_ctm_models)


        fig, ax = sru.plot_maxlike_and_maxprob(spec, param, parabestfit, ind_maxprob,
                                               ind_bestfit, fitparanames, flat_st_ctm_models, pad=params["pad"])

        if params["save_fit"]:
            fig.savefig(results_folder+"stctm_bestfit_model_with_obs_"+runname+".png")

        sru.save_bestfit_stats(spec, ind_bestfit, fitparanames, flat_st_ctm_models,
                               results_folder, runname, save_fit=params["save_fit"])

        t_defaultparam = table.Table([param])
        if params["save_fit"]:
            aio.ascii.write(t_defaultparam, results_folder + "stctm_defaultparams_" + runname + '.csv', format='csv',
                            overwrite=True)

        st_ctm_models = np.array([flat_st_ctm_models[i] for i in range(flat_st_ctm_models.size)])

    ## plot 1,2,3 sigma with blobs
    if 1:
        fig, ax = sru.plot_stctm_blobs(spec, st_ctm_models,
                                  ind_bestfit,
                                  bestfit_color = 'k', color="coral",
                                  plot3sig=True,plot2sig=True, plot1sig=True, plotmedian=True,
                                  plotbestfit=True, legend_loc=4, save_csv=True,
                                  results_folder=results_folder, runname=runname)

        ax.set_xlim(np.min(spec["waveMin"])-params["pad"]/2, np.max(spec["waveMax"])+params["pad"])
        ax.set_ylim(0.8*np.median(spec['yval']), 1.15*np.median(spec['yval']))

        if params["save_fit"]:
            fig.savefig(results_folder+"stctm_1_2_3_sigma_"+runname+".pdf")

        if params["save_fit"]:
            fig.savefig(results_folder+"stctm_1_2_3_sigma_"+runname+".png")


    ## plot 1,2,3 sigma with the amplitude of stellar contamination
    if 1: # v1: for the paper

        fig = plt.figure()
        gs = gridspec.GridSpec(3, 1)
        ax = fig.add_subplot(gs[:2, :])

        sru.plot_stctm_blobs(spec, st_ctm_models,
                                  ind_bestfit,ax=ax,
                                  bestfit_color = 'k', color="coral",
                                  plot2sig=True,plot3sig=True, plot1sig=True, plotmedian=True,
                                  plotbestfit=True, legend_loc=4)
        ax.set_ylabel(r'Transit Depth [ppm]')
        sru.xspeclog(ax,level=1)
        ax.axhline(np.median(spec["yval"]), ls="--", color="k", label="Flat spectrum")
        ax.legend(ncol=2)

        ax2 = fig.add_subplot(gs[2, :], sharex=ax)
        sru.plot_stctm_amplitude(spec, st_ctm_models,
                                  ax=ax2,color="coral")
        ax2.axhline(100., color="grey", ls=":", label="100 ppm")
        sru.xspeclog(ax2,level=1)

        ax2.legend()
        fig.tight_layout()

        ax2.set_ylim(ymin=0)
        if params["save_fit"]:
            fig.savefig(results_folder+"stctm_1_2_3_sigma_withamplitude_noprint"+runname+".pdf")


##
if __name__ == "__main__":
    main(sys.argv)


































