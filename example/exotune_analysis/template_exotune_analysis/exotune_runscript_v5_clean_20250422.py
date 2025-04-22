# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:40:05 2019

@author: caroline

Fit a stellar spectrum using combination of stellar models
v4: with ini file for inputs
v4: implement custom pool for multiprocessing
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
import sys
import astropy.io as aio
from copy import deepcopy
from multiprocessing import Pool
from matplotlib.gridspec import GridSpec


## The main code starts here

def main(argv):
    '''
        Example:
     	python exotune_runscript_v5_clean_20250422.py template_ini_exotune.ini
        '''

    if len(argv) > 1:
        iniFile = argv[1]
    else:
        iniFile = 'template_ini_exotune.ini'

    if not os.path.exists(iniFile):
        raise FileNotFoundError('USER ERROR: iniFile does not exist.')

    params, label_run = xtu.parse_all_input_exotune_iniFile(iniFile)
    print("\nStarting Run:", label_run)

    ## Pre-processing
    if params["start_from_timeseries"]:
        #  Reading in the observed spectrum
        timeUnique, waveUnique, optspec2D = xtu.read_in_2d_spec(params["path_to_stellar_spec_ts"])

        # Remove time intervals from 2D array to calculate median
        fig1, timeUnique, optspec2D = xtu.calc_and_plot_timeMasked_mednorm_lc(timeUnique, optspec2D, params)

        print("\nCalculating median out-of-transit spectrum and use to initialize StellSpec object...")
        spec = xtu.get_StellSpec_from_spec_timeseries(waveUnique,optspec2D, params, flux_calibrated=True)

    else:
        print("\nRead in median out-of-transit spectrum and use to initialize StellSpec object...")
        t = aio.ascii.read(path_to_spec) # read in the spectrum that was provided
        spec = psspec.StellSpec(t, inputtype=spec_format)

    ##  Apply operations to the StellSpec object

    spec, fig2 = xtu.prep_spec_waveMasked(spec, params)

    print("\nObtain values from spec object to save as picklable object for MCMC...")
    # get the values we need from the spec object while keeping it picklable (only numpy arrays)
    waveMin = np.array(spec.waveMin)
    waveMax = np.array(spec.waveMax)
    yval = np.array(spec.yval)
    yerrUpp = np.array(spec.yerrUpp)
    spec_pickleformat = (waveMin, waveMax, yval, yerrUpp )

    ## Stellar params and models
    # Get default log g for the heterogeneity and photosphere
    xtu.get_default_logg(params)

    # grid of stellar models
    wv_template_thisR, wv_template_edges_thisR, Teffs_grid, loggs_grid, models_grid_fixedR = sru.get_stellar_model_grid(params,preprocessed=True)

    ##  Visualize nominal model in grid on top of observations (scaled)

    param, fitparanames = xtu.init_default_and_fitted_param_fromDict(params)
    param = xtu.get_derived_param(param)

    Fscale_guess, fl_phot_spot_fac, fig3 = xtu.calc_and_plot_Fscale_guess(param, Teffs_grid, loggs_grid,  wv_template_thisR, models_grid_fixedR, spec)

    ##  Stop here if you just wanted to update your data selection
    res_suffix = params["res_suffix"]
    if params["optimize_param"]:
        this_dir = os.getcwd()+"/"
        print("\nSaving previously-created figures...")
        if params["start_from_timeseries"]:
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
        param, fitparanames = xtu.init_default_and_fitted_param_fromDict(params,
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
        if params["start_from_timeseries"]:
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


        xtu.save_ref_files(this_dir, this_script, iniFile, utils_script, res_dir)


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


        other_args = (spec_pickleformat, param, fitparanames, params["gaussparanames"], params["hyperp_gausspriors"],
                      params["fitLogfSpotFac"], params["hyperp_logpriors"],
                      fl_phot_spot_fac, Teffs_grid, loggs_grid,
                      wv_template_thisR, models_grid_fixedR, Fscale_guess)


    ##
    if 1:
        print("\nCreating sampler...")

        if params["parallel"] is False:
            print("Running serial version of the MCMC fit!")
            sampler = emcee.EnsembleSampler(nwalkers, ndim, xtu.lnprob,
                                            args=other_args,
                                            blobs_dtype=dtype)
            print("\nRunning MCMC...")

            sampler.run_mcmc(pos, params["nsteps"], progress=True, store=True)
        else:
            print("Running parallel version of the MCMC fit on", params["ncpu"], "CPUs!")


            with xtu.MyPool(processes=params["ncpu"], lnprob=xtu.lnprob, lnprob_args=other_args) as pool:
                sampler = emcee.EnsembleSampler(
                    nwalkers, ndim, None,  # lnprob is handled internally via the pool
                    pool=pool,
                    blobs_dtype=dtype
                )

                print("\nRunning MCMC...")
                sampler.run_mcmc(pos, params["nsteps"], progress=True, store=True)



    ##  Post-process + Get blobs
    if 1:
        # MCMC: where the burn-in starts in the chains
        burnin = int(params["nsteps"] * params["frac_burnin"])

        print("\nRecall -- MCMC setup:")
        print("\n** Fitparanames:", fitparanames)
        print("\n** Param: ", param)

    ##  ----- Post-processing plots ----
    if 1: # chainplot

        print("\nGenerating chain plot...")

        fig, axes = xtu.chainplot(sampler.chain[:, burnin:, :], labels=fitparanames)
        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_chainplot_noburnin_"+runname+'.png')
            plt.close(fig)

        fig, axes = xtu.chainplot(sampler.chain[:, :, :], labels=fitparanames)
        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_chainplot_"+runname+'.png')
            plt.close(fig)

    ##
    if 1: #saving
        print("\nSaving results to pandas DataFrame...")

        rs = xtu.save_mcmc_to_pandas(res_dir, runname, sampler, burnin, ndim,
                                     fitparanames, params["save_fit"])
        bestfit, ind_bestfit, ind_maxprob, parabestfit, samples, t_res = rs

    ##
    ##  Setup for plotting

    pad = params["pad"]
    target_resP = params["target_resP"]
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

        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.png')
            fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.pdf')
        xtu.xspeclog(ax,level=1)
        if params["save_fit"]:
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

        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.png')
            fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.pdf')
        xtu.xspeclog(ax,level=1)
        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.png')
            fig.savefig(res_dir + "exotune_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.pdf')
            plt.close(fig)

    ##  Save blobs + plot best fit
    if 1:
        print("\nSaving blobs...")
        # Blobs
        blobs = sampler.get_blobs()
        if params["ncpu"]>1:
            # Close the thread pool
            pool.close()
            pool.join()
        flat_oot_spec_models = blobs.T[:, burnin:]["oot_spec_model"].reshape((-1))
        if params["save_fit"]:
            np.save(res_dir+"oot_spec_model_blobs_"+runname+".npy", flat_oot_spec_models)

        print("\nPlotting best fit model with observations...")
        fig, ax = xtu.plot_maxlike_and_maxprob(spec, param, parabestfit, ind_maxprob,
                                               ind_bestfit, fitparanames, flat_oot_spec_models, pad=pad)

        if params["save_fit"]:
            fig.savefig(res_dir+"exotune_bestfit_model_with_obs_"+runname+".png")
            plt.close(fig)

        print("\nComputing statistics (chi-squared, BIC) on the run results...")
        xtu.save_bestfit_stats(spec, ind_bestfit, fitparanames, flat_oot_spec_models,
                               res_dir, runname, save_fit=save_fit)


        print("\nSaving default parameters to file...")
        t_defaultparam = table.Table([param])
        if params["save_fit"]:
            aio.ascii.write(t_defaultparam, res_dir+"exotune_defaultparams_"+runname+'.csv', format='csv', overwrite=True)

        oot_spec_models = np.array([flat_oot_spec_models[i] for i in range(flat_oot_spec_models.size)])

    ## plot 1,2,3 sigma with blobs

    if 1:
        fig, ax = xtu.plot_exotune_blobs(spec, oot_spec_models,
                                       ind_bestfit,
                                       bestfit_color='k', color="coral",
                                       plot3sig=True, plot2sig=True, plot1sig=True, plotmedian=True,
                                       plotbestfit=True, legend_loc=4, save_csv=True,
                                       results_folder=res_dir, runname=runname)

        ax.set_xlim(np.min(spec.waveMin) - pad / 2, np.max(spec.waveMax) + pad)
        ax.set_ylim(0.8 * np.median(spec.yval), 1.15 * np.median(spec.yval))

        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_1_2_3_sigma_" + runname + ".pdf")

        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_1_2_3_sigma_" + runname + ".png")

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
        if params["ncpu"]>1:
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

        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.png')
            fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_1_2_3sigma_"+runname+'.pdf')
        xtu.xspeclog(axspec,level=1)
        if params["save_fit"]:
            fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.png')
            fig.savefig(res_dir + "exotune_combo_resP"+str(target_resP)+"_logwave_1_2_3sigma_"+runname+'.pdf')
            plt.close(fig)

    ##  Corner plot
    # try: #corner plot
    print("\nCreating corner plot...")
    fig = xtu.plot_custom_corner(samples, fitparanames, parabestfit, param,
                       gaussparanames,hyperp_gausspriors,fitLogfSpotFac,hyperp_logpriors,Teffs_grid,loggs_grid)


    suffix = "_custom"
    if params["save_fit"]:
        fig.savefig(res_dir+"exotune_corner_bestfit_"+runname+suffix+".pdf")
        plt.close(fig)
    print("\nPost-processing done!")

##
if __name__ == "__main__":
    main(sys.argv)


