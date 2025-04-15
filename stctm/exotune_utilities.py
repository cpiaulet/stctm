# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:34:45 2023

@author: caroline

Stellar model fitting utility functions
"""

## ------------------- Import modules ------------------#
import os

import numpy as np
import matplotlib.pyplot as plt
import os

os.environ["OMP_NUM_THREADS"] = "1"
import astropy.constants as const
import pdb
import emcee
import astropy.table as table
from copy import deepcopy
import astropy.io as aio

import pandas as pd
import matplotlib.ticker as ticker

from astropy.convolution import convolve,  Gaussian1DKernel
import corner



## Pre-processing

def get_median_spectrum_and_unc(wave_um, spec_ts, flux_calibrated=False):
    """
    Compute the median stellar spectrum and associated uncertainties.

    This function calculates the median spectrum across a time series of spectra
    and estimates the uncertainty in each wavelength bin using the standard deviation
    across the time series. Optionally, the result can be converted to physical flux
    units if the input is flux-calibrated.

    Parameters
    ----------
    wave_um : array_like
        Wavelength array in microns.
    spec_ts : array_like
        2D array of spectral time series with shape (n_times, n_wavelengths).
    flux_calibrated : bool, optional
        If True, applies flux calibration to the output spectra. Default is False.

    Returns
    -------
    median_spectrum : ndarray
        The median stellar spectrum across the time series.
    err_median_spectrum : ndarray
        The estimated uncertainties in each wavelength bin.
    """

    # calculate median stellar spectrum
    print("Calculating median spectrum...")
    median_spectrum = np.nanmedian(spec_ts, axis=0)

    # calculate uncertainties on median stellar spectrum from ptp scatter in each wavelength band
    print("Calculating uncertainties in each wavelength bin from point to point scatter...")
    err_median_spectrum = np.nanstd(spec_ts, axis=0)
    if flux_calibrated:
        median_spectrum = median_spectrum * 1e-6 * 1e-23 * const.c.value * 1e+6 * 1e-4 / wave_um ** 2 * 1e4
        err_median_spectrum = err_median_spectrum * 1e-6 * 1e-23 * const.c.value * 1e+6 * 1e-4 / wave_um ** 2 * 1e4

    return median_spectrum, err_median_spectrum


def plot_median_spectrum_and_unc(wave_arr, median_spectrum, err_median_spectrum,
                                 label="Median from entire time series"):
    """
    Plot median stellar spectrum with uncertainty shading.

    Parameters
    ----------
    wave_arr : array_like
        Wavelength array in microns.
    median_spectrum : array_like
        Median spectrum values corresponding to each wavelength.
    err_median_spectrum : array_like
        Uncertainties associated with each value in the median spectrum.
    label : str, optional
        Label for the spectrum in the plot legend. Default is
        "Median from entire time series".

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    ax : matplotlib.axes._subplots.AxesSubplot
        The matplotlib axes object containing the plot.
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.plot(wave_arr, median_spectrum, color="k", zorder=0, label=label)
    ax.fill_between(wave_arr, median_spectrum - err_median_spectrum, median_spectrum + err_median_spectrum,
                    color="gray", zorder=-9)

    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel("Abs. calibrated flux [mJy]")
    ax.legend(loc=1)

    return fig, ax

def get_scaling_guess(param, T_grid, logg_grid, model_wv, models_grid, spec,
                      overwrite_param=True, wave_min_match_um=1.0, wave_max_match_um=2.0):
    """
    Estimate an initial flux scaling factor (Fscale) between an observed spectrum and a model spectrum
    by comparing their median flux values over a specified wavelength range.

    Parameters
    ----------
    param : dict
        Dictionary of parameters (e.g., for model fitting). May be updated with 'logFscale' and 'Fscale'.
    T_grid : array-like
        Grid of effective temperatures used for selecting the closest model.
    logg_grid : array-like
        Grid of surface gravities used for selecting the closest model.
    model_wv : array-like
        Wavelength grid corresponding to the model spectra.
    models_grid : ndarray
        Grid of model spectra to compare against the observed data.
    spec : object
        Observed spectrum object. Must have attributes: `wave`, `waveMin`, `waveMax`, `yval`, and `yerrUpp`.
    overwrite_param : bool, optional
        If True, updates `param` with the estimated Fscale value. Default is True.
    wave_min_match_um : float, optional
        Lower bound (in microns) of the wavelength range used for matching spectra. Default is 1.0.
    wave_max_match_um : float, optional
        Upper bound (in microns) of the wavelength range used for matching spectra. Default is 2.0.

    Returns
    -------
    Fscale_guess : float
        Estimated flux scaling factor based on the median ratio of observed to model fluxes.
    model_int : array-like
        Model spectrum integrated over the observation bins.
    """

    waveMin = np.array(spec.waveMin)
    # pdb.set_trace()
    waveMax = np.array(spec.waveMax)
    yval = np.array(spec.yval)
    yerrUpp = np.array(spec.yerrUpp)
    spec_pickleformat = (waveMin, waveMax, yval, yerrUpp)

    model_phot_fixR, _, _ = get_closest_models_from_grid(param, models_grid, T_grid, logg_grid)

    # identify wavelengths in common
    ind_obs = np.where((spec.wave > wave_min_match_um) * (spec.wave < wave_max_match_um))

    model_int = calc_spectrum_model_and_integrate(1, 0, 0,
                                                  spec_pickleformat, model_wv, model_phot_fixR, model_phot_fixR,
                                                  model_phot_fixR)

    median_model = np.nanmedian(model_int[ind_obs])
    median_obs = np.nanmedian(spec.yval[ind_obs])
    Fscale_guess = median_obs / median_model

    print("\nGuessing Fscale=", Fscale_guess)
    if overwrite_param:
        param["logFscale"] = np.log10(Fscale_guess)
        param["Fscale"] = Fscale_guess
    return Fscale_guess, model_int

## Modeling

def Fscale_obs(Fscale, fspot, ffac, model_phot_fixR, model_spot_fixR, model_fac_fixR):
    """
    Compute the model spectrum scaled by a flux factor.

    Parameters
    ----------
    Fscale : float
        Flux scaling factor to apply to the combined model spectrum.
    fspot : float
        Spot covering fraction on the stellar surface.
    ffac : float
        Facula covering fraction on the stellar surface.
    model_phot_fixR : array_like
        Flux values of the photosphere model at each wavelength.
    model_spot_fixR : array_like
        Flux values of the spot model at each wavelength.
    model_fac_fixR : array_like
        Flux values of the facula model at each wavelength.

    Returns
    -------
    ndarray
        Scaled full stellar model spectrum.
    """

    phot = (1 - ffac - fspot) * model_phot_fixR
    full_model = phot

    if fspot > 0:
        spot = fspot * model_spot_fixR
        full_model = full_model + spot
    if ffac > 0:
        fac = ffac * model_fac_fixR
        full_model = full_model + fac
    return Fscale * full_model

def integ_model(lam_min,lam_max,wave,F):
    """
    Integrate a model over specified wavelength ranges.

    Computes the average flux within multiple wavelength intervals by integrating
    the flux density and normalizing by the width of each interval.

    Parameters
    ----------
    lam_min : array_like
        Array of lower bounds for the wavelength intervals.
    lam_max : array_like
        Array of upper bounds for the wavelength intervals.
    wave : array_like
        Wavelength array corresponding to the model flux values.
    F : array_like
        Model flux values corresponding to each wavelength in `wave`.

    Returns
    -------
    F_int : ndarray
        Array of average integrated fluxes over the specified wavelength ranges.
    """

    F_int = []
    dwave=wave[1:]-wave[:-1]

    for i in range(lam_min.size):
        ind = np.where((wave>=lam_min[i])*(wave<lam_max[i]))
        numerator = np.trapz(F[ind] * dwave[ind],x=wave[ind])
        denominator = np.trapz(dwave[ind],x=wave[ind])
        F_int.append(numerator/denominator)
    return np.array(F_int)


def calc_spectrum_model_and_integrate(Fscale, fspot, ffac, spec, model_wv, model_phot_fixR, model_spot_fixR,
                                      model_fac_fixR):
    """
    Compute a stellar model and integrate it over spectral bandpasses.

    Parameters
    ----------
    Fscale : float
        Flux scaling factor to apply to the integrated model.
    fspot : float
        Spot covering fraction on the stellar surface.
    ffac : float
        Facula covering fraction on the stellar surface.
    spec : tuple
        Tuple containing:
            - waveMin (array_like): Lower bounds of wavelength bins.
            - waveMax (array_like): Upper bounds of wavelength bins.
            - yval (array_like): Observed flux values (not used here).
            - yerrLow (array_like): Lower error values (not used here).
    model_wv : array_like
        Wavelength array for the model spectra.
    model_phot_fixR : array_like
        Flux values of the photosphere model at each wavelength.
    model_spot_fixR : array_like
        Flux values of the spot model at each wavelength.
    model_fac_fixR : array_like
        Flux values of the facula model at each wavelength.

    Returns
    -------
    ndarray
        Flux-calibrated model integrated over the specified wavelength bins.
    """

    waveMin, waveMax, yval, yerrLow = spec
    # calculate stellar contamination model

    phot = (1 - ffac - fspot) * model_phot_fixR
    full_model = phot

    if fspot > 0:
        spot = fspot * model_spot_fixR
        full_model = full_model + spot
    if ffac > 0:
        fac = ffac * model_fac_fixR
        full_model = full_model + fac

    # integrate in bandpass
    model_int = integ_model(waveMin, waveMax, model_wv, full_model)

    return Fscale * model_int


def get_closest_models_from_grid(param, models_grid_fixedR, Teffs_grid, loggs_grid):
    """
    Retrieve the closest-matching model spectra for the photosphere, spot, and faculae
    components based on input stellar parameters.

    Parameters
    ----------
    param : dict
        Dictionary containing stellar parameters. Must include:
    models_grid_fixedR : ndarray
        Grid of model spectra, indexed by [Teff index, logg index], assumed to be at fixed radius.
    Teffs_grid : array-like
        Grid of effective temperatures corresponding to the first dimension of the model grid.
    loggs_grid : array-like
        Grid of surface gravities corresponding to the second dimension of the model grid.

    Returns
    -------
    model_phot_fixR : array-like
        Model spectrum corresponding to the closest match for the photosphere.
    model_spot_fixR : array-like
        Model spectrum corresponding to the closest match for the spot.
    model_fac_fixR : array-like
        Model spectrum corresponding to the closest match for the faculae.

    """

    # Spot model
    ind_Tspot = np.argmin(np.abs(Teffs_grid - param["Tspot"]))
    ind_logg_het = np.argmin(np.abs(loggs_grid - param["logghet"]))
    model_spot_fixR = models_grid_fixedR[ind_Tspot, ind_logg_het]

    # Faculae model
    ind_Tfac = np.argmin(np.abs(Teffs_grid - param["Tfac"]))
    ind_logg_het = np.argmin(np.abs(loggs_grid - param["logghet"]))
    model_fac_fixR = models_grid_fixedR[ind_Tfac, ind_logg_het]

    # Photosphere model
    ind_Tphot = np.argmin(np.abs(Teffs_grid - param["Tphot"]))
    ind_logg_phot = np.argmin(np.abs(loggs_grid - param["loggphot"]))
    model_phot_fixR = models_grid_fixedR[ind_Tphot, ind_logg_phot]

    return model_phot_fixR, model_spot_fixR, model_fac_fixR

## MCMC setup

def init_default_and_fitted_param(Tphot, met, logg_phot, fitspot=True, fitfac=True,
                                  fitThet=False, fitTphot=False, fitFscale=True,
                                  fitlogg_phot=False, fitdlogg_het=False,
                                  fitLogfSpotFac=[0, 0],
                                  fiterrInfl=False,
                                  Fscale_guess=1, dlogg_het_guess=0.0):
    """
    Initialize default stellar parameters and define which are to be fitted.

    Parameters
    ----------
    Tphot : float
        Photosphere temperature in Kelvin.
    met : float
        Logarithmic stellar metallicity [dex].
    logg_phot : float
        Logarithmic surface gravity of the stellar photosphere.
    fitspot : bool, optional
        If True, include starspots in the model. Default is True.
    fitfac : bool, optional
        If True, include faculae in the model. Default is True.
    fitThet : bool, optional
        If True, fit the temperature differences of spots and faculae. Default is False.
    fitTphot : bool, optional
        If True, fit the photospheric temperature. Default is False.
    fitFscale : bool, optional
        If True, fit the flux scaling factor (logFscale). Default is True.
    fitlogg_phot : bool, optional
        If True, fit the photospheric surface gravity (logg_phot). Default is False.
    fitdlogg_het : bool, optional
        If True, fit the difference in log(g) between photosphere and heterogeneities. Default is False.
    fitLogfSpotFac : list of bool, optional
        Two-element list indicating whether to fit log10 of spot and facula covering fractions,
        respectively. Default is [0, 0].
    fiterrInfl : bool, optional
        If True, fit an error inflation parameter (logErrInfl). Default is False.
    Fscale_guess : float, optional
        Initial guess for the flux scaling factor. Default is 1.
    dlogg_het_guess : float, optional
        Initial guess for the log(g) difference of heterogeneities. Default is 0.0.

    Returns
    -------
    param : dict
        Dictionary containing initialized model parameters and their default values.
    fitparanames : list of str
        List of parameter names to be varied during fitting.
    """

    param = dict()
    fitparanames = []
    param["Tphot"] = Tphot
    param["met"] = met
    param['loggphot'] = logg_phot
    param["dlogghet"] = 0.
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
    param with derived parameters updated
    """
    param["Tspot"] = param["Tphot"] + param["deltaTspot"]
    param["Tfac"] = param["Tphot"] + param["deltaTfac"]
    param["logghet"] = param["loggphot"] + param["dlogghet"]
    param["Fscale"] = 10 ** (param["logFscale"])
    param["errInfl"] = 10 ** (param["logErrInfl"])
    if "log_fspot" in param.keys():
        param["fspot"] = 10 ** param["log_fspot"]
    if "log_ffac" in param.keys():
        param["ffac"] = 10 ** param["log_ffac"]
    return param


def get_param_priors(param, gaussparanames, hyperp_gausspriors=[],
                     fitLogfSpotFac=[False, False], hyperp_logpriors=[], T_grid=[2300, 10000],
                     logg_grid=[2.5, 5.5], Fscale_guess=1.):
    """
    Generate prior bounds for model parameters based on default ranges, user-specified Gaussian priors,
    and optional log-space fitting options for specific parameters.

    Parameters
    ----------
    param : dict
        Dictionary of initial/default parameter values (e.g., Tphot, logg_phot).
    gaussparanames : list of str
        List of parameter names for which Gaussian priors are to be applied.
    hyperp_gausspriors : list of list, optional
        Each entry should be [mean, std] for the corresponding parameter in `gaussparanames`.
        Used to narrow the prior range to within ±5σ around the mean.
    fitLogfSpotFac : list of bool, optional
        Two-element list indicating whether to use log-uniform priors for 'fspot' and 'ffac', respectively.
    hyperp_logpriors : list, optional
        Lower and upper bounds to apply when using log-space priors for 'fspot' or 'ffac'.
    T_grid : list of float
        Temperature grid defining the allowable range for temperature-related parameters (e.g., Tphot).
    logg_grid : list of float, optional
        Surface gravity grid defining the bounds for logg-related parameters.
    Fscale_guess : float, optional
        Initial guess for the flux scaling factor, used to set the bounds for 'logFscale'.

    Returns
    -------
    dict
        Dictionary mapping each parameter name to a [lower_bound, upper_bound] list representing its prior range.
        For parameters fitted in log-space, the key is prefixed with 'log_'.
    """

    defaultpriors = dict()

    defaultpriors['ffac'] = [0., 0.9]
    defaultpriors['fspot'] = [0., 0.9]
    defaultpriors['deltaTfac'] = [100, T_grid[-1] - param["Tphot"]]
    defaultpriors['deltaTspot'] = [T_grid[0] - param["Tphot"], -100.]
    defaultpriors["Tphot"] = [T_grid[0], T_grid[-1]]
    defaultpriors["logFscale"] = [np.log10(Fscale_guess) - 5, np.log10(Fscale_guess) + 5]
    defaultpriors["logErrInfl"] = [0, np.log10(50)]

    defaultpriors["loggphot"] = [np.max([logg_grid[0], 2.5]), np.min([5.5, logg_grid[-1]])]
    defaultpriors["dlogghet"] = [logg_grid[0] - param["loggphot"], 0.]

    parampriors = dict()

    # check parameters for log priors
    allparanames = ['ffac', "fspot", "deltaTfac", "deltaTspot", "Tphot", "logFscale", "loggphot", "dlogghet",
                    "logErrInfl"]
    for par in allparanames:
        if par not in ["fspot", "ffac"]:
            parampriors[par] = defaultpriors[par]
        else:
            if par == "fspot":
                if fitLogfSpotFac[0]:
                    lowlim = hyperp_logpriors[0]
                    upplim = hyperp_logpriors[1]
                    parampriors["log_" + par] = [lowlim, upplim]
                else:
                    parampriors[par] = defaultpriors[par]
            elif par == "ffac":
                if fitLogfSpotFac[1]:
                    lowlim = hyperp_logpriors[0]
                    upplim = hyperp_logpriors[1]
                    parampriors["log_" + par] = [lowlim, upplim]
                else:
                    parampriors[par] = defaultpriors[par]
            else:
                parampriors[par] = defaultpriors[par]

    for par in param:
        if par in gaussparanames:
            ind_gaussparam = np.where(gaussparanames == par)[0][0]
            mean_gaussparam = hyperp_gausspriors[ind_gaussparam][0]
            std_gaussparam = hyperp_gausspriors[ind_gaussparam][1]
            new_prior = [np.max([mean_gaussparam - 5. * std_gaussparam, parampriors[par][0]]),
                         np.min([mean_gaussparam + 5. * std_gaussparam, parampriors[par][1]])]
            parampriors[par] = new_prior
    return parampriors


def lnprior(theta, param, fitparanames, gaussparanames, hyperp_gausspriors,
            fitLogfSpotFac, hyperp_logpriors, T_grid, logg_grid, Fscale_guess):
    """
    Get log prior for a set of model parameters.

    This function evaluates whether the given parameter vector `theta` falls within the
    allowed prior bounds (uniform or Gaussian priors) and applies additional constraints such as
    ensuring the sum of 'fspot' and 'ffac' does not exceed 1.

    Parameters
    ----------
    theta : array-like
        Array of current values for the fitted model parameters
    param : dict
        Dictionary of default parameters used to compute derived parameters.
        This will be updated with values from `theta`.
    fitparanames : list of str
        Names of the parameters being fit.
    gaussparanames : list of str
        Subset of `fitparanames` for which Gaussian priors should be applied.
    hyperp_gausspriors : list of list
        List of [mean, std] values for each parameter in `gaussparanames`, defining the Gaussian priors.
    fitLogfSpotFac : list of bool
        Two-element list indicating whether 'fspot' and 'ffac' are being fit in log-space.
    hyperp_logpriors : list
        Lower and upper bounds for log-space priors on 'fspot' or 'ffac', if applicable.
    T_grid : array-like
        Grid of allowed temperatures for model evaluation; used to constrain 'Tphot', 'deltaTfac', and 'deltaTspot'.
    logg_grid : array-like
        Grid of allowed surface gravities; used to constrain 'loggphot' and 'dlogghet'.
    Fscale_guess : float
        Initial estimate for the flux scaling factor, used to set bounds on 'logFscale'.

    Returns
    -------
    float
        Logarithm of the prior probability (ln prior). Returns `-np.inf` if any parameter falls outside its allowed prior range,
        or if the sum of 'fspot' and 'ffac' exceeds 1.
    """

    lp = 0.

    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta[i]
    param = get_derived_param(param)

    parampriors = get_param_priors(param, gaussparanames, hyperp_gausspriors=hyperp_gausspriors,
                                   fitLogfSpotFac=fitLogfSpotFac, hyperp_logpriors=hyperp_logpriors,
                                   T_grid=T_grid, logg_grid=logg_grid, Fscale_guess=Fscale_guess)
    for i, paraname in enumerate(fitparanames):
        if parampriors[paraname][0] < theta[i] < parampriors[paraname][1]:
            if paraname in gaussparanames:
                ind_gaussparam = np.where(gaussparanames == paraname)[0][0]
                mean_gaussparam = hyperp_gausspriors[ind_gaussparam][0]
                std_gaussparam = hyperp_gausspriors[ind_gaussparam][1]
                lp = lp - 0.5 * ((theta[i] - mean_gaussparam) ** 2. / (std_gaussparam ** 2.))

            else:
                lp = lp + 0.0
        else:
            lp = - np.inf

    if param["fspot"] + param["ffac"] > 1:
        lp = - np.inf

    return lp


def lnlike(theta, param, fitparanames, spec, models_baseline, T_grid, logg_grid, model_wv, models_grid):
    """
    log-likelihood function

    Parameters
    ----------
    theta : array
        Array of fitted parameter values
    param : dict
        Dictionary of all model parameters, including both default (fixed) and fitted para.
    fitparanames : list of str
        Names of the parameters in `theta` that are being optimized.
    spec : tuple
        Observed spectrum as a tuple: (waveMin, waveMax, yval, yerrLow).
    models_baseline : list of arrays
        Precomputed baseline model spectra (photosphere, spot, faculae) at fixed radius.
    T_grid : array-like
        Grid of effective temperatures corresponding to the model grid.
    logg_grid : array-like
        Grid of surface gravities corresponding to the model grid.
    model_wv : array-like
        Wavelength grid of the model spectra.
    models_grid : ndarray
        3D grid of model spectra indexed by temperature and gravity (fixed radius).

    Returns
    -------
    lnlk : float
        Log-likelihood value computed from the squared residuals between observed and model fluxes.
    model_int : array
        Integrated composite model spectrum matched to the observed wavelength bins.


    """

    waveMin, waveMax, yval, yerrLow = spec

    # expand models from the baseline models tuple
    model_phot_baseline, model_spot_baseline, model_fac_baseline = models_baseline

    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta[i]
    param = get_derived_param(param)

    if "deltaTspot" in fitparanames:
        ind_Tspot = np.argmin(np.abs(T_grid - param["Tspot"]))
        ind_logg_het = np.argmin(np.abs(logg_grid - param["logghet"]))
        model_spot_fixR = models_grid[ind_Tspot, ind_logg_het]
    else:
        model_spot_fixR = model_spot_baseline

    if "deltaTfac" in fitparanames:
        ind_Tfac = np.argmin(np.abs(T_grid - param["Tfac"]))
        ind_logg_het = np.argmin(np.abs(logg_grid - param["logghet"]))
        model_fac_fixR = models_grid[ind_Tfac, ind_logg_het]

    else:
        model_fac_fixR = model_fac_baseline

    if "Tphot" in fitparanames:
        ind_Tphot = np.argmin(np.abs(T_grid - param["Tphot"]))
        ind_logg_phot = np.argmin(np.abs(logg_grid - param["loggphot"]))
        model_phot_fixR = models_grid[ind_Tphot, ind_logg_phot]

    else:
        model_phot_fixR = model_phot_baseline

    if "logFscale" in fitparanames:
        Fscale = param["Fscale"]
    else:
        Fscale = 1

    model_int = calc_spectrum_model_and_integrate(Fscale, param['fspot'], param['ffac'], spec,
                                                  model_wv, model_phot_fixR, model_spot_fixR,
                                                  model_fac_fixR)

    yerrLow = yerrLow * param["errInfl"]
    lnlk = -0.5 * (np.sum((yval - model_int) ** 2. / yerrLow ** 2.))

    return lnlk, model_int


def lnprob(theta, spec, param, fitparanames, gaussparanames, hyperp_gausspriors,
           fitLogfSpotFac, hyperp_logpriors, models_baseline_fixR, T_grid, logg_grid,
           model_wv, models_grid_fixR, Fscale_guess):
    """
    Compute ln prob for a given parameter vector.

    Parameters
    ----------
    theta : array-like
        Array of current values for the fitted model parameters
    spec : StellSpec object
        Stellar spectrum to be fitted
    param : dict
        Dictionary of all model parameter values, with defaults; updated with values from `theta`.
    fitparanames : list of str
        List of parameter names being optimized (order must match `theta`).
    gaussparanames : list of str
        Subset of `fitparanames` for which Gaussian priors are applied.
    hyperp_gausspriors : list of list
        List of [mean, std] for each parameter in `gaussparanames` defining the Gaussian priors.
    fitLogfSpotFac : list of bool
        Flags indicating whether to apply log-space priors to 'fspot' and 'ffac', respectively.
    hyperp_logpriors : list
        Log-space prior bounds for 'fspot' or 'ffac', if applicable.
    models_baseline_fixR : array-like
        Precomputed baseline model spectra (e.g., unspotted models) at fixed radius.
    T_grid : array-like
        Temperature grid used to define valid parameter ranges.
    logg_grid : array-like
        Surface gravity grid used to define valid parameter ranges.
    model_wv : array-like
        Wavelength array corresponding to the model spectra.
    models_grid_fixR : array-like
        Full grid of model spectra with fixed radius.
    Fscale_guess : float
        Initial guess for the flux scaling factor, used for prior construction.

    Returns
    -------
    tuple
        - float : Log-posterior probability (log-prior + log-likelihood). Returns `-np.inf` if prior is invalid.
        - array : Model spectrum corresponding to the input parameters. Returned only if prior is valid.
    """

    lp = lnprior(theta, param, fitparanames, gaussparanames, hyperp_gausspriors,
                 fitLogfSpotFac, hyperp_logpriors, T_grid, logg_grid, Fscale_guess)

    if not np.isfinite(lp):
        return -np.inf, np.zeros(len(spec[2])) * np.nan
        # return -np.inf
    else:

        lnlk, model = lnlike(theta, param, fitparanames, spec, models_baseline_fixR, T_grid,
                             logg_grid, model_wv, models_grid_fixR)

        lnprb = lp + lnlk
        return lnprb, model

##  Saving and basic diagnostics

def save_mcmc_to_pandas(results_folder, runname, sampler, burnin, ndim, fitparanames, save_fit):
    # write csv file and astropy table with samples outside of burn-in
    """
    Save MCMC sampling results to disk and extract best-fit parameters.

    This function processes the MCMC sampler output by removing the burn-in phase,
    reshaping and storing the flattened samples as a Pandas DataFrame and an Astropy Table.
    It also computes and saves the best-fit parameters based on the maximum likelihood.

    Parameters
    ----------
    results_folder : str
        Path to the directory where output files should be saved.
    runname : str
        Identifier to append to output filenames for labeling different runs.
    sampler : emcee.EnsembleSampler
        The MCMC sampler object containing the chains and log-probabilities.
    burnin : int
        Number of initial steps (per walker) to discard as burn-in.
    ndim : int
        Number of parameters in the fit.
    fitparanames : list of str
        Names of the fitted parameters.
    save_fit : bool
        If True, saves the flattened samples and best-fit parameters to CSV files.
        -> Writes two CSV files if `save_fit=True`:
        1. Flattened sample chains with log-probabilities (`exotune_pandas_<runname>.csv`)
        2. Best-fit parameter values (`exotune_bestfit_<runname>.csv`)
    Returns
    -------
    bestfit : pandas.DataFrame
        DataFrame containing best-fit values including maximum likelihood estimate.
    ind_bestfit : int
        Index of the best-fit sample (based on median or quantiles).
    ind_maxprob : int
        Index of the sample with the highest posterior probability.
    parabestfit : ndarray
        Array of parameter values corresponding to the maximum likelihood sample.
    samples : pandas.DataFrame
        Flattened MCMC samples after burn-in, with associated log-probabilities.
    t_res : astropy.table.Table
        Astropy table version of the `samples` DataFrame.
    """

    samples = pd.DataFrame(sampler.chain[:, burnin:, :].reshape((-1, ndim)),
                           columns=[x for x in fitparanames])
    lnprobability = pd.Series(sampler.lnprobability[:, burnin:].reshape(-1),
                              name='lnprobability')
    lnlikesave = pd.Series(sampler.lnprobability[:, burnin:].reshape(-1),
                           name='lnlike')
    panda = pd.concat([lnprobability, lnlikesave, samples], axis=1)
    t_res = table.Table.from_pandas(panda)
    if save_fit:
        aio.ascii.write(t_res, results_folder + "exotune_pandas_" + runname + '.csv', format='csv', overwrite=True)

    # save best fit parameters and quantiles
    bestfit, ind_bestfit, ind_maxprob = samp2bestfit(panda)
    print(bestfit)
    if save_fit:
        bestfit.to_csv(results_folder + 'exotune_bestfit_' + runname + '.csv')

    print("\nMax Likelihood:\n")
    print(bestfit["MaxLike"])

    parabestfit = np.array(bestfit["MaxLike"][2:])
    return bestfit, ind_bestfit, ind_maxprob, parabestfit, samples, t_res


def samp2bestfit(samp):
    """
    Compute summary statistics and best-fit values from MCMC samples.

    Calculates key quantiles (median, 1 and 2 sigma) and identifies the samples
    corresponding to the maximum posterior probability and maximum likelihood.

    Parameters
    ----------
    samp : pandas.DataFrame
        DataFrame containing flattened MCMC samples, including at least the columns
        'lnprobability' and 'lnlike', along with one column for each parameter.

    Returns
    -------
    bestfit : pandas.DataFrame
        DataFrame summarizing posterior statistics for each parameter. Includes:
            - 50th percentile (median)
            - 16th and 84th percentiles (for +/-1 sigma)
            - 2.5th and 97.5th percentiles (for +/-2 sigma)
            - Parameter values at maximum likelihood (MaxLike)
            - Parameter values at maximum posterior probability (MaxProb)
            - Explicit +/-1 sigma and +/-2 sigma intervals based on the above quantiles
    ind_bestfit : int
        Index of the sample with the maximum likelihood.
    ind : int
        Index of the sample with the maximum posterior probability.
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

    """
    Calculate the Bayesian Information Criterion (BIC) for model comparison.

    Parameters
    ----------
    chi2 : float
        The chi-squared statistic of the model fit.
    nDataPoints : int
        The number of independent data points used in the model fitting.
    nPara : int
        The number of free parameters in the model.

    Returns
    -------
    float
        The Bayesian Information Criterion value. Lower values indicate a more
        favorable model under the BIC framework.

    References
    ----------
    Adapted from auxbenneke/utilities.py
    Liddle, A.R. (2007). Information criteria for astrophysical model selection.
    Monthly Notices of the Royal Astronomical Society: Letters, 377, L74–L78.
    https://doi.org/10.1111/j.1745-3933.2007.00306.x
    """

    return chi2 + nPara*np.log(nDataPoints)

def save_bestfit_stats(spec, ind_bestfit, fitparanames, flat_oot_spec_models, results_folder, runname, save_fit=True):
    """
    Compute and optionally save summary statistics (red/chi2, BIC) for the best-fit model.

    Parameters
    ----------
    spec : object
        Observed spectrum object, expected to contain attributes:
        `yval` (observed flux values) and `yerrLow` (associated uncertainties).
    ind_bestfit : int
        Index of the best-fit model (typically from maximum likelihood).
    fitparanames : list of str
        Names of the parameters used in the fit, used to compute degrees of freedom.
    flat_oot_spec_models : ndarray
        Array of modeled spectra (post-burn-in), from which the best-fit model is selected.
    results_folder : str
        Directory path where the output file should be saved.
    runname : str
        Label for identifying this run; appended to output filenames.
    save_fit : bool, default=True
        If True, saves the best-fit statistics to a CSV file.

    Returns
    -------
    Astropy Table containing the computed fit statistics.

    Writes a CSV file if `save_fit=True`. Prints diagnostic info to console.
    Output file: 'exotune_bestfit_stats_<runname>.csv'
    """

    print("Saving stats on the best fit...")
    bestfit_stats = dict()
    best_model = flat_oot_spec_models[ind_bestfit]
    nPara = len(fitparanames)
    nDataPoints = len(spec.yval)
    n_dof = nDataPoints - nPara
    bestfit_stats["ind_postburnin"] = ind_bestfit
    bestfit_stats["chi2"] = np.sum((spec.yval - best_model) ** 2. / spec.yerrLow ** 2.)
    bestfit_stats["redchi2"] = bestfit_stats["chi2"] / n_dof
    bestfit_stats["BIC"] = BIC(bestfit_stats["chi2"], nDataPoints, nPara)
    t_bestfit_stats = table.Table([bestfit_stats])

    if save_fit:
        print("Writing to file...")
        aio.ascii.write(t_bestfit_stats, results_folder + "exotune_bestfit_stats_" + runname + '.csv', format='csv',
                        overwrite=True)
    return t_bestfit_stats


def get_exotune_blobs(exotune_models_blobs, percentiles=[0.2, 2.3, 15.9, 84.1, 97.7, 99.8]):
    """
    Compute percentiles of ExoTUNE model blobs across MCMC samples.

    Parameters
    ----------
    exotune_models_blobs : ndarray
        A 3D array of shape (nsamples, ...), where each entry corresponds to
        a model output stored during MCMC sampling.
    percentiles : list of float, optional
        List of percentiles to compute (in the range 0–100). Default corresponds
        to +/- 1,2,3 sigma intervals.

    Returns
    -------
    exotune_percentiles : ndarray
        Array of the same shape as `exotune_models_blobs` without the first dimension,
        but with a new first dimension equal to the number of percentiles.
        For example, shape will be (n_percentiles, ...) corresponding to each computed percentile.
    """

    exotune_percentiles = np.nanpercentile(exotune_models_blobs, percentiles, axis=0)
    return exotune_percentiles



##  Plotting utilities

def xspeclog(ax, xlim=None, level=1, fmt="%2.1f"):
    """
    Configure logarithmic x-axis ticks for spectrum plots.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axis object on which to set the log-scale x-axis ticks.
    xlim : tuple of float, optional
        Optional (xmin, xmax) tuple to explicitly set the x-axis limits.
    level : int or float, default=1
        Preset level for tick density and appearance:
            - 0.5 : Dense major ticks (0.2 to 6, step 0.2) and finer minor ticks (step 0.1)
            - 1   : Standard major ticks at typical wavelengths (e.g., 1.0, 1.5, 2.0, ...)
                    and minor ticks at 0.1 intervals
            - 2   : Integer-based major ticks from 2 to 8, with no minor ticks
            - 3   : Sparse major ticks with minor ticks at selected in-between values
            - 4   : Very sparse major and minor ticks for compact plots
    fmt : str, default="%2.1f"
        Format string for major tick labels.

    Notes
    -----
    This function modifies the axis in-place and is adapted from `auxbenneke/utilities.py`.
    """

    if level == 0.5:
        majorticks = np.arange(0.2, 6, 0.2)
        minorticks = np.arange(0.2, 6, 0.1)
    if level == 1:
        majorticks = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0]
        minorticks = np.arange(0.3, 6, 0.1)
    if level == 2:
        majorticks = np.r_[1.0, 1.5, np.arange(2, 9)]
        minorticks = None
    if level == 3:
        majorticks = np.r_[1.0, 1.5, 2, 3, 4, 6, 8]
        minorticks = np.r_[5, 7, 9.0]
    if level == 4:
        majorticks = np.r_[1, 2, 3, 4, 6, 8]
        minorticks = np.r_[5, 7, 9]

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


def setAxesFontSize(ax, fontsize):
    """
    Set font size for all axis labels and tick labels of a matplotlib axis.

    Updates the font size of the title, x-axis label, y-axis label,
    and all tick labels associated with the given axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axis object whose fonts will be updated.
    fontsize : int or float
        Font size to apply to all axis text elements.

    Notes
    -----
    This function modifies the axis in-place and is adapted from `auxbenneke/utilities.py`.
    """

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)



    
## Post-processing plots

def get_labels_from_fitparanames(fitparanames):
    """
    Generate LaTeX-formatted plot labels for a list of fitted parameter names.

    Parameters
    ----------
    fitparanames : list of str
        List of parameter names as used in the model and MCMC fitting.

    Returns
    -------
    labels : list of str
        List of corresponding LaTeX-formatted labels suitable for plot axes.

    Notes
    -----
    Parameters not explicitly listed are ignored and not included in the output.
    """

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

def format_param_str(param, fitparanames):
    """
    Generate a formatted summary string of stellar and fitted parameters.

    Parameters
    ----------
    param : dict
        Dictionary of default model parameters.
    fitparanames : list of str
        List of parameter names that were varied in the fit.

    Returns
    -------
    str1 : str
        A formatted string summary of the best fit parameters

    """


    param = get_derived_param(param)

    str1 = "Fitted params: "+str(fitparanames)+"\n"
    str1 = str1 + "Stellar params: Tphot="+str(int(param["Tphot"]*10.)/10.)+" met="+str(param["met"]) + "\n"
    str1 = str1 + "logg_phot="+str(param["loggphot"]) + " logg_het="+str(param["logghet"])+ "\n"
    str1 = str1 + "Tspot="+str(int(param["Tspot"]*10.)/10.)+" Tfac="+str(int(param["Tfac"]*10.)/10.)+"\n"
    str1 = str1 + "fspot=" + str(int(param["fspot"]*10000.)/10000.) + " ffac="+str(int(param["ffac"]*10000)/10000.)
    return str1


def chainplot(samples, labels=None, nwalkers=None, fontsize=None, lw=1, stepRange=None):
    """
    Plot parameter chains from MCMC sampling results.

    Parameters
    ----------
    samples : ndarray
        MCMC samples, either 2D (nsamples, nparameters) or 3D (nwalkers, nsteps, nparameters).
    labels : list of str, optional
        List of parameter names to label each subplot's y-axis. Should match the number of parameters.
    nwalkers : int, optional
        Number of walkers used in the MCMC sampling. Required if reshaping 2D samples to 3D.
    fontsize : int or float, optional
        Font size to apply to subplot labels and tick labels.
    lw : float, default=1
        Line width for the chain plots.
    stepRange : tuple of float, optional
        If provided, defines the (start, end) of the x-axis range for each walker chain plot.
        Used to simulate a time/step axis for 3D samples.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object containing the chain plots.
    axes : ndarray of matplotlib.axes.Axes
        2D array of subplot axes, each showing the chain for one parameter.

    Notes
    -----
    This function is adapted from `auxbenneke/utilities.py`.
    """

    if samples.ndim == 2:
        nsamp = samples.shape[0]
        npara = samples.shape[1]
    elif samples.ndim == 3:
        npara = samples.shape[2]

    if nwalkers is not None:
        samples = samples.reshape([nwalkers, int(nsamp / nwalkers), npara])

    nx = int(np.sqrt(npara)) + 1
    ny = int(npara * 1.0 / nx) + 1

    fig, axes = plt.subplots(nx, ny, sharex=True, figsize=(20, 12))
    if axes.ndim == 1:
        axes = np.array([axes])

    if stepRange is not None:
        steps = np.linspace(stepRange[0], stepRange[1], samples.shape[1])

    for i in range(npara):
        setAxesFontSize(axes[int(i / ny), i % ny], fontsize)
        if samples.ndim == 2:
            axes[int(i / ny), i % ny].plot(samples[:, i].T, color="k", alpha=0.4, lw=lw, rasterized=False)
        elif samples.ndim == 3:
            if stepRange is not None:
                axes[int(i / ny), i % ny].plot(steps, samples[:, :, i].T, color="k", alpha=0.4, lw=lw, rasterized=False)
            else:
                axes[int(i / ny), i % ny].plot(samples[:, :, i].T, color="k", alpha=0.4, lw=lw, rasterized=False)

        if labels is not None:
            axes[int(i / ny), i % ny].set_ylabel(labels[i])

    return fig, axes


def plot_exotune_blobs(spec, oot_models_blobs,
                     ind_bestfit, ax=None,
                     bestfit_color='k', color="b",
                     plot3sig=True, plot2sig=True, plot1sig=True, plotmedian=True,
                     plotbestfit=True, legend_loc=1, save_csv=False,
                     results_folder="", runname=""):
    """
    Plot spectra obtained from the retrieval with uncertainty from blobs.

    Parameters
    ----------
    spec: TransSpec object
        (planet atmosphere observations in transmission)
    oot_models_blobs : array-like
        array of model outputs sampled from the posterior distribution.
    ind_bestfit : int
        Index of the max-like model in `stctm_models_blobs`
    ax : matplotlib.axes.Axes, optional
        Matplotlib axis object to plot on. If None, a new axis is created.
    bestfit_color : str, default 'k'
        Color used to plot the best-fit model.
    color : str, default 'b'
        Color used for median, 1, 2, 3 sigma contours.
    plot3sig : bool, default True
        If True, plot the 3-sigma uncertainty region.
    plot2sig : bool, default True
        If True, plot the 2-sigma uncertainty region.
    plot1sig : bool, default True
        If True, plot the 1-sigma uncertainty region.
    plotmedian : bool, default True
        If True, plot the median model prediction.
    plotbestfit : bool, default True
        If True, plot the best-fit model.
    legend_loc : int, default 1
        Location of the legend in the plot (uses matplotlib legend location codes).
    save_csv : bool, default False
        If True, save a CSV file with the quantiles and best-fit predictions.
    results_folder : str, optional
        Directory path where the CSV file should be saved (if `save_csv=True`).
    runname : str, optional
        Base name used for the saved CSV file (if `save_csv=True`).

    Returns
    -------
    ax : matplotlib.axes.Axes
        The axis object with the plot drawn on it.

    """

    # percentiles = [2.5, 16., 50., 84., 97.5]
    percentiles = [0.2, 2.3, 15.9, 50., 84.1, 97.7, 99.8]

    # for each epoch of each planet, get the median, 1 and 2sigma pred time
    oot_percentiles = np.nanpercentile(oot_models_blobs, percentiles, axis=0)

    # calc for best fit (need ind_bestfit)
    oot_best = oot_models_blobs[ind_bestfit, :]
    if ax is None:
        fig, ax = spec.plot()
    else:
        fig = None
        spec.plot(ax=ax)

    if plot1sig:
        ax.fill_between(spec["wave"], oot_percentiles[2], oot_percentiles[4], color=color, alpha=0.5, zorder=-9,
                        label=r'1$\sigma$')

    if plot2sig:
        ax.fill_between(spec["wave"], oot_percentiles[1], oot_percentiles[5], color=color, alpha=0.3, zorder=-10,
                        label=r'2$\sigma$')

    if plot3sig:
        ax.fill_between(spec["wave"], oot_percentiles[0], oot_percentiles[6], color=color, alpha=0.2, zorder=-10,
                        label=r'3$\sigma$')

    if plotmedian:
        ax.plot(spec["wave"], oot_percentiles[3], color=color, zorder=-8, label=r'Median')

    if plotbestfit:
        ax.plot(spec["wave"], oot_best, color=bestfit_color, zorder=-8, label=r'Max. likelihood')

    ax.legend(loc=legend_loc)

    if save_csv:
        dct = dict()
        dct["wave"] = deepcopy(np.array(spec["wave"]))
        dct["bestfit"] = oot_best
        dct["median"] = oot_percentiles[3]
        dct["+1 sigma"] = oot_percentiles[4]
        dct["-1 sigma"] = oot_percentiles[2]
        dct["+2 sigma"] = oot_percentiles[5]
        dct["-2 sigma"] = oot_percentiles[1]
        dct["+3 sigma"] = oot_percentiles[6]
        dct["-3 sigma"] = oot_percentiles[0]

        t = table.Table(data=dct)
        t.write(results_folder + "blobs_1_2_3_sigma" + runname + ".csv", overwrite=True)
    return fig, ax

def plot_exotune_samples_res(spec, param, fitparanames,
                             ind_bestfit, post_burnin_samples, T_grid, logg_grid,
                             modelgrid_wave,
                             grid_models_fixedresP,
                             sample_spectra=None, modelgrid_resP=10000,
                             target_resP=100, N_samp=1000, ax=None,
                             bestfit_color='k', color="b", plot3sig=True,
                             plot2sig=True, plot1sig=True, plotmedian=True,
                             plotbestfit=True, legend_loc=1, save_csv=True,
                             results_folder="", runname=""):
    """
    Plot ExoTUNE model spectra with uncertainty intervals from MCMC samples.

    Visualizes the model spectra generated from posterior samples
    of the fit parameters, including the maximum likelihood model and confidence intervals
    (1, 2, 3 sigma) derived from the sample distribution. It also optionally saves these results
    to a CSV file.

    Parameters
    ----------
    spec : StellSpec object
        Observed spectrum object.
    param : dict
        Dictionary of default model parameters (updated for each sample).
    fitparanames : list of str
        Names of the fitted parameters.
    ind_bestfit : int
        Index of the maximum likelihood sample within `post_burnin_samples`.
    post_burnin_samples : pandas.DataFrame
        Flattened MCMC samples after burn-in.
    T_grid : array-like
        Grid of temperatures used in the model.
    logg_grid : array-like
        Grid of surface gravities used in the model.
    modelgrid_wave : array-like
        Wavelength array for the spectra in the model grid.
    grid_models_fixedresP : ndarray
        2D array of model spectra at fixed resolution, indexed by (temperature, logg).
    sample_spectra : ndarray, optional
        Precomputed model spectra for posterior samples. If None, they are computed.
    modelgrid_resP : float, default=10000
        Resolving power of the model spectra.
    target_resP : float, default=100
        Target resolving power for visualization and saving.
    N_samp : int, default=1000
        Number of posterior samples to draw and compute spectra from.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None, a new figure and axis are created.
    bestfit_color : str, default='k'
        Color for the maximum likelihood spectrum line.
    color : str, default='b'
        Color for uncertainty bands and median spectrum line.
    plot3sig : bool, default=True
        Whether to plot the +/- 3 sigma envelope.
    plot2sig : bool, default=True
        Whether to plot the +/- 2 sigma envelope.
    plot1sig : bool, default=True
        Whether to plot the +/- 1 sigma envelope.
    plotmedian : bool, default=True
        Whether to plot the median spectrum.
    plotbestfit : bool, default=True
        Whether to plot the maximum likelihood spectrum.
    legend_loc : int, default=1
        Location of the plot legend.
    save_csv : bool, default=True
        Whether to save the resulting model spectra (median, percentile intervals) to a CSV file.
    results_folder : str, optional
        Folder path where CSV will be saved if `save_csv=True`.
    runname : str, optional
        Identifier to append to output filenames.

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        The matplotlib figure object. None if `ax` is provided.
    ax : matplotlib.axes.Axes
        The axis object with the plotted spectra.
    sample_spectra : ndarray
        Array of model spectra generated from posterior samples.
    """

    # get the target wavelength array
    kernel = Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)  # *2 because FWHM = 2 standard deviation
    l = int(kernel.shape[0] / 2)
    wv_array = modelgrid_wave[l:-l]
    # ind = np.where((wv_array>spec["waveMin"][0])*(wv_array<spec["waveMax"][-1]))

    if sample_spectra is None:  # unless the grid of sample spectra is provided
        # draw random samples from the post-burnin samples
        if len(post_burnin_samples) < N_samp:
            N_samp = len(post_burnin_samples)
            print("Too few samples, using N_samp=", N_samp)

        ind_samples = np.random.choice(np.arange(len(post_burnin_samples)), size=N_samp, replace=False)

        # calculate array of sample spectra
        sample_spectra = np.zeros((N_samp, wv_array.size))

        for i in range(N_samp):
            if i % 100 == 0:
                print("i=", i + 1, "/", N_samp)
            theta = post_burnin_samples.iloc[ind_samples[i]]
            for j, paraname in enumerate(fitparanames):
                param[paraname] = theta.loc[paraname]
            param = get_derived_param(param)

            ind_Tspot = np.argmin(np.abs(T_grid - param["Tspot"]))
            ind_Tfac = np.argmin(np.abs(T_grid - param["Tfac"]))
            ind_Tphot = np.argmin(np.abs(T_grid - param["Tphot"]))
            ind_logg_phot = np.argmin(np.abs(logg_grid - param["loggphot"]))
            ind_logg_het = np.argmin(np.abs(logg_grid - param["logghet"]))

            model_spot = grid_models_fixedresP[ind_Tspot, ind_logg_het]
            model_phot = grid_models_fixedresP[ind_Tphot, ind_logg_phot]
            model_fac = grid_models_fixedresP[ind_Tfac, ind_logg_het]

            model_oot = Fscale_obs(param["Fscale"], param['fspot'], param['ffac'],
                                   model_phot, model_spot, model_fac)

            # convolve the model to the target resolving power
            kernel = Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)  # *2 because FWHM = 2 standard deviation
            l = int(kernel.shape[0] / 2)
            # x = wave[l:-l]
            model_oot = convolve(model_oot, kernel)[l:-l]
            # pdb.set_trace()
            sample_spectra[i] = deepcopy(model_oot)

    # calculate spectrum for best-fit parameters
    theta = post_burnin_samples.iloc[ind_bestfit]
    for i, paraname in enumerate(fitparanames):
        param[paraname] = theta.loc[paraname]
    param = get_derived_param(param)

    ind_Tspot = np.argmin(np.abs(T_grid - param["Tspot"]))
    ind_Tfac = np.argmin(np.abs(T_grid - param["Tfac"]))
    ind_Tphot = np.argmin(np.abs(T_grid - param["Tphot"]))
    ind_logg_phot = np.argmin(np.abs(logg_grid - param["loggphot"]))
    ind_logg_het = np.argmin(np.abs(logg_grid - param["logghet"]))

    model_spot = grid_models_fixedresP[ind_Tspot, ind_logg_het]
    model_phot = grid_models_fixedresP[ind_Tphot, ind_logg_phot]
    model_fac = grid_models_fixedresP[ind_Tfac, ind_logg_het]

    model_oot = Fscale_obs(param["Fscale"], param['fspot'], param['ffac'],
                           model_phot, model_spot, model_fac)

    # convolve the model to the target resolving power
    kernel = Gaussian1DKernel(modelgrid_resP / target_resP / 2.35)  # *2 because FWHM = 2 standard deviation
    l = int(kernel.shape[0] / 2)
    # x = wave[l:-l]
    oot_spec_best = convolve(model_oot, kernel)[l:-l]

    # calculate percentiles from sample spectra
    perc_spectra = get_exotune_blobs(sample_spectra, percentiles=[0.2, 2.3, 15.9, 50., 84.1, 97.7, 99.8])

    lowest_z = 1000

    if ax is None:
        fig, ax = spec.plot(marker=".", markersize=1, ls="", color="gray", zorder=lowest_z + 3,
                            label="Observed spectrum")

    else:
        fig = None
        spec.plot(ax=ax, marker=".", markersize=1, color="gray", zorder=lowest_z + 3, label="Observed spectrum")

    if plot1sig:
        ax.fill_between(wv_array, perc_spectra[2], perc_spectra[4], color=color, alpha=0.5, zorder=lowest_z + 1,
                        label=r'1$\sigma$')

    if plot2sig:
        ax.fill_between(wv_array, perc_spectra[1], perc_spectra[5], color=color, alpha=0.3, zorder=lowest_z,
                        label=r'2$\sigma$')

    if plot3sig:
        ax.fill_between(wv_array, perc_spectra[0], perc_spectra[6], color=color, alpha=0.2, zorder=lowest_z,
                        label=r'3$\sigma$')

    if plotmedian:
        ax.plot(wv_array, perc_spectra[3], color=color, zorder=lowest_z + 2, label=r'Median')

    if plotbestfit:
        ax.plot(wv_array, oot_spec_best, color=bestfit_color, zorder=lowest_z + 2, label=r'Max. likelihood')

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
        t.write(results_folder + "fixedR_1_2_3_sigma" + runname + ".csv", overwrite=True)

    return fig, ax, sample_spectra


def plot_corner(samples, plotparams, plot_datapoints=False, smooth=1.5,
                        quantiles=[0.16, 0.5, 0.84], title_kwargs={'fontsize':14},
                        hist_kwargs={"linewidth":3}, rg=None, 
                        show_titles=True, levels=(0.393,0.865,0.989), **kwargs):
    """
    Generate a corner plot to visualize posterior distributions from MCMC samples.

    This function wraps the `corner.corner` plotting with customization of styling, smoothing, quantiles, and
    credible intervals.

    Parameters
    ----------
    samples : ndarray
        MCMC samples array of shape (nsamples, nparameters).
    plotparams : dict
        Dictionary containing plot settings. Must include:
            - "labels" : list of str, axis labels for each parameter
            - "hist_color" : str, color for the histograms and contour plots
    plot_datapoints : bool, default=False
        Whether to overplot individual samples in the 2D contour plots.
    smooth : float, default=1.5
        Smoothing parameter passed to `corner.corner` for KDE.
    quantiles : list of float, default=[0.16, 0.5, 0.84]
        Quantiles to display as vertical lines on the 1D histograms.
    title_kwargs : dict, default={'fontsize': 14}
        Keyword arguments for customizing title appearance in each subplot.
    hist_kwargs : dict, default={"linewidth": 3}
        Keyword arguments for customizing the histogram appearance.
    rg : list or tuple, optional
        Plotting range for each parameter.
    show_titles : bool, default=True
        Whether to show titles with median and quantile values for each parameter.
    levels : tuple of float, default=(0.393, 0.865, 0.989)
        Contour levels for 2D histograms, corresponding to 1, 2, and 3 sigma for a 2D Gaussian.
    **kwargs : dict
        Additional keyword arguments passed to `corner.corner`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object containing the corner plot.
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

def plot_custom_corner(samples, fitparanames, parabestfit, param,gaussparanames,
                       hyperp_gausspriors,fitLogfSpotFac,hyperp_logpriors,T_grid,logg_grid):
    """
    Generate a customized corner plot with parameters ordered for clarity and axes ranges informed by priors.

    Parameters
    ----------
    samples : ndarray
        Array of MCMC samples with shape (n_samples, n_parameters), post burn-in.
    fitparanames : list of str
        Names of the fitted parameters in the order they appear in `samples`.
    parabestfit : ndarray
        Best-fit values for the parameters.
    param : dict
        Dictionary containing default parameter values, used for computing prior ranges.
    gaussparanames : list of str
        Names of parameters with Gaussian priors.
    hyperp_gausspriors : list of lists
        Hyperparameters defining the Gaussian priors for relevant parameters.
    fitLogfSpotFac : list of bool
        Flag indicating whether to fit the log of fspot/ffac parameters.
    hyperp_logpriors : dict
        Hyperparameters defining log-uniform priors for specific parameters.
    T_grid : ndarray
        Grid of stellar effective temperatures used to inform priors.
    logg_grid : ndarray
        Grid of stellar surface gravities used to inform priors.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the customized corner plot.
    """



    # reorder samples
    # if "logErrInfl" in fitparanames:
        # ordered_samples = np.zeros((samples.shape[0],samples.shape[1]-1))
    # else:
    ordered_samples = np.zeros((samples.shape[0],samples.shape[1]))
    ordered_fitparanames_all = ["fspot", "log_fspot","deltaTspot", "ffac", "log_ffac", "deltaTfac", "dlogghet", "Tphot", "loggphot", "logFscale","logErrInfl"]
    ordered_fitparanames = []
    ind = np.where(np.array(fitparanames)=="logFscale")
    Fscale_guess =  10**parabestfit[ind[0][0]]
    parampriors = get_param_priors(param, gaussparanames, hyperp_gausspriors=hyperp_gausspriors,
    fitLogfSpotFac=fitLogfSpotFac, hyperp_logpriors=hyperp_logpriors,
    T_grid=T_grid, logg_grid=logg_grid, Fscale_guess=Fscale_guess)

    ordered_rgs = [] # ranges for each parameter
    for p in ordered_fitparanames_all:
        if p in fitparanames:
            ordered_fitparanames.append(p)
            
    ordered_parabestfit = np.zeros_like(parabestfit)
    for i, p in enumerate(ordered_fitparanames):
        ind = np.where(np.array(fitparanames)==p)
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

    plotparams = dict()
    plotparams["hist_color"] = "coral"
    plotparams["labels"] = get_labels_from_fitparanames(ordered_fitparanames)

    print("\nSanity check:")
    print("Shape of ordered_samples:",ordered_samples.shape)
    print("len of labels:",len(plotparams["labels"]))

    print("\nProceed to plotting:")
    fig = plot_corner(ordered_samples, plotparams, smooth=0.8, fill_contours=True,
                      truths = ordered_parabestfit, truth_color="k")
    # fig = plot_corner(ordered_samples[:,ind], plotparams, smooth=0.8, fill_contours=True, 
    #                   truths = ordered_parabestfit[ind], truth_color="k")    
                                    

    fig.set_dpi(50)  
    fig.set_figheight(11)
    fig.set_figwidth(15)   
    
    return fig


def plot_maxlike_and_maxprob(spec, param, parabestfit, ind_maxprob, ind_bestfit, fitparanames, flat_st_ctm_models, pad=0.25):
    """
    Plot the best-fit stellar models for maximum likelihood and maximum posterior samples.


    Parameters
    ----------
    spec : StellSpec object
        Observed spectrum object
    param : dict
        Dictionary of default model parameters to be updated with best-fit values.
    parabestfit : ndarray
        Array of parameter values corresponding to the best-fit model (typically max likelihood).
    ind_maxprob : int
        Index of the sample with the highest posterior probability.
    ind_bestfit : int
        Index of the sample with the maximum likelihood.
    fitparanames : list of str
        Names of the parameters used in the fit; must match order of `parabestfit`.
    flat_st_ctm_models : ndarray
        Array of model predictions corresponding to all samples.
    pad : float, default=0.25
        Padding to apply on the x-axis when setting plot limits.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object with the plotted spectrum and models.
    ax : matplotlib.axes.Axes
        The axis object containing the plot.
    """
    fig, ax = spec.plot()
    param_bestfit = deepcopy(param)
    for i, p in enumerate(fitparanames):
        param_bestfit[p] = parabestfit[i]
    ax.scatter(spec.wave, flat_st_ctm_models[ind_maxprob], label="Max. Probability", color="slateblue", alpha=1.)
    ax.scatter(spec.wave, flat_st_ctm_models[ind_bestfit], label="Max. Likelihood", color="r", alpha=0.5, marker=".", s=50)

    ymin, ymax = ax.get_ylims()
    ax.text(np.max(spec.waveMax)+pad-pad/4, 0.95*ymax, format_param_str(param_bestfit, fitparanames), fontsize=8, color="k",ha="right", va="top")

    ax.set_xlim(np.min(spec.waveMin)-pad/2, np.max(spec.waveMax)+pad)

    ax.set_title("Best-fit model")
    ax.legend(loc=4)
    return fig, ax  