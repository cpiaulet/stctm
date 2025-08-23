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
import stctm.pystellspec as psspec

import h5py
import scipy.signal as sig

import pandas as pd
import matplotlib.ticker as ticker

from astropy.convolution import convolve,  Gaussian1DKernel
import corner
import multiprocessing as mp
import shutil

import configparser
import argparse
from . import stellar_retrieval_utilities as sru
try:
    import tomllib as tomli  # Py 3.11+
except Exception:
    import tomli             # backwards-compatible

## Parsing

def convert_to_bool(x):
    """
    Convert a value into a boolean.

    Parameters
    ----------
    x : Any
        Value to convert. Can be bool, str, int, or None.

    Returns
    -------
    bool
        Converted boolean value.
    """
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    s = str(x).strip().lower()
    if s in ("true", "t", "1", "yes", "y"):
        return True
    if s in ("false", "f", "0", "no", "n"):
        return False
    return bool(x)


def parse_range_pair(value, as_type=float):
    """
    Normalize a 2-value range into a list [min, max].

    Parameters
    ----------
    value : str, list, tuple, or None
        The range to parse. Can be "a_b", [a, b], (a, b), or None/"".

    as_type : type, optional
        The type to cast values to (default: float).

    Returns
    -------
    list
        A 2-element list [min, max] of the given type, or [] if empty.
    """
    if value is None or value == "":
        return []
    if isinstance(value, (list, tuple)):
        if len(value) != 2:
            raise ValueError(f"Range must have 2 items, got {value}")
        return [as_type(value[0]), as_type(value[1])]
    parts = str(value).split("_")
    if len(parts) != 2:
        raise ValueError(f"Could not parse range '{value}'. Expected 'min_max' or [min,max].")
    return [as_type(parts[0]), as_type(parts[1])]

# table aliases defined for transition from .ini file to .toml file as input
_TABLE_ALIASES = {
    "choice of inputs": ("choice of inputs", "choice_of_inputs"),
    "preprocessing": ("preprocessing",),
    "saving options": ("saving options", "saving_options"),
    "stellar params": ("stellar params", "stellar_params"),
    "stellar models": ("stellar models", "stellar_models"),
    "MCMC params": ("MCMC params", "MCMC_params"),
    "priors": ("priors",),
    "plotting": ("plotting",),
}

def parse_mask_ranges(v):
    """
    Normalize a mask specification into a list of ranges.

    Parameters
    ----------
    v : str, list, or None
        Mask specification. Can be:
        - "" or None -> []
        - "a_b|c_d"  -> [[a, b], [c, d]]
        - [[a, b], [c, d]] -> returned as list of floats

    Returns
    -------
    list
        A list of [min, max] float pairs.
    """
    if v in (None, ""):
        return []
    if isinstance(v, str):
        out = []
        for pair in v.split("|"):
            pair = pair.strip()
            if not pair:
                continue
            out.append(parse_range_pair(pair, float))
        return out
    return [[float(p[0]), float(p[1])] for p in v]


def preprocess_gaussians(params):
    """
    Pre-normalize Gaussian prior fields so the sru.parse_gaussian_inputs() forwarder
    receives the shapes it expects.

    Parameters
    ----------
    params : dict
        Parameter dictionary to be mutated in place.

    Returns
    -------
    None
        Mutates 'gaussparanames' and 'hyperp_gausspriors' in params.
    """
    # gaussparanames: "", None -> []; "A_B" -> ["A","B"]; list/array -> list
    names = params.get("gaussparanames")
    if names in (None, ""):
        params["gaussparanames"] = []
    elif isinstance(names, str):
        params["gaussparanames"] = [n for n in names.split("_") if n != ""]
    elif isinstance(names, (list, tuple, np.ndarray)):
        params["gaussparanames"] = list(names)
    else:
        raise TypeError(f"Unsupported gaussparanames type: {type(names)}")

    # hyperp_gausspriors: "", None -> [];
    # "m1_s1|m2_s2" -> [[m1,s1],[m2,s2]];
    # [m,s] -> [[m,s]];
    # [[m1,s1],[m2,s2]] -> as-is
    pri = params.get("hyperp_gausspriors")
    if pri in (None, ""):
        params["hyperp_gausspriors"] = []
    elif isinstance(pri, str):
        chunks = [c for c in pri.split("|") if c.strip() != ""]
        params["hyperp_gausspriors"] = [[float(a), float(b)] for (a, b) in
                                        [p.split("_") for p in chunks]]
    elif isinstance(pri, (list, tuple)):
        if len(pri) == 2 and all(isinstance(x, (int, float)) for x in pri):
            params["hyperp_gausspriors"] = [[float(pri[0]), float(pri[1])]]
        else:
            params["hyperp_gausspriors"] = [[float(p[0]), float(p[1])] for p in pri]
    else:
        raise TypeError(f"Unsupported hyperp_gausspriors type: {type(pri)}")



def table_keys(name):
    """
    Get all table alias names given a section of the TOML file.

    Parameters
    ----------
    name : str
        Section name (e.g., "saving options").

    Returns
    -------
    tuple
        All possible alias names for that table.
    """
    return _TABLE_ALIASES.get(name, (name,))


def config_get(cfg, logical_table, key, default=None, required=False):
    """
    Unified accessor for config values from INI or TOML.

    Parameters
    ----------
    cfg : dict
        Dictionary of parsed configuration tables.

    logical_table : str
        Logical table name (e.g., "stellar params"). Underscore/space variants are handled.

    key : str
        Key to fetch from the table.

    default : Any, optional
        Default value if the key is missing (default: None).

    required : bool, optional
        If True, raises an error if the key is missing (default: False).

    Returns
    -------
    Any
        The value from the config or the default.
    """
    for t in table_keys(logical_table):
        if t in cfg and key in cfg[t]:
            return cfg[t][key]
    if required:
        raise KeyError(f"Missing key '{key}' in [{logical_table}]")
    return default


def load_as_dict(path):
    """
    Load an INI or TOML config file into a uniform dictionary.

    Parameters
    ----------
    path : str
        Path to the configuration file.

    Returns
    -------
    tuple
        (format, dict) where format is "ini" or "toml", and dict contains sections with keys/values.
    """
    ext = (os.path.splitext(path)[1].lower() or "")
    use_toml = ext in {".toml", ".tml", ".tom"}

    if use_toml:
        with open(path, "rb") as f:
            data = tomli.load(f)
        return "toml", data

    cp = configparser.ConfigParser()
    cp.read(path)
    sections = {}
    for sec in cp.sections():
        sections[sec] = {k: v for k, v in cp.items(sec)}
    return "ini", sections


def parse_all_input_exotune_iniFile(config_path):
    """
    Parse an INI or TOML configuration file and command-line style arguments into a parameter dictionary.

    This function provides backward compatibility for INI files while supporting
    the newer TOML format. Command-line arguments override values from the config file.

    Parameters
    ----------
    config_path : str
        Path to the configuration file (.ini or .toml).

    Returns
    -------
    tuple
        - params : dict
            Dictionary of parameter names and values, normalized for downstream use.
        - label_run : str
            A label string describing the run, built from selected parameters.
    """
    fmt, cfg = load_as_dict(config_path)

    # ---- Argparse with defaults from config ----
    parser = argparse.ArgumentParser(description='Inputs for stellar contamination retrieval.')

    # [choice of inputs]
    parser.add_argument('-label_obs', default=config_get(cfg, 'choice of inputs', 'label_obs'))
    parser.add_argument('-start_from_timeseries',
                        type=convert_to_bool, default=convert_to_bool(config_get(cfg, 'choice of inputs', 'start_from_timeseries')))
    parser.add_argument('-save_median_spectrum',
                        type=convert_to_bool, default=convert_to_bool(config_get(cfg, 'choice of inputs', 'save_median_spectrum')))
    parser.add_argument('-path_save_median_spectrum',
                        default=config_get(cfg, 'choice of inputs', 'path_save_median_spectrum', default=""))
    parser.add_argument('-path_to_stellar_spec_ts',
                        default=config_get(cfg, 'choice of inputs', 'path_to_stellar_spec_ts', default=""))
    parser.add_argument('-path_to_spec', default=config_get(cfg, 'choice of inputs', 'path_to_spec'))
    parser.add_argument('-spec_format', default=config_get(cfg, 'choice of inputs', 'spec_format'))
    parser.add_argument('-stmodfile', default=config_get(cfg, 'choice of inputs', 'stmodfile'))

    # [preprocessing]
    parser.add_argument('-optimize_param', type=convert_to_bool, default=convert_to_bool(config_get(cfg, 'preprocessing', 'optimize_param')))
    parser.add_argument('-obsmaskpattern', default=config_get(cfg, 'preprocessing', 'obsmaskpattern'))
    parser.add_argument('-kern_size', default=config_get(cfg, 'preprocessing', 'kern_size', default=""))
    parser.add_argument('-jd_range_mask', default=config_get(cfg, 'preprocessing', 'jd_range_mask', default=""))
    parser.add_argument('-wave_range_mask', default=config_get(cfg, 'preprocessing', 'wave_range_mask', default=""))

    # [saving options]
    parser.add_argument('-save_fit', type=convert_to_bool, default=convert_to_bool(config_get(cfg, 'saving options', 'save_fit')))
    parser.add_argument('-res_suffix', default=config_get(cfg, 'saving options', 'res_suffix'))

    # [stellar params]
    parser.add_argument('-Teffstar', type=float, default=float(config_get(cfg, 'stellar params', 'Teffstar')))
    parser.add_argument('-feh', type=float, default=float(config_get(cfg, 'stellar params', 'feh')))
    parser.add_argument('-loggstar', type=float, default=float(config_get(cfg, 'stellar params', 'loggstar')))
    parser.add_argument('-logg_phot_source', default=config_get(cfg, 'stellar params', 'logg_phot_source'))
    parser.add_argument('-logg_phot_value', type=float, default=float(config_get(cfg, 'stellar params', 'logg_phot_value')))

    # [stellar models]
    parser.add_argument('-label_grid', default=config_get(cfg, 'stellar models', 'label_grid'))
    parser.add_argument('-logg_range', default=config_get(cfg, 'stellar models', 'logg_range'))
    parser.add_argument('-loggstep', type=float, default=float(config_get(cfg, 'stellar models', 'loggstep')))
    parser.add_argument('-Teff_range', default=config_get(cfg, 'stellar models', 'Teff_range'))
    parser.add_argument('-Teffstep', type=float, default=float(config_get(cfg, 'stellar models', 'Teffstep')))
    parser.add_argument('-resPower_target', type=int, default=int(config_get(cfg, 'stellar models', 'resPower_target')))
    parser.add_argument('-wave_range', default=config_get(cfg, 'stellar models', 'wave_range'))

    # [MCMC params]
    parser.add_argument('-parallel', type=convert_to_bool, default=convert_to_bool(config_get(cfg, 'MCMC params', 'parallel')))
    parser.add_argument('-ncpu', type=int, default=int(config_get(cfg, 'MCMC params', 'ncpu')))
    parser.add_argument('-nsteps', type=int, default=int(config_get(cfg, 'MCMC params', 'nsteps')))
    parser.add_argument('-frac_burnin', type=float, default=float(config_get(cfg, 'MCMC params', 'frac_burnin')))
    # fitFscale is your original; accept fitDscale alias from your other module
    fitFscale_default = config_get(cfg, 'MCMC params', 'fitFscale')
    if fitFscale_default is None:
        fitFscale_default = config_get(cfg, 'MCMC params', 'fitDscale')  # alias
    parser.add_argument('-fitFscale', type=convert_to_bool, default=convert_to_bool(fitFscale_default))
    for param in ['fitspot', 'fitfac', 'fitThet', 'fitTphot', 'fitlogg_phot', 'fitlogg_het', 'fiterrInfl']:
        parser.add_argument(f'-{param}', type=convert_to_bool, default=convert_to_bool(config_get(cfg, 'MCMC params', param)))

    # [priors]
    parser.add_argument('-gaussparanames', default=config_get(cfg, 'priors', 'gaussparanames', default=""))
    parser.add_argument('-hyperp_gausspriors', default=config_get(cfg, 'priors', 'hyperp_gausspriors', default=""))
    parser.add_argument('-fitLogfSpotFac', default=config_get(cfg, 'priors', 'fitLogfSpotFac', default=""))
    parser.add_argument('-hyperp_logpriors', default=config_get(cfg, 'priors', 'hyperp_logpriors', default=""))

    # [plotting]
    parser.add_argument('-pad', type=float, default=float(config_get(cfg, 'plotting', 'pad')))
    parser.add_argument('-target_resP', type=int, default=int(config_get(cfg, 'plotting', 'target_resP')))

    # ---- CLI override parse ----
    args, _unknown = parser.parse_known_args()
    params = deepcopy(vars(args))

    # ---------- Normalization (works for INI strings and TOML arrays) ----------
    # masks
    for k in ("jd_range_mask", "wave_range_mask"):
        params[k] = parse_mask_ranges(params[k])

    # ranges
    params["logg_range"] = parse_range_pair(params["logg_range"], float)
    if isinstance(params["Teff_range"], str) and params["Teff_range"].lower() == "default":
        params["Teff_range"] = [np.max([params["Teffstar"] - 1000.0, 2300.0]), params["Teffstar"] + 1000.0]
    else:
        params["Teff_range"] = parse_range_pair(params["Teff_range"], float)
    params["wave_range"] = parse_range_pair(params["wave_range"], float)

    # gaussian priors: minimal pre-normalization, then delegate to sru.parse_gaussian_inputs
    preprocess_gaussians(params)
    names, priors = parse_gaussian_inputs(params)
    params["gaussparanames"], params["hyperp_gausspriors"] = names, priors

    # fitLogfSpotFac / hyperp_logpriors -> lists of length 2
    params["fitLogfSpotFac"] = parse_range_pair(params["fitLogfSpotFac"], int)
    params["hyperp_logpriors"] = parse_range_pair(params["hyperp_logpriors"], float)

    # kern_size
    if params["kern_size"] in (None, "", []):
        params["kern_size"] = None
    else:
        params["kern_size"] = int(params["kern_size"])

    # parallel guard
    if params["parallel"] is False:
        params["ncpu"] = 1

    # label
    label_run = f"Obs_{params['label_obs']}_Mod_{params['label_grid']}_Mask{params['obsmaskpattern']}"

    print("\nDone with parsing!")
    return params, label_run


def parse_gaussian_inputs(input_para):
    """
    Forwarder to the master implementation in stellar_retrieval_utilities.

    Parameters
    ----------
    input_para : dict

    Returns
    -------
    (gaussparanames: np.ndarray, hyperp_gausspriors: list[list[float]])
    """
    return sru.parse_gaussian_inputs(input_para)



## Pre-processing

def read_in_2d_spec(path):
    """
    Read a 2D time-series spectrum from an HDF5 file.

    This function reads an HDF5 file containing time-series spectral data
    (e.g., from a pipeline like Eureka!), extracts the 1D wavelength and time
    arrays, and reshapes the optimal extraction spectrum into a 2D array
    of shape (time, wavelength).

    Parameters
    ----------
    path : str
        Path to the HDF5 file containing spectral time series data.

    Returns
    -------
    timeUnique : ndarray
        1D array of unique time values.

    waveUnique : ndarray
        1D array of unique wavelength values.

    optspec2D : ndarray
        2D array of optimally extracted spectra with shape (n_times, n_wavelengths).
    """

    print("\nReading in time-series of spectra: ", path)

    s3output = h5py.File(path, 'r')
    waveUnique = np.unique(np.array(s3output["wave_1d"]))
    timeUnique = np.unique(np.array(s3output["time"]))

    print("\nReading in optimal extraction results...")
    optspec = np.array(s3output["optspec"])
    optspec2D = optspec.reshape(timeUnique.size, waveUnique.size)
    print("\nShape (time, wave) of spectra time series:", optspec2D.shape)
    return timeUnique, waveUnique, optspec2D


def calc_and_plot_timeMasked_mednorm_lc(timeUnique, optspec2D, params):
    """
    Calculate and plot a median-normalized light curve from a time-series spectrum,
    optionally applying smoothing and masking specified time ranges.

    Parameters
    ----------
    timeUnique : ndarray
        1D array of time values (e.g., BJD) for the spectral time series.

    optspec2D : ndarray
        2D array of spectral data with shape (n_times, n_wavelengths).

    params : dict
        Dictionary of parameters including:
            - "kern_size": (int or None) Kernel size for median filtering.
            - "mask_range": (tuple or list) One or more (start, end) time ranges to mask.
            - "jd_range_mask": list of lists or strings defining time intervals to be masked.

    Returns
    -------
    fig1 : matplotlib.figure.Figure
        The matplotlib figure containing the light curve plot.

    timeUnique : ndarray
        The time array after applying any specified masks.

    optspec2D : ndarray
        The 2D spectral array after masking the specified time intervals.
    """
    # create median-normalized light curve
    lc_mednorm = np.nansum(optspec2D / np.nanmedian(optspec2D, axis=0)[None, :], axis=1)
    lc_mednorm = lc_mednorm / np.nanmedian(lc_mednorm)

    if params["kern_size"] is not None:
        lc_medfilt = sig.medfilt(lc_mednorm, kernel_size=params["kern_size"])
        t_medfilt = sig.medfilt(timeUnique, kernel_size=params["kern_size"])

    fig1, ax = plt.subplots(1, 1)
    if params["kern_size"] is not None:
        bjd0 = int(t_medfilt[0])
        ax.plot(t_medfilt - bjd0, lc_medfilt, color="k")
    else:
        bjd0 = int(timeUnique[0])
        ax.plot(timeUnique - bjd0, lc_mednorm, color="k")

    if len(params["jd_range_mask"][0]):
        for i, mask_range in enumerate(params["jd_range_mask"]):
            ind = np.where((timeUnique < mask_range[0]) + (timeUnique > mask_range[1]))
            ax.axvspan(mask_range[0] - bjd0, mask_range[1] - bjd0, zorder=-1000, color="gray",
                       alpha=0.5, label=mask_range)
            timeUnique = timeUnique[ind]
            optspec2D = np.squeeze(optspec2D[ind, :])
            print("Masking time range ", mask_range, "new shape (time, wave):", optspec2D.shape)

    ax.legend()
    ax.set_xlabel("Time [BJD-" + str(bjd0) + "]")
    ax.set_ylabel("Flux (median-filt.)")
    ax.set_title("Masked regions (TIME)")
    return fig1, timeUnique, optspec2D

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


def get_StellSpec_from_spec_timeseries(waveUnique, optspec2D, params, flux_calibrated=True):
    """
    Compute the median spectrum from a spectral time series and initialize a StellSpec object.

    Parameters
    ----------
    waveUnique : ndarray
        1D array of unique wavelength values corresponding to the spectrum.

    optspec2D : ndarray
        2D array of optimal extraction spectral data with shape (n_times, n_wavelengths).

    params : dict
        Dictionary of parameters including:
            - "save_median_spectrum": (bool) Whether to save the median spectrum to file.
            - "path_save_median_spectrum": (str) Path to save the median spectrum file.

    flux_calibrated : bool, optional
        If True, assumes the spectral flux is already flux-calibrated. Default is True.

    Returns
    -------
    spec : psspec.StellSpec
        A `StellSpec` object containing the median spectrum.

    Notes
    -----
    - The spectral values and uncertainties are scaled by 1e10 for plotting purposes.
    """

    print("/nCalculate median spectrum and use to initialize StellSpec object...")
    median_spectrum, err_median_spectrum = get_median_spectrum_and_unc(waveUnique, optspec2D,
                                                                           flux_calibrated=flux_calibrated)

    specdict = dict()
    specdict["wave"] = waveUnique
    specdict["yval"] = median_spectrum * 1e10
    specdict["yerrLow"] = err_median_spectrum * 1e10
    specdict["yerrUpp"] = err_median_spectrum * 1e10
    t = table.Table(data=specdict)
    if params["save_median_spectrum"]:
        t.write(params["path_save_median_spectrum"], format="ascii.ecsv")
    spec = psspec.StellSpec(t, inputtype="basic")

    return spec

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

def prep_spec_waveMasked(spec, params):
    """
    Preprocess a spectrum object by masking invalid or unwanted wavelength regions.

    This function modifies a `StellSpec`-like spectrum object by:
    - Removing wavelengths beyond the red limit of the PHOENIX model grid.
    - Removing data points where the observed flux is NaN.
    - Applying custom wavelength range masks defined in the `params` dictionary.
    - Plotting the resulting masked spectrum with shaded regions indicating removed wavelength intervals.

    Parameters
    ----------
    spec : psspec.StellSpec
        A StellSpec-like object containing spectral data and methods for masking.

    params : dict
        Dictionary of parameters including:
            - "wave_range_mask": list of lists, where each item defines a wavelength range to mask
              (e.g., [[1.2, 1.3], [2.8, 3.0]]).

    Returns
    -------
    spec : StellSpec object
        The updated spectrum object with masked data removed.

    fig2 : matplotlib.figure.Figure
        The matplotlib figure showing the masked spectrum and shaded regions.

    Notes
    -----
    - The wavelength red limit (5.39 µm) is hardcoded based on the PHOENIX model grid.
    - NaN flux values are automatically identified and removed.
    - The function assumes `spec` provides `.wave`, `.yval`, `.remDataByIndex()`, `.waveMin`, `.waveMax`, and `.plot()` attributes/methods.
    """

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
    for i, mask_wv_range in enumerate(params["wave_range_mask"]):
        ind_to_keep = np.where((spec.waveMax<mask_wv_range[0])+(spec.waveMin>mask_wv_range[1]))
        ind_to_mask = np.where((spec.waveMax>=mask_wv_range[0])*(spec.waveMin<=mask_wv_range[1]))
        spec.remDataByIndex(ind_to_mask[0])
        ax.axvspan(mask_wv_range[0],mask_wv_range[1],zorder=-1000,color="gray",alpha=0.5)
        print("Masking wave range ",mask_wv_range, "in median spectrum prior to fit!")
        print("New wave shape:", spec.wave.shape)

    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel(r"Flux (median in time)")
    ax.set_title("Masked regions (WAVELENGTH)")
    return spec, fig2

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

def calc_and_plot_Fscale_guess(param, Teffs_grid, loggs_grid, wv_template_thisR, models_grid_fixedR, spec):
    print("\nGuessing Fscale value from the data...")
    Fscale_guess, model_int = get_scaling_guess(param, Teffs_grid, loggs_grid, wv_template_thisR,
                                                models_grid_fixedR, spec, wave_min_match_um=2.5,
                                                wave_max_match_um=3.0)

    # save the fixed R spectra as baseline specs to do the integration after the fact
    fl_phot_spot_fac = get_closest_models_from_grid(param, models_grid_fixedR, Teffs_grid, loggs_grid)

    print("\nPlotting default model with inferred initial guess Fscale applied...")
    # Plot model closest to the data
    fig3, ax = plt.subplots(1, 1, figsize=(10, 4))
    ax.plot(spec.wave, spec.yval, color="k", zorder=0)  # plot the observed spectrum
    ax.plot(wv_template_thisR, fl_phot_spot_fac[0] * Fscale_guess, color="r", label="Model", zorder=-1, alpha=0.2)
    ax.plot(spec.wave, model_int * Fscale_guess, color="r", label="Model in data bandpass", zorder=-1)
    ax.set_xlabel(r"Wavelength [$\mu$m]")
    ax.set_ylabel(r'Stellar flux [$\times$ 10$^{-10}$ erg/s/cm$^2$/$\mu$m]')

    ax.legend(loc=1)  # get nominal model from default parameters

    return Fscale_guess, fl_phot_spot_fac, fig3

## Modeling


def Fscale_obs(Fscale, fspot, ffac, model_phot_fixR, model_spot_fixR, model_fac_fixR):
    """
    Compute the model spectrum scaled by a flux factor.
    Wrapper around stellar_retrieval_utilities' combine_components function.

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
    return sru.combine_components("flux", Fscale, fspot, ffac,
                              model_phot_fixR, model_spot_fixR, model_fac_fixR)

def integ_model(lam_min,lam_max,wave,F):
    """
    Integrate a model over specified wavelength ranges.

    Computes the average flux within multiple wavelength intervals by integrating
    the flux density and normalizing by the width of each interval.

    Here simply a compatibility wrapper that defers to the other utilities function.
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

    return sru.integ_model(lam_min, lam_max, wave, F)


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

    Wrapper around the function compatible with both retrieval modalities in stellar_retrieval_utilities.

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

    return sru.get_closest_models_from_grid(param, models_grid_fixedR, Teffs_grid, loggs_grid)


## MCMC setup

def save_ref_files(this_dir, this_script, iniFile, utils_script, res_dir):
    """
    Save reference files used in a run to the results directory for reproducibility.

    This function copies the main script, the INI configuration file, and a utility
    script into a specified results directory, renaming them with standard names.

    Wrapper around `sru.save_ref_files` adapted to ExoTUNE's naming conventions.

    Parameters
    ----------
    this_dir : str
        Path to the directory where the current scripts and INI file are located.

    this_script : str
        Filename of the main run script (e.g., 'exotune_runscript.py').

    iniFile : str
        Filename of the INI configuration file used for the run.

    utils_script : str
        Path to the utility script used by the main script (e.g., 'exotune_utilities.py').

    res_dir : str
        Path to the directory where the reference copies will be saved.
    """


    names = dict(
        script_name="exotune_runfile_thisrun.py",
        utils_script_name="exotune_utilities_thisrun.py",
        ini_name="iniFile_exotune_thisrun.py",
    )
    return sru.save_ref_files(this_dir, this_script, iniFile, utils_script, res_dir, names=names)


def get_default_logg(params):
    """
    Get the default values for the stellar log g

    Parameters
    ----------
    params: dict
    Dictionary obtained when reading in the ini file
    This dictionary is updated with logg_phot_default

    Returns
    -------
    None

    """

    # logg for MCMC
    if params["logg_phot_source"] == "loggstar":
        logg_phot = deepcopy(params["loggstar"])
    elif params["logg_phot_source"] == "value":
        logg_phot = deepcopy(params["logg_phot_value"])

    params["logg_phot_default"] = logg_phot

def init_default_and_fitted_param(Tphot, met, logg_phot, fitspot=True, fitfac=True,
                                  fitThet=False, fitTphot=False, fitFscale=True,
                                  fitlogg_phot=False, fitlogg_het=False,
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
    fitlogg_het : bool, optional
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
    if fitlogg_het:
        fitparanames.append("dlogghet")

    return param, fitparanames

def init_default_and_fitted_param_fromDict(params, Fscale_guess=1, dlogg_het_guess=0.0):
    """
    Initialize default stellar parameters and define which are to be fitted,
    using a parameter dictionary input (e.g. parsed from an INI file).

    Parameters
    ----------
    params : dict
        Dictionary containing model parameters and fit flags

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
    fitLogfSpotFac = params["fitLogfSpotFac"]

    return init_default_and_fitted_param(
        Tphot=params["Teffstar"],
        met=params["feh"],
        logg_phot=params["logg_phot_default"],
        fitspot=params.get("fitspot", True),
        fitfac=params.get("fitfac", True),
        fitThet=params.get("fitThet", True),
        fitTphot=params.get("fitTphot", True),
        fitFscale=params.get("fitFscale", True),
        fitlogg_phot=params.get("fitlogg_phot", False),
        fitlogg_het=params.get("fitlogg_het", False),
        fitLogfSpotFac=fitLogfSpotFac,
        fiterrInfl=params.get("fiterrInfl", False),
        Fscale_guess=Fscale_guess,
        dlogg_het_guess=dlogg_het_guess
    )

def get_derived_param(param):
    """
    Wrapper around stellar_retrieval_utilities.get_derived_param to calculate derived parameters
    Parameters
    ----------
    param : dict
        dictionary of default param values .

    Returns
    -------
    param with derived parameters updated
    """
    return sru.get_derived_param(param)


def get_param_priors(param, gaussparanames, hyperp_gausspriors=[],
                         fitLogfSpotFac=[False, False], hyperp_logpriors=[],
                         T_grid=[2300., 10000.], logg_grid=[2.5, 5.5],
                         Fscale_guess=1.0):
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
        For parameters fitted in log-space, the key is prefixed with 'log' followed by underscore.
    """
    return sru.get_param_priors_master(
        param, gaussparanames, hyperp_gausspriors,
        fitLogfSpotFac, hyperp_logpriors,
        T_grid, logg_grid,
        mode="flux", Fscale_guess=Fscale_guess)

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




def save_mcmc_to_pandas(*args, **kwargs):
    """
    Save MCMC sampling results and extract best-fit parameters.

    This function processes the MCMC sampler output by removing the burn-in phase,
    reshaping and storing the flattened samples as a Pandas DataFrame and an Astropy Table.
    It also computes and saves the best-fit parameters based on the maximum likelihood.

    This is simply a wrapper around `stellar_retrieval_utilities.save_mcmc_to_pandas`

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
    return sru.save_mcmc_to_pandas(*args, **kwargs, prefix="exotune")

def samp2bestfit(samp):
    """
    Compute summary statistics and best-fit values from MCMC samples.

    Calculates key quantiles (median, 1 and 2 sigma) and identifies the samples
    corresponding to the maximum posterior probability and maximum likelihood.

    Wrapper around `stellar_retrieval_utilities.samp2bestfit` (identical behavior).

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

    return sru.samp2bestfit(samp)

def BIC(chi2,nDataPoints,nPara):

    """
    Calculate the Bayesian Information Criterion (BIC) for model comparison.

    Simple wrapper around `stellar_retrieval_utilities.BIC` to compute the BIC value
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

    return sru.BIC(chi2, nDataPoints, nPara)

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

    return sru.save_bestfit_stats(spec, ind_bestfit, fitparanames,
                   flat_oot_spec_models, results_folder, runname, prefix="exotune",save_fit=save_fit)


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

    return sru.get_percentile_blobs(exotune_models_blobs, percentiles=percentiles)



##  Plotting utilities

def xspeclog(ax, xlim=None, level=1, fmt="%2.1f"):
    """
    Wrapper around `stellar_retrieval_utilities.xspeclog` to set logarithmic x-axis ticks
    Configure logarithmic x-axis ticks for spectrum plots.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matplotlib axis object on which to set the log-scale x-axis ticks.
    xlim : tuple of float, optional
        Optional (xmin, xmax) tuple to explicitly set the x-axis limits.
    level : int or float, default=1. Preset level for tick density and appearance:
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

    sru.xspeclog(ax, xlim=xlim, level=level, fmt=fmt)


def setAxesFontSize(ax, fontsize):
    """
    Wrapper around `stellar_retrieval_utilities.setAxesFontSize` to set font sizes for axis labels and ticks.

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

    sru.setAxesFontSize(ax, fontsize)



    
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


    return sru.get_labels_from_fitparanames_master(fitparanames, include_units=True)

def wrap_text(text, nchar=60):
    """
    Splits a long string (text, str.) by inserting a line break every nchar (int.) characters,
    ideally breaking at spaces for better readability.
    Parameters
    ----------
    text : str
        The input string to be wrapped.
    nchar : int, default=60
        Maximum number of characters per line before inserting a line break.

    Returns
    -------
    str
        The input string with line breaks inserted for wrapping.
    """
    words = text.split()
    lines = []
    current_line = ""

    for word in words:
        if len(current_line + ' ' + word) <= nchar:
            current_line += ' ' + word if current_line else word
        else:
            lines.append(current_line)
            current_line = word

    if current_line:
        lines.append(current_line)

    return '\n'.join(lines)

def format_param_str(param, fitparanames):
    """
    Generate a formatted summary string of stellar and fitted parameters.

    Wrapper around `stellar_retrieval_utilities.format_param_str`
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


    return sru.format_param_str(
        param,
        fitparanames,
        mode="flux",
        include_fit_list=True,
        wrap_at=50,   # matches prior XTU wrap_text(…, nchar=50)
    )


def chainplot(samples, labels=None, nwalkers=None, fontsize=None, lw=1, stepRange=None):
    """
    Plot parameter chains from MCMC sampling results.

    Wrapper around `stellar_retrieval_utilities.chainplot`
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
    return sru.chainplot(
            samples,
            labels=labels,
            nwalkers=nwalkers,
            fontsize=fontsize,
            lw=lw,
            stepRange=stepRange,
        )


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

    return sru.plot_blobs_spectrum(
        spec,
        oot_models_blobs,
        ind_bestfit,
        ax=ax,
        kind="exotune",
        bestfit_color=bestfit_color,
        color=color,
        plot3sig=plot3sig,
        plot2sig=plot2sig,
        plot1sig=plot1sig,
        plotmedian=plotmedian,
        plotbestfit=plotbestfit,
        legend_loc=legend_loc,
        save_csv=save_csv,
        results_folder=results_folder,
        runname=runname)

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
                        quantiles=[0.16, 0.5, 0.84], title_kwargs={'fontsize':10},
                        hist_kwargs={"linewidth":3}, rg=None, 
                        show_titles=True, levels=(0.393,0.865,0.989), **kwargs):
    """
    Generate a corner plot to visualize posterior distributions from MCMC samples.

    This function wraps the `corner.corner` plotting and uses the master function from stellar_retrieval_utilities.py

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

    return sru.plot_corner_master(
        samples,
        labels=plotparams["labels"],
        color=plotparams["hist_color"],
        plot_datapoints=plot_datapoints,
        smooth=smooth,
        quantiles=quantiles,
        title_kwargs=title_kwargs,
        hist_kwargs=hist_kwargs,
        rg=rg,
        show_titles=show_titles,
        levels=levels,
        **kwargs
    )

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

    print("\nProceed to plotting corner plot...")
    fig = plot_corner(ordered_samples, plotparams, smooth=0.8, fill_contours=True,
                      truths = ordered_parabestfit, truth_color="k")
    # fig = plot_corner(ordered_samples[:,ind], plotparams, smooth=0.8, fill_contours=True, 
    #                   truths = ordered_parabestfit[ind], truth_color="k")    
                                    

    fig.set_dpi(50)
    fig.set_figheight(11)
    fig.set_figwidth(15)
    # pdb.set_trace()


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
    plt = ax.figure  # avoid flake warnings
    # get y-range after base plot
    ymin, ymax = ax.get_ylim()
    # where to write: slightly inside the right padding, near top
    text_x = np.max(spec.waveMax) + pad - pad / 4
    text_y = 0.95 * ymax

    return sru.plot_maxlike_and_maxprob(
        spec, param, parabestfit, ind_maxprob, ind_bestfit, fitparanames, flat_st_ctm_models,
        pad=pad,legend_loc=4,title="Best-fit model",best_label="Max. Likelihood",prob_label="Max. Probability",
        best_color="r",prob_color="slateblue",best_alpha=0.5,prob_alpha=1.0,best_marker=".",prob_marker="o",
        best_ms=50,prob_ms=30,show_param_text=True,formatter=format_param_str,
        include_units=True,  text_xy=(text_x, text_y), text_fontsize=8,text_color="k",text_bbox=None,text_ha="right",
        text_va="top")


## Multiprocessing tools

class FunctionCache:
    """
    A global container for storing a target function and its associated arguments,
    intended to minimize pickling overhead when using multiprocessing.

    Attributes
    ----------
    lnprob : callable
        The target function to be evaluated in parallel.
    lnprob_args : tuple
        Positional arguments to pass to the function.
    lnprob_kwargs : dict
        Keyword arguments to pass to the function.
    """
    lnprob = None
    lnprob_args = ()
    lnprob_kwargs = {}

def initializer(lnprob, lnprob_args=(), lnprob_kwargs={}):
    """
    Initializes the global FunctionCache with a function and its arguments.

    Parameters
    ----------
    lnprob : callable
        The function to evaluate.
    lnprob_args : tuple, optional
        Positional arguments for the function.
    lnprob_kwargs : dict, optional
        Keyword arguments for the function.
    """
    FunctionCache.lnprob = lnprob
    FunctionCache.lnprob_args = lnprob_args
    FunctionCache.lnprob_kwargs = lnprob_kwargs

def lnprob_cache(theta):
    """
    Evaluates the cached function with the provided input.

    Parameters
    ----------
    theta : any
        The input to pass as the first argument to the cached function.

    Returns
    -------
    The result of evaluating the cached function.
    """
    return FunctionCache.lnprob(theta, *FunctionCache.lnprob_args, **FunctionCache.lnprob_kwargs)

class FunctionWrapper:
    """
    A simple function wrapper to bundle a function and its shared arguments,
    useful for serialization-safe function calls.

    Attributes
    ----------
    func : callable
        The function to wrap.
    args : tuple
        The shared arguments to pass to the function.
    """
    def __init__(self, func, args):
        """
        Parameters
        ----------
        func : callable
            The function to be called.
        args : tuple
            The shared arguments to apply with the function.
        """
        self.func = func
        self.args = args

    def __call__(self, theta):
        """
        Calls the wrapped function with the provided input and shared arguments.

        Parameters
        ----------
        theta : any
            The input to the function.

        Returns
        -------
        The result of the function call.
        """
        return self.func(theta, *self.args)

class MyPool:
    """
    A lightweight multiprocessing pool that avoids pickling large argument objects
    by using a shared global cache (`FunctionCache`).

    Methods
    -------
    map(func, iterable)
        Applies the cached function over an iterable using worker processes.
    close()
        Easy stop of pool.
    terminate()
        Force stops the pool.
    """
    def __init__(self, processes, lnprob, lnprob_args=(), lnprob_kwargs={}):
        """
        Parameters
        ----------
        processes : int
            Number of worker processes to use.
        lnprob : callable
            The function to be applied in parallel.
        lnprob_args : tuple, optional
            Positional arguments for the function.
        lnprob_kwargs : dict, optional
            Keyword arguments for the function.
        """
        self.pool = mp.Pool(
            processes=processes,
            initializer=initializer,
            initargs=(lnprob, lnprob_args, lnprob_kwargs)
        )

    def map(self, func, iterable):
        """
        Applies the cached function to every item in the iterable using the pool.

        Parameters
        ----------
        func : callable
            Ignored; required to match map() signature.
        iterable : iterable
            Inputs to evaluate the function on.

        Returns
        -------
        list
            List of function results.
        """
        return self.pool.map(lnprob_cache, iterable)

    def close(self):
        """
        Closes the pool and waits for all worker processes to exit.
        """
        self.pool.close()
        self.pool.join()

    def join(self):
        """
        Waits for the worker processes to exit.
        """
        self.pool.join()

    def terminate(self):
        """
        Immediately terminates all worker processes.
        """
        self.pool.terminate()

    def __enter__(self):
        """
        Enables use of `with` statement for context management.

        Returns
        -------
        MyPool
            The pool instance itself.
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Ensures that the pool is closed when exiting a `with` block.
        """
        self.close()
