.. _running_exotune:

Running stellar spectrum retrievals with *exotune*
==================================================

Now that you have installed *stctm* and set up your environment, you are ready to run a retrieval on stellar spectra!

You can fit any stellar spectrum (loaded into the code as a ``StellSpec`` object) assuming it can be represented by a linear combination of up to three components: the photosphere, cooler regions (spots), and hotter regions (faculae). The code allows you to fit the contributions of spots and/or faculae, vary or fix the temperatures of the heterogeneity components, fit the surface gravity used for photosphere and/or heterogeneity models, and control priors on all fitted parameters. Gaussian and uniform priors are available, including log-uniform priors for heterogeneity covering fractions. Posterior samples, model spectra, diagnostic statistics, and publication-ready plots are generated.

The analysis also supports parallel processing, allowing you to run the retrieval efficiently on multiple cores or large clusters, which is especially useful for high-resolution or large-coverage spectral datasets.


The data
--------

The example dataset provided under ``example/observations/`` for *exotune* is a spectrum of TRAPPIST-1 (out-of-transit) observed with JWST NIRSpec/PRISM (Piaulet-Ghorayeb et al., 2025). The units are microns for wavelength and :math:`10^{-10}` erg/cm\ :sup:`2`/:math:`\mu m` for flux.

Your data file will need to be readable using one of the supported formats in ``stctm/pystellspec.py`` (see the ``StellSpec`` class definition below). The formats implemented so far are ``basic`` or ``MR_csv``— but you can add a new handler if necessary by extending the ``StellSpec`` constructor.


For reference, here is how ``StellSpec`` objects are initialized from a table in ``stctm/pystellspec.py``::

    elif inputtype=='MR_csv':
        inputtable["yval"] = inputtable["spec"]*1e10
        inputtable["waveMin"] = inputtable["wave_low"]
        inputtable["waveMax"] = inputtable["wave_high"]
        inputtable["yerrLow"] = inputtable["err"]*1e10
        inputtable["yerrUpp"] = inputtable["err"]*1e10
        inputtable["wave"] = 0.5*(inputtable["waveMin"] + inputtable["waveMax"])
        super(StellSpec,self).__init__(inputtable)
        for c in self.colnames:
            setattr(self, c, self[c])



* To use a custom format, add a new ``elif inputtype=='myformat':`` block, and transform your table columns accordingly, making sure that the ``yval`` (flux value in :math:`10^{-10}` erg/cm\ :sup:`2`/:math:`\mu m`), ``yerrLow`` and ``yerrUpp`` (lower and upper flux errors) are specified, as well as each bin's central wavelength (``wave`` in microns) and the minimum/maximum (wavelength edges) for each bin: ``waveMin`` and ``waveMax``.
* You may read your data using any compatible astropy table format (see the `Astropy Table I/O documentation <https://docs.astropy.org/en/stable/table/io.html>`_).
* Ensure your wavelength array, ``waveMin``, and ``waveMax`` are in microns (convert if necessary).

The recommended approach is to provide your spectrum file in one of these formats and then instantiate your object by specifying the format (via the ``spec_format`` parameter in your ``.ini`` file, see below).

In summary, as long as your input file fits one of these supported formats or you extend ``StellSpec`` for your format, you will be able to analyze your observed stellar spectrum with *exotune*.

Setting up an *exotune* retrieval: Run instructions
---------------------------------------------------

The ``.ini`` file structure largely mirrors that of TLS retrievals (see :ref:`runningTLS` for folder layout). The key exotune-specific modifications are:

* The analysis folder is ``exotune_analysis/any_analysis_folder_name/``. Create a new folder for each project under ``exotune_analysis/``, e.g., ``template_exotune_analysis/``.
* The default INI file is ``template_ini_exotune.ini`` which you can use as a starting point to build your own.
* Results are saved in ``exotune_results/``, a sibling folder to ``exotune_analysis/``. Each run produces a subfolder in ``exotune_results/``.
* Ensure all environment variables and paths are correctly set (refer to comments at the top of the run script for requirements), copied here::

    # make sure your environment variables are set up (see example below)
    # os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
    # os.environ['CRDS_PATH'] = "/home/caroline/crds_cache"
    # os.environ['PYSYN_CDBS'] = "/home/caroline/trds"

Example command to run exotune using the CLI helper, after navigating to your analysis folder::

    stctm_exotune template_ini_exotune.ini

Inputs can be edited in the INI file or passed as command-line arguments. For instance, to switch off ``fitspot``::

    stctm_exotune template_ini_exotune.ini -fitspot=0

Modifying the INI file to make it your own
------------------------------------------

The sections and parameters below correspond directly to the exotune ``.ini`` file.

Choosing inputs and starting format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Contrary to ``stctm``, for ``exotune`` you can either start from a pre-processed existing flux-calibrated stellar spectrum or from a time series of spectra. If the latter is the starting point, then *exotune* will proceed to some pre-processing steps to compute the "master" out-of-transit stellar spectrum given the pre-processing options detailed in the subsection below.

Under ``[choice of inputs]``::

    [choice of inputs]

    label_obs = test_visit
    start_from_timeseries = False
    save_median_spectrum = False
    path_save_median_spectrum = ../../observations/planetname/planet_outoftransit_spectrum.csv
    path_to_stellar_spec_ts =
    path_to_spec = ../../observations/TRAPPIST_1_NIRSpec/exotune_templatespectrum.csv
    spec_format = basic
    stmodfile = ../../R10000_model_grids/TRAPPIST_1_pymsg.h5

* ``label_obs``: Short label for this dataset for tracking results and outputs.
* ``start_from_timeseries``: ``True`` if starting from a time series of spectra, ``False`` for a pre-processed single spectrum file (see below for pre-processing).
* ``save_median_spectrum``: If starting from a timeseries, set ``True`` to save the computed median spectrum.
* ``path_save_median_spectrum``: Output path for the median spectrum CSV (used only if ``save_median_spectrum`` is ``True``).
* ``path_to_stellar_spec_ts``: Path to the time series file (used if ``start_from_timeseries`` is ``True``).
* ``path_to_spec``: Path to single spectrum file (used if ``start_from_timeseries`` is ``False``).
* ``spec_format``: Spectrum format string for loading into ``StellSpec``. See ``pystellspec.py`` for supported formats or to add a custom format.
* ``stmodfile``: Path to the stellar models grid file (HDF5).

Preprocessing options
^^^^^^^^^^^^^^^^^^^^^

Under ``[preprocessing]``::

    [preprocessing]
    optimize_param = False
    obsmaskpattern= nomask
    kern_size = 19
    jd_range_mask =
    wave_range_mask =

* ``optimize_param``: ``True`` to only preprocess and visualize (no MCMC, just diagnostic plots).
* ``obsmaskpattern``: Label used for the specific mask pattern (will be used as a string when saving the run).
* ``kern_size``: Kernel size for median filtering the plotted light curve (for visualization only).
* ``jd_range_mask``: Custom time-domain mask. To make sure that some intervals of time are ignored, e.g. in-transit, or during a stellar flare, enter their time stamps as ``start1_end1|start2_end2|...``.
* ``wave_range_mask``: Custom wavelength-domain mask, same format as above.

Saving options
^^^^^^^^^^^^^^

Under ``[saving options]``::

    [saving options]
    save_fit = True
    res_suffix = test_for_GitHub

* ``save_fit``: ``True`` to save results in the output directory after completion.
* ``res_suffix``: Suffix tagging the output files for identification; change for each new run.

Stellar parameters
^^^^^^^^^^^^^^^^^^

Under ``[stellar params]``::

    [stellar params]
    Teffstar = 2566
    feh = 0.040
    loggstar = 5.2396
    logg_phot_source = value
    logg_phot_value = 2.5

* ``Teffstar``: Effective temperature of the star in Kelvin.
* ``feh``: Metallicity [Fe/H] in dex.
* ``loggstar``: Surface gravity log(g) in cgs.
* ``logg_phot_source``: ``value`` to use ``logg_phot_value`` for the photosphere log(g) default, ``loggstar`` to use the star value.

Reading in the grid of stellar models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under ``[stellar models]``::

    [stellar models]
    label_grid = PHOENIX_TRAPPIST_1
    logg_range = 2.5_5.5
    loggstep = 0.1
    Teff_range = default
    Teffstep = 20.
    resPower_target = 10000
    wave_range = 0.2_5.4

* ``label_grid``: Name/label of the stellar model grid (used as a string to save the run).

At this stage, refer to your ``create_fixedR_grid_pymsg_template.py`` file (or the equivalent file you used to create your grid of stellar models).
In that file, you will find the setup of the grid in a block such as::

    # range of params for the grid

    logg_range = [2.5,5.5]
    Teff_range = [np.max([param["Tphot"]-1000, 2300.]), param["Tphot"]+1000.]
    loggstep = 0.1 #cgs
    Teffstep = 20. #K
    resPower_target = 10000
    wv_min_um = 0.2
    wv_max_um = 5.4

Returning to the ``.ini`` file:
* ``logg_range``: Range of log(g) covered in the grid(format ``minlogg_maxlogg``).
* ``loggstep``: Grid step in log(g).
* ``Teff_range``: Temperature range; ``default`` uses values calculated from ``Teffstar``: it assumes the default grid calculation setup, with`` min = np.max([Teffstar-1000, 2300.]) `` and ```max=Teffstar+1000``.
* ``Teffstep``: Grid step in temperature.
* ``resPower_target``: Resolving power at which the grid was created.
* ``wave_range``: Wavelength range for fitting (microns, ``min_max`` format).

MCMC sampling parameters
^^^^^^^^^^^^^^^^^^^^^^^^

Under ``[MCMC params]``::

    [MCMC params]
    parallel = True
    ncpu = 30
    nsteps = 3000
    frac_burnin = 0.6
    fitspot = True
    fitfac = True
    fitThet = True
    fitTphot = True
    fitlogg_phot = True
    fitlogg_het = True
    fitFscale = True
    fiterrInfl = True

* ``parallel``: Use multiprocessing (``True`` recommended).
* ``ncpu``: Number of CPUs for parallel MCMC run.
* ``nsteps``: Number of steps for each MCMC chain (recommend 5000+ for analysis).
* ``frac_burnin``: Fraction of chain steps discarded as burn-in (e.g., ``0.6``).
* ``fitspot`` / ``fitfac``: Whether to fit spot/faculae covering fractions.
* ``fitThet`` / ``fitTphot``: Whether to fit spots/faculae/photosphere temperature.
* ``fitlogg_phot`` / ``fitlogg_het``: Whether to fit log(g) for photosphere and/or heterogeneity.
* ``fitFscale``: Fit a flux scaling factor to match observed/model spectra.
* ``fiterrInfl``: Fit an error inflation factor to relax the provided data error bars if model/data mismatch is large.

Priors on the fitted parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under ``[priors]``::

    [priors]
    gaussparanames = Tphot
    hyperp_gausspriors = 2566_70
    fitLogfSpotFac = 0_0
    hyperp_logpriors = -5_0

* ``gaussparanames``: List of parameters to apply a Gaussian prior (separated by underscores, e.g. ``Tphot_ffac``).
* ``hyperp_gausspriors``: Mean and std for each Gaussian prior. For multiple parameters separate with a vertical line: e.g. ``mean1_std1|mean2_std2``

* ``fitLogfSpotFac``: Specifies if spot/faculae priors are uniform in linear (toggle ``0``) or log space (toggle ``1``).
* ``hyperp_logpriors``: Bounds for log-priors (``lowerBound_upperBound``).

Beyond the flexibility provided in the ``.ini`` file, you can look up the logic in ``get_param_priors()`` in ``stctm/exotune_utilities.py``.

Plotting
^^^^^^^^

Under ``[plotting]``::

    [plotting]
    pad = 0.25
    target_resP = 300

* ``pad``: Padding in microns to adjust spectra plot axis boundaries.
* ``target_resP``: Resolving power model spectra are downgraded to when plotted.

Post-processing
---------------

By default, *exotune* generates and saves the following files to a custom directory created under ``exotune_results/``, starting with the prefix "fit". If only preprocessing is run (i.e., ``optimize_param = True``), or if the starting point is a time-series of spectra, the results from this preprocessing step will be in a folder starting with "preprocessOnly".

Inputs and recordkeeping:

- Copy of run script, INI file, and ``exotune_utilities.py`` used
- Figure of the fitted spectrum
- ``defaultparams`` CSV file with fit initial values

Pre-processing steps:

- ``select_time``: Median-filtered light curve marking masked intervals used in generating the median spectrum
- ``select_wave``: Median spectrum before masking, with masked wavelength intervals shaded
- ``get_fscale``: Initial model/data comparison used for scavenging the starting value of ``Fscale``

Outputs (CSV files):

- ``pandas``: Fitted parameters from chain, with log-likelihood and log-probability
- ``bestfit``: Best-fit value (maximum likelihood), max-probability, and percentiles for quoting
- ``bestfit_stats``: Model comparison statistics: best-fit model index, reduced chi-squared, and BIC
- ``fixedR_1_2_3_sigma``: Model spectra at the plotted resolving power (default ``target_resP``) for max-likelihood, max-probability, and percentile intervals
- ``blobs_1_2_3_sigma``: Model spectra integrated in observed data bins for max-likelihood, max-probability, and percentiles

Calculated models:

- NPY file containing "blobs": the series of models from MCMC sampling

Diagnostics figures:

- ``chainplot``: Chain plots, before and after burn-in
- ``bestfit_model``: Plot of the best-fit model over data

Publication-ready figures:

- ``resP..._1_2_3_sigma``: Fitted spectra with 1/2/3 sigma intervals at high resolution (resolving power ``target_resP``), log or lin wavelength axis
- ``combo_resP..._1_2_3_sigma``: Top: fitted spectrum and intervals; Bottom: marginalized posterior distributions for component parameters
- ``1_2_3_sigma``: Fitted spectrum with intervals using data bin integration
- Corner plot of post-burnin samples

Please let me know (or create a pull request!) if there are additional outputs that would be useful defaults.