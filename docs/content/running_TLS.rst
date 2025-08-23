.. _runningTLS:

Running your first TLS retrieval
================================


Now that you have properly installed *stctm* and run the test case, you are ready to run your first TLS retrieval!

You can fit any transmission spectrum (loaded into the code as a ``TransSpec`` object) assuming any variations from a straight line are due to the transit light source effect. In its current configuration, you can fit the contributions of spots and/or faculae, varying or fixing the temperatures of the heterogeneity components, as well as fit the surface gravity used for the photosphere and/or heterogeneity models. The code has the flexibility to handle Gaussian priors on any fitted parameter as well as linear or log-uniform priors on the heterogeneity covering fractions. You obtain a range of outputs including posterior samples (parameters, spectra), run statistics for model comparison, and publication-ready plots.

I will walk you through how to set up your ``.toml`` file from the example file provided in order to fit any of your needs, how to run the retrieval, and what to expect in terms of baseline outputs (that can be used for later post-processing).



The data
--------


The example dataset provided under ``example/observations/`` for *stctm* is a transmission spectrum of TRAPPIST-1 b (visit 1 from Lim et al., 2023 where I ran *stctm*) observed with JWST NIRISS/SOSS. The units in that file are millimeters for the wavelength (hence the choice of ``wavemm`` for ``spec_format``, to convert mm to the microns unit used throughout *stctm*), and ppm for the transit depth.

Your data file will have to be readable by one of the options in ``stctm/pytransspec.py``.

The best way to proceed is to copy a code block in the ``__init__`` function in ``stctm/pytransspec.py``::

    elif inputtype=='wavemm':
        table=Table.read(inputpath,format='ascii.ecsv',header_start=header_start)

        super(TransSpec,self).__init__(table)
        self["wave"] = self["wave"]*1e3 # convert to microns from millimeters
        self["waveMin"] = self["waveMin"]*1e3
        self["waveMax"] = self["waveMax"]*1e3

        self.sort(keys=["waveMin"])
        self['wave'].unit=waveunit
        self['waveMin'].unit=waveunit
        self['waveMax'].unit=waveunit

        #Metadata
        self.meta['waveunit']=waveunit
        self.meta['label']=label
        self.meta['color']=color

* replace "wavemm" in the example with the way you want to refer to your format (e.g. ``elif inputtype=='myformat':``)

* read in the data file using the appropriate ``Table.read()`` command (you can check the `astropy.table` documentation for more information on how to read in different formats: https://docs.astropy.org/en/stable/table/io.html)

The code will then create an attribute to the object for each column in the table (``super()`` command).

* make sure the wavelength array, as well as ``waveMin`` and ``waveMax`` (the min and max wavelength for each wavelength bin) are in microns (convert if necessary as above)
* you can simply copy-paste the other lines with units and color.

Run instructions
----------------

You should follow the following folder structure, starting with a copy of the ``stctm/example/`` directory anywhere in your installation you may want to run the code.

* ``stellar_contamination_analysis/any_analysis_folder_name/``: my advice is to create a new analysis folder name under the ``stellar_contamination_analysis`` folder for each project you work on. In the example, that folder is called ``template_analysis/``.
* In that folder, you'll need to have at least one ``.toml`` file, where you specify all the inputs following the instructions in the comments (see more information below).
* At the same level as ``stellar_contamination_analysis/``, create a folder called ``stellar_contamination_results/``. For each of your runs, a subfolder will be added under this results directory and the run results will be saved there.

Here is an example one-liner to run a *stctm* retrieval from a planet spectrum, after navigating to ``stellar_contamination_analysis/template_analysis`` and in the command line::

    stctm_TLS template_ini_stctm.toml

A few additional tips:

- The "TOML file" (.toml) contains all the human-readable user input necessary to run the code, including stellar and mcmc parameters, saving paths and plotting options (see below for details of setting it up)
- The path to a TOML file needs to be provided as an argument if the script is run from the command line
- Any parameter in the TOML file can alternatively be modified from the default using the command line (instead of modifying the file). For instance, if you want to run the same fit as above, but only modify the suffix used for creating the output directory, you can do it as follows::

    stctm_TLS template_ini_stctm.toml -res_suffix=second_test

Modifying the TOML file to make it your own
-------------------------------------------

All the inputs you need to provide are specified in the ``.toml`` file.

Setting up labels and path to spectrum file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under ``[paths_and_labels]``  the TOML file will look like this::

    [paths_and_labels]

    path_to_spec = "../../observations/TRAPPIST_1_b_NIRISS_SOSS/Visit1_Order1And2.spec"
    spec_format = "wavemm"
    stmodfile = "../../R10000_model_grids/TRAPPIST_1_pymsg.h5"
    save_fit = true
    res_suffix = "TRAPPIST_1_b_testINIexample"


You can edit the set up as follows:

* ``path_to_spec`` (path to your spectrum file) as well as ``spec_format`` (your spectrum is read in from your data file as a ``TransSpec`` object using the ``spec_format`` setting you choose). You can choose `basic` for `spec_format` if your spectrum already has all the right column names and wavelength in microns, or `wavemm` if the wavelengths are in millimeters - if you are not sure which option to choose, or need to add another option to read in your specific format, you can do so in ``stctm/pytransspec.py`` as documented above!
* ``stmodfile``: the path to your stellar models grid file, in the HDF5 format
* ``save_fit``: ``True`` to save files to the results directory during the post-processing steps.
* ``res_suffix``: a suffix used for all the files that will be saved as a result of this run, in the results folder. This is the identifier you can use to record information on the spectrum, the setup of the fit, etc: make sure it is unique to avoid overwriting the contents of your results folder!

Setting up the stellar parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under ``[stellar_params]``, you have the following options::

    [stellar_params]

    Teffstar = 2566
    feh = 0.040
    loggstar = 5.2396

    logg_phot_source = "loggstar"
    logg_phot_value = 5
    logg_het_default_source = "logg_phot"
    logg_het_value = 5


* Enter the parameters of the star (effective temperature, metallicity Fe/H, log g in cgs) to set the defaults for the fit.

This is how to set up (potentially distinct) default values for the stellar and heterogeneity log g:

* ``logg_phot_source``: ``value`` to use the value of ``logg_phot_value`` as the stellar photosphere log g, otherwise ``loggstar`` to use the value provided in the code block below containing the stellar parameters;
* ``logg_het_default_source``: ``value`` to use the value of ``logg_het_value`` as the heterogeneities (default, if fitted) log g, otherwise ``logg_phot`` to set it to the same value as the stellar photosphere log g.

Reading in the grid of stellar models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under ``[stellar_models]``, you will find the following options::

    [stellar_models]

    logg_range = [2.5,5.5]
    loggstep = 0.1

    # options are default or min_max. default assumes the default grid calculation setup, with
    # min = np.max([Teffstar-1000, 2300.]) and max=Teffstar+1000.
    Teff_range = "default"

    Teffstep = 20.0
    resPower_target = 10000
    wave_range = [0.2,5.4]

As indicated in the comment, setting Teff_range to "default" means ``min = np.max([Teffstar-1000, 2300.])`` and  ``max=Teffstar+1000.``. For custom values, set the parameter to ``[min, max]``.
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

This is the information you need to take and paste into your ``.toml`` file under the ``[stellar models]`` section.
In particular, make sure to modify the range and spacing of the grid in the log g and Teff dimensions to match those of the grid you generated. You also need to match the resolving power, and wavelength edges you picked when setting up the grid.

Choosing the setup of your retrieval
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under ``[MCMC_params]`` you will see the following options::

    [MCMC_params]

    # whether to run the MCMC in parallel
    parallel = true
    # number of cpu to run the MCMC on
    ncpu = 30

    # number of MCMC steps
    nsteps=5000
    # fraction of the chains to be discarded as burn-in [0-1]
    frac_burnin = 0.6

    ## Which parameters to fit

    # whether to fit heterogeneity fractions
    fitspot = true
    fitfac = true

    # whether to fit temperatures of spectral components
    fitThet = true
    fitTphot = true

    # whether to fit log g of spectral components
    fitlogg_phot = true
    fitlogg_het = true

    # whether to marginalize over/fit the bare-rock transit depth
    fitDscale = true

Here is how to choose your setup for each of these parameters:

* ``parallel``: if set to ``true``, then the MCMC will be run in parallel on a number of CPUs specified by the ``ncpu`` parameter right below (by default, 30)
* ``ncpu``: Number of CPUs to use for the parallel MCMC run.
* ``nsteps``: the number of steps for each of the MCMC chains. I recommend at least 5000.
* ``frac_burnin``: the fraction of steps discarded as burn-in to obtain the posterior. By default, set to 60% (value of 0.6).
* ``fitspot``: ``true`` if you want to fit for the fraction of unocculted spots, ``false`` otherwise.
* ``fitfac``: ``true`` if you want to fit for the fraction of unocculted faculae, ``false`` otherwise.
* ``fitThet``: ``true`` if you want to fit for the temperature of unocculted spots and/or faculae, ``false`` otherwise.
* ``fitTphot``: ``true`` if you want to fit for the temperature of the photosphere, ``false`` otherwise.
* ``fitlogg_phot``: ``true`` if you want to fit the photosphere log g, ``false`` otherwise.
* ``fitlogg_het``: ``true`` if you want to fit a different log g for the spectrum of the heterogeneity component compared to that of the photosphere, ``false`` otherwise.
* ``fitDscale``: ``true`` if you want to fit for the bare-rock transit depth (recommended), ``false`` otherwise.

Priors
^^^^^^

Under ``[priors]``, the ``.toml`` file should look like this::

    [priors]

    # list of parameters with Gaussian priors. For multiple params add to the list: e.g. ["Tphot","ffac"]. Otherwise leave empty.
    gaussparanames = ["Tphot"]
    # mean and std of the Gaussian prior. For multiple parameters create a list of lists with [[mean1, std1],[mean2,std2]]. Length must match gaussparanames.
    # leave empty if no gaussparanames.
    hyperp_gausspriors = [[2566,70]]

    # specify whether we want to fit fspot/ffac with prior uniform in log or lin space
    # (e.g. [0,0]: both in lin space; [1,0]: fspot in log space, ffac in lin space)
    fitLogfSpotFac = [0,0]
    # lower and upper bound of the log(prior) on the heterogeneity fraction(s).
    hyperp_logpriors = [-5, 0]

First, you can set a Gaussian prior on any of your fitted parameters, using the ``gaussparanames`` and ``hyperp_gausspriors`` variables.

By default (uniform priors on all fitted parameters)::

    gaussparanames = ""
    hyperp_gausspriors = ""

Otherwise, you can add the name of the parameter(s) for which you want to use a Gaussian prior to ``gaussparanames``, and add a component to ``hyperp_gausspriors`` that specifies the mean and standard deviation of the gaussian parameter to adopt (looks like ``[[mean,std]]`` or `[[mean1,std1],[mean2,std2]]`` for multiple paramters in ``gaussparanames``).
Here's an example when using a Gaussian prior on the photosphere temperature (recommended, since it is not constrained by the TLSE)::

    gaussparanames = Tphot
    hyperp_gausspriors = 2566_70

The spot/faculae covering fractions can also be fitted with priors that are uniform in linear space (default) or in log space. This is dictated by the ``fitLogfSpotFac`` parameter.
* Use ``fitLogfSpotFac = [0,0]`` for the default settings of both parameters fitted with linear-uniform priors
* Set the first/second element to 1 instead to use a log-uniform prior on ``fspot``/``ffac``.
* If you choose to fit either parameter in log space, the boundaries of the prior on log(fhet) will be set by ``hyperp_logpriors = [lowerBound,upperBound]``.

If you wish to change the way the prior is set up on any of the fitted parameters, you can do it by changing the dictionary created by the function ``get_param_priors()`` in ``stellar_retrieval_utilities.py``.

Plotting
^^^^^^^^

I am providing some flexibility on how your output plots will look under ``[plotting]``::

    [plotting]

    # amount of padding in microns (unit used for spectrum plots)
    pad = 0.25

    # resolving power to smooth model spectra to (when plotting them)
    target_resP = 100

* The ``pad`` parameter roughly regulates the padding in microns added to the left and right of the spectrum plots compared to the extent of the observed spectrum
* ``target_resP`` specifies the resolving power at which you wish your stellar contamination spectra to be plotted.

Post-processing
---------------

By default, the code will produce (and save to the results folder):

Inputs to the code:

1. Input records:

* a copy of the run file that was used and of the .toml file with the specified inputs
* a copy of the version of ``stellar_retrieval_utilities.py`` that was used
* a figure displaying the spectrum being fitted
* ``defaultparams``: CSV file with the default parameters used to initialize the fit

Outputs of the code:

2. CSV files:

* ``pandas`` file: fitted parameters from the chain, with the associated log likelihood and log probability values
* ``bestfit`` file: for each parameter, the best-fit value (maximum likelihood), the max-probability values, as well as percentiles which can be used for quoting in tables
* ``bestfit_stats`` file: model comparison statistics: index of the best-fit model (in the post-burnin samples), the corresponding (reduced) chi-squared value, and BIC
* ``fixedR_1_2_3_sigma`` file: a csv file containing a set of models at the resolving power ``target_resP`` (R=100 by default) corresponding to the max-likelihood, max-probability samples, and percentiles
* ``blobs_1_2_3_sigma`` file: a csv file containing a set of models integrated within the bins of the observed spectrum corresponding to the max-likelihood, max-probability samples, and percentiles

3. NPY file: contains the "blobs": the series of models computed by the MCMC.

4. Diagnostics figures:

* ``chainplot``: chain plots, with and without the burn-in steps
* ``bestfit_model`` file: a plot of the best-fit model, integrated to match the bins in the observed spectrum, with the best-fit parameter values quoted

5. Publication-ready figures:

* ``1_2_3_sigma_withamplitude`` file: same as ``1_2_3_sigma`` but with a lower panel showing the amplitude of the stellar contamination signature across wavelength in the spectrum (in absolute terms)
* ``resP..._1_2_3_sigma`` files: fitted spectrum with the results of the fit (max-likelihood, max-probability samples, and +/- 1, 2, 3 sigma), with stellar models at higher resolution (resolving power ``target_resP``), with a log or lin scale for the wavelength axis.
* ``1_2_3_sigma`` files: fitted spectrum with the results of the fit (max-likelihood, max-probability samples, and +/- 1, 2, 3 sigma), with stellar models all integrated within the same bins as the data, with a log or lin scale for the wavelength axis.
* a corner plot of post-burnin samples

Please let me know if other things would be useful for you to have as default outputs, or feel free to create pull requests with your nice additions!