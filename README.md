# stctm
*stctm* (STellar ConTamination Modeling) performs spectral retrievals on exoplanet transmission spectra assuming that the spectrum can be explained by stellar contamination alone.


If you use this code, please cite the associated Zenodo repository, with the following DOI: 10.5281/zenodo.13153251 (see Citation section below). 
If you use it, please make sure to also cite MSG (https://doi.org/10.21105/joss.04) and the source of the stellar models grid you use.
First public release associated to: Piaulet-Ghorayeb et al., submitted.

Previous uses of the code:
* Lim et al., 2023 (https://ui.adsabs.harvard.edu/abs/2023ApJ...955L..22L/abstract), on TRAPPIST-1 b
* Roy et al., 2023 (https://ui.adsabs.harvard.edu/abs/2023ApJ...954L..52R/abstract), on GJ 9827 d
* Benneke et al., submitted, on TRAPPIST-1 g

Bare-bones skeleton of the code released for now - full user instructions and user-friendly setup to come soon!

I provide as an example for the application of *stctm* to Visit 1 of TRAPPIST-1 b (Lim et al., 2023, https://ui.adsabs.harvard.edu/abs/2023ApJ...955L..22L/abstract), fitting for two populations of heterogeneities: spots, faculae, and allowing the log g of the heterogeneities to vary relative to that of the star. 

## Installation

You can install *stctm* from GitHub:

    git clone https://github.com/cpiaulet/stctm.git
    cd stctm
    pip install -e .
    
### Dependencies
The dependencies of *stctm* are 
* *NumPy*
* *scipy*
* *emcee*
* *corner*
* *astropy*
* *h5py*
* *matplotlib*
* *pandas*
* *tqdm*

#### Stellar models

You may also need 
* *pysynphot* (if you want to use their version of the stellar models interpolator), and/or
* *pymsg* (my personal favorite - needed to run ```create_fixedR_grid_pymsg_template.py```)

To install *pymsg*, you can find instructions at https://msg.readthedocs.io/en/stable/ and then download the grid(s) of your choice from http://user.astro.wisc.edu/~townsend/static.php?ref=msg-grids.

#### Create your own grid of stellar models using MSG

If you choose to use the *MSG* module for stellar models, the code requires a pre-computed grid of stellar models for the planet of interest.
I provide a template code snippet for how to go about computing this stellar models grid in ```create_fixedR_grid_pymsg_template.py```. Here are a few things to pay attention to.

1. Make sure that your paths are set up properly.
Specifically, you need to have the ```MESASDK_ROOT``` and ```MSG_DIR``` environment variables defined.
You can do this via the command-line:

    export MESASDK_ROOT=~/mesasdk
   
or in the code itself:

    import os
    os.environ['MESASDK_ROOT'] = "/home/caroline/mesasdk"

3. Choose your stellar parameters.
You will need to edit the star effective temperature, Fe/H, and log g. The cleanest way to do this is to add another code block corresponding to the name of your star

4. Choose the stellar grid you already downloaded.
In my case, I downloaded the ```'sg-Goettingen-HiRes.h5'```, but you can edit this to match the grid of your choice from the sample available at http://user.astro.wisc.edu/~townsend/static.php?ref=msg-grids.

5. Edit the grid model parameters to match your needs
The template I provide sets a range of log g values (```logg_range```), stellar effective temperature values (```Teff_range```), and a grid spacing (defined by ```loggstep``` and ```Teffstep```) that matches the default settings of the main code. You can however change these depending on your needs for the specific star-planet case. If you edit these, make sure to pay attention to the section "Setting up the stellar parameters and reading in the grid of stellar models" in the retrieval run instructions!

I also compute the grid at a resolving power of 10,000 (```resPower_target```), and over a wavelength range from 0.2 to 5.4 microns (```wv_min_um``` and ```wv_max_um```), which you can also change to fit your needs.

To calculate a grid of models, navigate to the folder where the run script resides, and simply run:
    python create_fixedR_grid_pymsg_template.py

### Setting up a retrieval
Copy the contents of ```stctm/example/``` wherever in your installation you want to run the code.

In the way it is currently released, ```stctm``` requires you to enter the planet and star information directly in the run script. I am working on a more user-friendly way of doing this using a setup file which should be added to the main branch soon.

#### File paths
The path to which files are saved does not depend on the name of the folder in which the main run file (in the example, ```stellar_retrieval_v10_generic_runfile.py``), resides. However, relative paths matter as the results folder will reside in ```../../stellar_contamination_results/``` (it is automatically created by the code). 

#### Setting up labels and path to spectrum file
Under ```User inputs: Spectrum and labels``` you can set up:
* ```label``` (used for plotting)
* ```path_to_spec``` (path to your spectrum file) as well as ```spec_format``` (your spectrum is read in using the ```spec_format``` setting you choose - if you are not sure or need to add another option to read in your specific format, you can do so in ```pytransspec.py```!)
* ```res_suffix```: a suffix used for all the files that will be saved as a result of this run, in the results folder. This is the identifier you can use to record information on the spectrum, the setup of the fit, etc.

#### Choosing the setup of your retrieval

Under ```User inputs: MCMC fitting params``` you can choose:
* ```nsteps```: the number of steps for each of the MCMC chains. I recommend at least 5000, but I chose 3000 to make the test run a bit quicker :)
* ```frac_burnin```: the fraction of steps discarded as burn-in to obtain the posterior. By default, set to 60%.
* ```fitspot```: ```True``` if you want to fit for the fraction of unocculted spots, ```False``` otherwise.
* ```fitfac```: ```True``` if you want to fit for the fraction of unocculted faculae, ```False``` otherwise.
* ```fitThet```: ```True``` if you want to fit for the temperature of unocculted spots and/or faculae, ```False``` otherwise.
* ```fitTphot```: ```True``` if you want to fit for the temperature of the photosphere, ```False``` otherwise.
* ```fitDscale```: ```True``` if you want to fit for the bare-rock transit depth (recommended), ```False``` otherwise.
* ```fitlogg_het```: ```True``` if you want to fit a different log g for the spectrum of the heterogeneity component compared to that of the photosphere, ```False``` otherwise.
* ```save_fit```: ```True``` to save files to the results directory during the post-processing steps

Priors on the stellar parameters:
* ```Teffstar_prior_std```: in Kelvin, standard deviation of the Gaussian prior adopted on the stellar effective temperature;
* ```logg_phot_source```: ```value``` to use the value of ```logg_phot_value``` as the stellar photosphere log g, otherwise ```loggstar``` to use the value provided in the code block below containing the stellar parameters;
* ```logg_het_default_source```: ```value``` to use the value of ```logg_het_value``` as the heterogeneities (default, if fitted) log g, otherwise ```logg_phot``` to set it to the same value as the stellar photosphere log g.

#### Setting up the stellar parameters and reading in the grid of stellar models

Under ```User inputs: Stellar parameters```, enter the parameters of the star to set the defaults for the fit. Make sure that you have an option ("if" statement) matching the ```which_star``` setting you entered above.

Under ```User inputs: Read in stellar models grid```, modify the range and spacing of the grid in the log g and Teff dimensions to match those of the grid you generated. You also need to match the resolving power, and wavelength edges you picked when setting up the grid.

#### Post-processing

By default, the code will produce (and save to the results folder):
* a copy of the run file that was used
* a copy of the version of ```stellar_retrieval_utilities.py``` that was used
* chain plots, with and without the burn-in steps
* a csv file containing the fitted parameters from the chain, with the associated log likelihood and log probability values
* a csv file containing, for each parameter, the best-fit value (maximum likelihood), the max-probability values, as well as percentiles which can be used for quoting in tables
* a set of models at the grid resolution (```target_resP```) also corresponding to the max-likelihood, max-probability samples, and percentiles
* figures showing the spectrum with the results of the fit (max-likelihood, max-probability samples, and +/- 1, 2, 3 sigma), with stellar models at higher resolution (resolving power of 100 in the plot) or integrated within the same bins as the data points), and a plot of the stellar contamination contribution as a function of wavelength
* a corner plot of post-burnin samples
* a plot of the best-fit model, integrated to match the bins in the observed spectrum

Please let me know if other things would be useful for you to have as default outputs, or feel free to create pull requests with your nice additions!

### Changing the prior setup

If you wish to change the way the prior is set up on any of the fitted parameters, you can do it by changing the dictionary created by the function ```get_param_priors()``` in ```stellar_retrieval_utilities.py```.

## Citation

The following entry to a bib file can be used to cite this code:
    @misc{piaulet_stctm_2024,
        author       = {Caroline Piaulet-Ghorayeb},
        title        = {{stctm: Stellar contamination retrievals and modeling for small planet transmission spectra}},
        month        = aug,
        year         = 2024,
        doi          = {10.5281/zenodo.13153251},
        version      = {1.0.0},
        publisher    = {Zenodo},
        url          = {https://doi.org/10.5281/zenodo.13153251}
        }

