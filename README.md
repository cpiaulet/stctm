# stctm
*stctm* (STellar ConTamination Modeling) performs spectral retrievals on exoplanet transmission spectra assuming that the spectrum can be explained by stellar contamination alone.


If you use this code, please cite Caroline Piaulet-Ghorayeb.
First public release: Piaulet-Ghorayeb et al., submitted.

Previous uses of the code (and for now, papers to cite as its best description):
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
    $ export MESASDK_ROOT=~/mesasdk
or in the code itself:
    import os
    os.environ['MESASDK_ROOT'] = "/home/caroline/mesasdk"

2. Choose your stellar parameters.
You will need to edit the star effective temperature, Fe/H, and log g. The cleanest way to do this is to add another code block corresponding to the name of your star

3. Choose the stellar grid you already downloaded.
In my case, I downloaded the ```'sg-Goettingen-HiRes.h5'```, but you can edit this to match the grid of your choice from the sample available at http://user.astro.wisc.edu/~townsend/static.php?ref=msg-grids.

4. Edit the grid model parameters to match your needs
The template I provide sets a range of log g values (```logg_range```), stellar effective temperature values (```Teff_range```), and a grid spacing (defined by ```loggstep``` and ```Teffstep```) that matches the default settings of the main code. You can however change these depending on your needs for the specific star-planet case. ** If you edit these, please pay attention to point X in the retrieval run instructions. **

I also compute the grid at a resolving power of 10,000 (```resPower_target```), and over a wavelength range from 0.2 to 5.4 microns (```wv_min_um``` and ```wv_max_um```), which you can also change to fit your needs.

To calculate a grid of models, navigate to the folder where the run script resides, and simply run:
    python create_fixedR_grid_pymsg_template.py

### Setting up a retrieval
Copy the contents of ```stctm/example/``` wherever in your installation you want to run the code.

In the way it is currently released, ```stctm``` requires you to enter the planet and star information directly in the run script.

The path to which files are saved depends on the name of the folder in which the main run file (in the example, ```stellar_retrieval_v10_generic_runfile.py``), resides.

#### MCMC fits

#### Post-processing
