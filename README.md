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
The dependencies of *smint* are 
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
* *pymsg* (my personal favorite - needed to run create_fixedR_grid_pymsg_template.py)

### Example

#### MCMC fits

#### Post-processing
