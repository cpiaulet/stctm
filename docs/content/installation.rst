Installation instructions
=========================


I recommend installing *stctm* in a clean conda environment, as with any new package.

Example command-line commands could look like this::

    conda create -n mystctm_env python=3.10.4
    conda activate mystctm_env

You can clone *stctm* directly from GitHub and place it in your installation at the path of your choice::

    git clone https://github.com/cpiaulet/stctm.git
    cd stctm
    pip install -e .

Dependencies
------------

The dependencies of *stctm* are: *NumPy*, *scipy*, *emcee*, *corner*, *astropy*, *h5py*, *matplotlib*, *pandas*, *pysynphot*, and *tqdm* (for progress bar with MCMC retrievals).

Stellar model package installation
----------------------------------

You may also need the additional dependency *pymsg* (my personal favorite).
To install *pymsg*, you can find instructions at `MSG documentation <https://msg.readthedocs.io/en/stable/>`_, and then download the grid(s) of your choice from `MSG grids <http://user.astro.wisc.edu/~townsend/static.php?ref=msg-grids>`_.

Required inputs: the stellar models grid
----------------------------------------

The code needs as input a grid of stellar models as a function of wavelength, with a regular spacing in log g and effective temperature.
I included in the repository an auxiliary file that allows you to create such a model grid (``create_fixedR_grid_pymsg_template.py``) using *pymsg*.

To make the first setup and testing of the installation easier, I provide on Zenodo the stellar models grid I generated for TRAPPIST-1. This grid is used in all the examples listed below, and contains a grid of models generated with the MSG package.
You can download it from the latest Zenodo link: `https://doi.org/10.5281/zenodo.15334399 <https://doi.org/10.5281/zenodo.15334399>`_ (``TRAPPIST_1_pymsg.h5``).

To run all the tests and examples in this documentation, please copy this file to the folder ``example/R10000_model_grids``. If you want to reproduce without any edits to the code all the tests and examples mentioned below, you will need to save this file with a relative path of ``../../R10000_model_grids/TRAPPIST_1_pymsg.h5`` relative to where the example run files for *stctm* and *exotune* are located.

To generate your own grid of interpolated models using MSG for any star of your choosing, following the instructions under `Create your own grid of stellar models using MSG <#create-your-own-grid-of-stellar-models-using-msg>`_.

Testing your installation
-------------------------

I created a dummy spectrum (with only 1 point) so you can run a few-second test of the code on your laptop for both serial and parallel runs.

#. Copy the contents of the ``example/`` directory wherever in your installation you want to run the code.
#. Make sure you have a suitable stellar models grid at the path recommended above.
#. Confirm that your environment paths are set up properly. Specifically, you need to have the ``CRDS_SERVER_URL``, ``CRDS_PATH``, and ``PYSYN_CDBS`` environment variables defined.
   You can do this via the command-line (see example below for ``CRDS_PATH``)::

     export CRDS_PATH=/home/yourusername/crds_cache

   or in the code of the analysis file itself::

     import os
     os.environ['CRDS_PATH'] = "/home/yourusername/crds_cache"

#. Navigate to your copy of ``stellar_contamination_analysis/template_analysis/`` and run the following test in the command line::

     stctm_TLS test_ini_stctm.toml

This should display print statements as the code is running, and create a directory with output files for this "mock" fit under ``stellar_contamination_results/``.

To test instead that the parallel version of the code works (with my custom wrapper around multiprocessing ``Pool``), you can simply run::

     stctm_TLS test_ini_stctm.toml -parallel=True -ncpu=2 -res_suffix=singlebin_testcode_parallel

... and that's it! You can read more in the following pages on how to customize what you do with *stctm*.