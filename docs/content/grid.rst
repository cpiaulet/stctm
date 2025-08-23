Create a grid of stellar models using MSG
=========================================


If you choose to use the *MSG* module for stellar models, the code requires a pre-computed grid of stellar models for the planet of interest.
I provide a template code snippet for how to compute this stellar models grid in ``create_fixedR_grid_pymsg_template.py``. Here are a few things to pay attention to:

1. **Make sure that your paths are set up properly.**
   Specifically, you need to have the ``MESASDK_ROOT`` and ``MSG_DIR`` environment variables defined.
   You can do this via the command-line::

       export MESASDK_ROOT=~/mesasdk

   or in the code itself::

       import os
       os.environ['MESASDK_ROOT'] = "/home/caroline/mesasdk"

2. **Choose your stellar parameters.**
   You will need to edit directly in the script (``create_fixedR_grid_pymsg_template.py``) the star's effective temperature, Fe/H, and log g. The cleanest way to do this is to add another code block corresponding to the name of your star.

e.g.::

    elif star_name =="TRAPPIST_1":
        # for TRAPPIST-1
        param = dict()
        param["Tphot"] = 2566. #K, Agol+2021
        param["met"] = 0.040 #[Fe/H], Gillon+2017
        param["logg_phot"] = 5.2 # cgs, Agol+2021

3. **Choose the stellar grid you already downloaded.**
   In my case, I downloaded the ``'sg-Goettingen-HiRes.h5'``, but you can edit this to match the grid of your choice from the sample available at `MSG grids <http://user.astro.wisc.edu/~townsend/static.php?ref=msg-grids>`_.

The way to specify this in the code is as follows::

    specgrid_file_name = os.path.join(GRID_DIR, 'sg-Goettingen-HiRes.h5')



4. **Edit the grid model parameters to match your needs.**
   The template I provide sets a range of log g values (``logg_range``), stellar effective temperature values (``Teff_range``), and a grid spacing (defined by ``loggstep`` and ``Teffstep``) that matches the default settings of the main code. You can however change these depending on your needs for the specific star-planet case. If you edit these, make sure to pay attention to the section "Setting up the stellar parameters and reading in the grid of stellar models" in the retrieval run instructions!

   I also compute the grid at a resolving power of 10,000 (``resPower_target``), and over a wavelength range from 0.2 to 5.4 microns (``wv_min_um`` and ``wv_max_um``), which you can also change to fit your needs by editing the following lines::

    logg_range = [2.5,5.5]
    Teff_range = [np.max([param["Tphot"]-1000, 2300.]), param["Tphot"]+1000.]
    loggstep = 0.1 #cgs
    Teffstep = 20. #K
    resPower_target = 10000
    wv_min_um = 0.2
    wv_max_um = 5.4

To calculate a grid of models, navigate to the folder where the run script resides, and simply run::

    python create_fixedR_grid_pymsg_template.py