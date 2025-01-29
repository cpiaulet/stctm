#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 12:34:45 2023

@author: caroline

Stellar model fitting utility functions (parallel runs)
"""

#%%------------------- Import modules ------------------#
import os
# os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
# os.environ['CRDS_PATH'] = "/Users/caroline/crds_cache"
# os.environ['PYSYN_CDBS'] = "/Users/caroline/trds"

os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
os.environ['CRDS_PATH'] = "/home/caroline/crds_cache"
os.environ['PYSYN_CDBS'] = "/home/caroline/trds"


os.environ["OMP_NUM_THREADS"] = "1"


# from emcee.utils import MPIPool

#%%
def init():
    global xxx
    xxx = 3
def define_global_parallel(param,fitparanames,spec_pickleformat,fl_phot_spot_fac,
                           Teffs_grid,loggs_grid,wv_template_thisR,models_grid_fixedR,
                           gaussparanames,mean_Gauss_para,std_Gauss_para):
    # For lnlike
    global param0
    global fitparanames0
    global spec0
    global models_baseline0
    global T_grid0
    global logg_grid0
    global model_wv0
    global models_grid0
    global gaussparanames0
    global mean_Gauss_para0
    global std_Gauss_para0
    
    param0 = param
    fitparanames0 = fitparanames
    spec0 = spec_pickleformat
    models_baseline0 = fl_phot_spot_fac
    T_grid0 = Teffs_grid
    logg_grid0 = loggs_grid
    model_wv0 = wv_template_thisR
    models_grid0 = models_grid_fixedR

    # For get_param_priors
    gaussparanames0 = gaussparanames
    mean_Gauss_para0 = mean_Gauss_para
    std_Gauss_para0 = std_Gauss_para  