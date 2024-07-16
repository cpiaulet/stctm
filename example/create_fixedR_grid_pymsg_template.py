#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 12:50:31 2023

@author: caroline
Create a grid of stellar models for stellar contamination analyses using pysynphot
"""
import numpy as np
import matplotlib.pyplot as plt
import os
os.environ['CRDS_SERVER_URL'] = "https://jwst-crds.stsci.edu"
os.environ['CRDS_PATH'] = "/home/caroline/crds_cache"
os.environ['PYSYN_CDBS'] = "/home/caroline/trds"


os.environ['MESASDK_ROOT'] = "/home/caroline/mesasdk"
os.environ['MSG_DIR'] = "/home/caroline/msg-1.2"

import pysynphot as S
import pymsg as msg
import astropy.constants as const
import astropy.units as u
import pdb
import h5py

import stellar_retrieval_utilities as sru 

#%% Stellar parameters

star_name = "TRAPPIST-1"
resPower_target = 10000

if 1:
    # for TRAPPIST-1
    param = dict()
    param["Tphot"] = 2566. #K, Agol+20
    param["met"] = 0.040 #[Fe/H], Gillon+2017
    param["logg_phot"] = 5.2 # cgs, Agol+20
    


use_pymsg = True

#%%
print("Loading pymsg SpecGrid...")
MSG_DIR = os.environ['MSG_DIR']
GRID_DIR = os.path.join(MSG_DIR, 'data', 'grids')

specgrid_file_name = os.path.join(GRID_DIR, 'sg-Goettingen-HiRes.h5')

specgrid = msg.SpecGrid(specgrid_file_name)

# Inspect grid parameters

print('Grid parameters:')

for label in specgrid.axis_labels:
    print(f'  {label} [{specgrid.axis_x_min[label]} -> {specgrid.axis_x_max[label]}]')

print(f'  lam [{specgrid.lam_min} -> {specgrid.lam_max}]')


#%% Generate a grid of PHOENIX models over the prior range for the star
# range of params for the grid

logg_range = [2.5,5.5] 
Teff_range = [np.min([2300.-param["Tphot"], -100.])+param["Tphot"], param["Tphot"]+1000.] 
loggstep = 0.1 #cgs
Teffstep = 20. #K
resPower_target = 10000
wv_min_um = 0.2
wv_max_um = 5.4

# stellar models grid file path
if star_name == "TRAPPIST_1":
    stmodfile = "R"+str(resPower_target)+"_model_grids/"+star_name+"_pymsg.h5" 
 
        # StellarModelGridFilePath = self.scarletpath+'/'+'../other_models/StellarModels/CustomGrids/'+os.path.basename(self.StellarModelGridFile)+'_pymsg_withTgrid_{:g}_logg{:g}_{:g}_{:g}_{:g}_{:g}_{:g}.h5'.format(params['Teffstar'],logg_range[0], logg_range[1],params['metfe'],waveRange[0],waveRange[1],resolution) 

print("Grid parameter ranges:")
print("logg:", logg_range)
print("Teff:", Teff_range)

wv_template_thisR, wv_template_edges_thisR = sru.make_wavelength_array(wv_min_um=wv_min_um, wv_max_um=wv_max_um, resPower=resPower_target, use_pymsg=use_pymsg)

# check if the grid file exists
if os.path.exists(stmodfile):
    print("Stellar model grid file existed ! Reading it in...")
    Teffs_grid = np.arange(Teff_range[0], Teff_range[1], Teffstep)
    loggs_grid = np.arange(logg_range[0], logg_range[1], loggstep)
    # loggs_grid = np.array([logg])
    
    h5f = h5py.File(stmodfile, 'r')

    models_grid_fixedR = h5f['model_grid'][:]
    h5f.close()    
    print("Grid loaded, file closed.")
# if not, generate it
else:
    print("Generating grid of fixed-resolution models...") 
    Teffs_grid, loggs_grid, models_grid_fixedR = sru.generate_PHOENIX_grid_fixedR(wv_template_thisR, wv_template_edges_thisR, feh=param["met"], 
                                                                                  Teff_range=Teff_range, logg_range=logg_range, 
                                                                                  Teffstep=Teffstep, loggstep=loggstep, fname_save=stmodfile, 
                                                                                  use_pymsg = True, pymsg_specgrid=specgrid)
    

