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
import numpy as np
import matplotlib.pyplot as plt
import os

os.environ["OMP_NUM_THREADS"] = "1"

import stctm.pystellspec as psspec
import stctm.exotune_utilities as xtu
import astropy.constants as const
import pdb
import emcee
import astropy.table as table
import h5py
from copy import deepcopy
import sys
import astropy.io as aio
import shutil

import matplotlib.gridspec as gridspec
import pandas as pd

from astropy.convolution import convolve, Box1DKernel,  Gaussian1DKernel, Trapezoid1DKernel
import corner
from multiprocessing import Pool
# from emcee.utils import MPIPool

from .exotune_parallel_setup import *
    
    
    
#%% Define required functions for parallel runs using global variables



    
    
    