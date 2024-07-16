# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13, 2024

@author: cpiaulet
"""

import numpy as np
import matplotlib.pyplot as plt
import pdb
import astropy.io as aio
import stctm.stellar_retrieval_utilities as sru

import pandas as pd
import os        
from astropy.table import Table
import astropy.io.fits as pf

from scipy.optimize import minimize
from scipy import stats

import collections
from copy import deepcopy

import warnings

jdref=2450000
big=1e10


#%%

class TransSpec(Table):  
    def __init__(self,inputpath,inputtype='basic',label="Transmission spectrum",color="k",
                 waveunit="um"):
        """
        Object containing the transmission spectrum

        Parameters
        ----------
        inputpath : str
            Path to the spectrum file
        inputtype : str, optional
            Type of input spectrum. Create your own to read in a custom 
            spectrum file format. The default is 'basic'.
        label : str, optional
            Label to be used for plotting. The default is "Transmission spectrum".
        color : str, optional
            Color for plotting. The default is "k".
        waveunit : str, optional
            Unit of the wavelength axis. The default is "um".

        Returns
        -------
        pyTransSpec object

        """

        
        if inputtype=='basic':
            table=Table.read(inputpath,format='ascii.ecsv',delimiter=',')
            super(TransSpec,self).__init__(table)
            self["wave"] = self["wave"]*1e3
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


            


        
    def plot(self,ax=None,title=None,label=None,
             xscale='linear',figsize=None,ylim=None,showxerr=True,xticks=None,xticklabels=None,
             markercolor='white',quantityToPlot=['yval','yerrLow','yerrUpp'],showcaps=True,color="k",ls="",
             marker="o",capsize=0, alpha=0.9,markeredgecolor="k", markerfacecolor="w",
             **kwargs):
        """
        Plot transmission spectrum

        Returns
        -------
        fig, ax

        """


        if label=='noshow':
            label=None
        elif label is None:
            label=self.meta['label']

                        
        if ax is None:
            fig,ax = plt.subplots(figsize=figsize)
            newFig=True
        else:
            newFig=False
        
        if xscale=='log':
            sru.xspeclog(ax,level=1)
        
        
            
        
        tbl = self
       
        xval=tbl["wave"]
        yval=tbl["yval"]
        
        xerr=np.abs(np.c_[tbl['waveMin'],tbl['waveMax']].T - xval)
        yerr = np.c_[tbl['yerrLow'],tbl['yerrUpp']].T
        
        ax.errorbar(xval,yval,xerr=xerr,yerr=yerr,color=color,ls=ls,
            marker=marker,capsize=capsize, alpha=alpha,markeredgecolor=markeredgecolor, 
            markerfacecolor=markerfacecolor,label=label,**kwargs)
        
        if newFig:
            ax.set_xlabel(r"Wavelength [$\mu$m]")
            ax.set_ylabel(r'Transit Depth [ppm]')
            ax.minorticks_on()

        if ylim is not None:
            ax.set_ylim(ylim)

        if title is not None:
            ax.set_title(title)

        if xticks is not None:
            ax.set_xticks(xticks)
            
        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
            for label in ax.get_xmajorticklabels():
                label.set_rotation(90)            

        if newFig:
           return fig,ax
       
