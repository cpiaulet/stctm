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

class StellSpec(Table):  
    def __init__(self,inputtable,inputtype='basic',label="Stellar spectrum",color="k",
                 waveunit="um", *args, **kwargs):
        """
        Object containing the transmission spectrum

        Parameters
        ----------
        inputtable : str
            Spectrum table (astropy table object)
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
        pyStellSpec object

        """

        
        if inputtype=='basic':
            super(StellSpec,self).__init__(inputtable)
            for c in self.colnames:
                setattr(self, c, self[c])
            self.waveMin = centersToEdges(self.wave)[:-1]
            self.waveMax = centersToEdges(self.wave)[1:]

        elif inputtype=='MR_csv':
            inputtable["yval"] = inputtable["spec"]*1e10
            inputtable["waveMin"] = inputtable["wave_low"]
            inputtable["waveMax"] = inputtable["wave_high"]
            inputtable["yerrLow"] = inputtable["err"]*1e10
            inputtable["yerrUpp"] = inputtable["err"]*1e10
            inputtable["wave"] = 0.5*(inputtable["waveMin"] + inputtable["waveMax"])
            super(StellSpec,self).__init__(inputtable)
            for c in self.colnames:
                setattr(self, c, self[c])
     
        self.label=label
        self.waveunit=waveunit


    def remDataByIndex(self,ind):
        self.wave = np.delete(self.wave, ind)
        self.yval = np.delete(self.yval, ind)
        self.yerrLow = np.delete(self.yerrLow, ind)
        self.yerrUpp = np.delete(self.yerrUpp, ind)
        self.waveMin = np.delete(self.waveMin, ind)
        self.waveMax = np.delete(self.waveMax, ind)

          


        
    def plot(self,ax=None,title=None,label=None,
             xscale='linear',figsize=None,ylim=None,showxerr=True,xticks=None,xticklabels=None,
             markercolor='white',quantityToPlot=['yval','yerrLow','yerrUpp'],showcaps=True,color="k",ls="",
             marker="o",capsize=0, alpha=0.9,markeredgecolor="k", zorder=0, markerfacecolor="w",plotError=False,
             **kwargs):
        """
        Plot stellar spectrum

        Returns
        -------
        fig, ax

        """


        if label=='noshow':
            label=None
        elif label is None:
            label=self.label

                        
        if ax is None:
            fig,ax = plt.subplots(figsize=figsize)
            newFig=True
        else:
            newFig=False
        
        if xscale=='log':
            sru.xspeclog(ax,level=1)
        
        
            
        
        
       
        xval=self.wave
        yval=self.yval
        
        xerr=np.abs(np.c_[self.waveMin,self.waveMax].T - xval)
        yerr = np.c_[self.yerrLow,self.yerrUpp].T
        
        # ax.errorbar(xval,yval,xerr=xerr,yerr=yerr,color=color,ls=ls,
        #     marker=marker,capsize=capsize, alpha=alpha,markeredgecolor=markeredgecolor, 
        #     markerfacecolor=markerfacecolor,label=label,**kwargs)
        
        ax.plot(xval,yval, color=color,marker=marker,label=label, alpha=alpha,zorder=zorder,ls=ls,**kwargs)
        if plotError:
            ax.fill_between(xval,yval-self.yerrLow,yval+self.yerrUpp,color="gray", zorder=-9)

        if newFig:
            ax.set_xlabel(r"Wavelength [$\mu$m]")
            ax.set_ylabel(r'Stellar flux [$\times$ 10$^{-10}$ erg/s/cm$^2$/$\mu$m]')

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

#%%

def centersToEdges(centers):
    edges=0.5*(centers[0:-1]+centers[1:])
    edges=np.r_[ edges[0]-(edges[1]-edges[0]), edges , edges[-1]+(edges[-1]-edges[-2])   ];
    return edges
