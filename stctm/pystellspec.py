"""
Created on Thu Jun 13, 2024

@author: cpiaulet
"""

import numpy as np
import matplotlib.pyplot as plt
import pdb
import stctm.stellar_retrieval_utilities as sru

import pandas as pd
import os        
from astropy.table import Table



## Class definition

class StellSpec(Table):  
    def __init__(self,inputtable,inputtype='basic',label="Stellar spectrum",color="k",
                 waveunit="um"):

        """
        Initialize a StellSpec object containing stellar spectrum data.

        Parameters
        ----------
        inputtable : astropy.table.Table
            Input spectrum table.
        inputtype : str, optional
            Format of the input data. Supported types are:
            - 'basic': Assumes a column named 'wave' is present.
            - 'MR_csv': Converts custom columns to standard format.
            Default is 'basic'.
        label : str, optional
            Label for plotting and identification. Default is "Stellar spectrum".
        color : str, optional
            Color used for plotting. Default is 'k' (black).
        waveunit : str, optional
            Unit of the wavelength axis. Default is 'um'.
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

    def remDataByIndex(self, ind):
        """
        Remove entries at the specified index/indices.
        Parameters
        ----------
        ind : int or array-like
            Index or indices of the data points to remove from the spectrum.
        Notes
        -----
        This method updates all relevant attributes, including wavelength,
        flux, errors, and wavelength bounds, by removing the specified entries.
        """
        self.wave = np.delete(self.wave, ind)
        self.yval = np.delete(self.yval, ind)
        self.yerrLow = np.delete(self.yerrLow, ind)
        self.yerrUpp = np.delete(self.yerrUpp, ind)
        self.waveMin = np.delete(self.waveMin, ind)
        self.waveMax = np.delete(self.waveMax, ind)

    def plot(self, ax=None, title=None, label=None,
             xscale='linear', figsize=None, ylim=None, xticks=None, xticklabels=None, color="gray", ls="",
             marker=".", alpha=0.9, zorder=0, markersize=1, plotError=False,
             **kwargs):
        """
        Plot the stellar spectrum.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, a new figure and axes are created.
        title : str, optional
            Title of the plot.
        label : str, optional
            Label for the plotted data. If 'noshow', no label is shown.
        xscale : str, optional
            Scale of the x-axis ('linear' or 'log'). Default is 'linear'.
        figsize : tuple, optional
            Size of the figure (width, height). Only used if ax is None.
        ylim : tuple, optional
            Limits for the y-axis.
        xticks : list, optional
            Custom x-axis tick locations.
        xticklabels : list, optional
            Custom x-axis tick labels.
        color : str, optional
            Color of the plotted line/markers. Default is "k".
        ls : str, optional
            Line style. Default is "" (no line).
        marker : str, optional
            Marker style.
        alpha : float, optional
            Transparency of the plot elements.
        zorder : int, optional
            Drawing order of plot elements.
        markersize : int, optional
            Size of the markers in the plot
        plotError : bool, optional
            If True, shaded error regions are plotted.
        **kwargs : dict
            Additional keyword arguments passed to the plotting function.
        Returns
        -------
        if ax was None, returns:
            fig : matplotlib.figure.Figure
                Figure object, only if a new figure was created.
            ax : matplotlib.axes.Axes
                Axes containing the plot.
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
        
        ax.plot(xval,yval, color=color,marker=marker,label=label, alpha=alpha,zorder=zorder,ls=ls,
                markersize=markersize,**kwargs)
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

## Spectrum-specific utility function

def centersToEdges(centers):
    """
    Convert an array of bin centers to bin edges.

    Parameters
    ----------
    centers : array-like
        1D array of bin center values.

    Returns
    -------
    edges : ndarray
        1D array of bin edges with length `len(centers) + 1`.
        The first and last edges are extrapolated to maintain uniform spacing.
    """
    edges=0.5*(centers[0:-1]+centers[1:])
    edges=np.r_[ edges[0]-(edges[1]-edges[0]), edges , edges[-1]+(edges[-1]-edges[-2])   ];
    return edges
