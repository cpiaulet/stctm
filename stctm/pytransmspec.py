"""
Created on Thu Jun 13, 2024

@author: cpiaulet
Definition of the TransSpec object
"""

import numpy as np
import matplotlib.pyplot as plt
import pdb
import stctm.stellar_retrieval_utilities as sru

from astropy.table import Table

from copy import deepcopy

import warnings


## Object definition

class TransSpec(Table):  
    def __init__(self,inputpath,inputtype='basic',label="Transmission spectrum",color="k",
                 waveunit="um",header_start=0):
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
        header_start : int, optional
            Starting line for header in input file if `inputtype` is set to 'wavemicrons'.

        Returns
        -------
        TransSpec object

        """

        
        if inputtype=='basic':
            table=Table.read(inputpath,format='ascii.ecsv',delimiter=',')

            super(TransSpec,self).__init__(table)

            self.sort(keys=["waveMin"])

            self['wave'].unit=waveunit
            self['waveMin'].unit=waveunit
            self['waveMax'].unit=waveunit
        
            #Metadata
            self.meta['waveunit']=waveunit
            self.meta['label']=label
            self.meta['color']=color

        elif inputtype=='wavemm':
            table=Table.read(inputpath,format='ascii.ecsv',header_start=header_start)
            
            super(TransSpec,self).__init__(table)
            self["wave"] = self["wave"]*1e3 # convert to microns from millimeters
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

    def binning(self,binFac=0,resPower=0,iwave=None):
        """
        adapted from auxbenneke/utilities.py binning function
        Bins the transmission spectrum.

        Parameters
        ----------
        binFac : int, optional
            If >0, bins by averaging every binFac points. Not used if `iwave` provided.
            resPower : float, optional
            Target spectral resolution power. If >0, spec is adaprively binned
            to achieve approximately this resolution
        iwave : list of float, optional
            Explicit wavelength bin edges. If provided, overrides binFac.
            Each bin is defined between iwave[i] and iwave[i+1].

        Returns
        -------
        spec : TransSpec
            A new TransSpec object containing the binned spectrum.
        """
        
        if iwave is not None:
        
            iwave=[0]+iwave
                
            
            specInput=deepcopy(self)
            spec=deepcopy(self[:0])
            
            for i in np.arange(len(iwave)-1):
                
                iwaveMin=iwave[i]
                iwaveMax=iwave[i+1]

                rows = specInput[np.logical_and(specInput['iwave']>=iwaveMin,specInput['iwave']<iwaveMax)]

                row=rows[0]
                row['waveMax']  = rows['waveMax'][-1]
                row['wave']     = 0.5*(row['waveMax']+row['waveMin'])
                row['resPower'] = 0.5*(row['waveMax']+row['waveMin']) / (row['waveMax']-row['waveMin'])                 
                row['yerrLow']  = 1.0 / len(rows) * np.sqrt( np.sum(rows['yerrLow']**2) )
                row['yerrUpp']  = 1.0 / len(rows) * np.sqrt( np.sum(rows['yerrUpp']**2) )
                row['yval']     = np.mean(rows['yval'])

                spec.add_row(row)


        
        elif binFac>0:
            spec=deepcopy(self)

            nBins=int(len(spec)/binFac)
            spec.remove_rows(np.arange(nBins*binFac,len(spec)))
    
            waveMin=spec['waveMin'].data[:nBins*binFac].reshape([-1,binFac])[:,0]
            waveMax=spec['waveMax'].data[:nBins*binFac].reshape([-1,binFac])[:,-1]
            yval=np.mean( spec['yval'].data.reshape([-1,binFac]) , axis=1)
            yerrLow= 1.0/binFac * np.sqrt(np.sum( spec['yerrLow'].data.reshape([-1,binFac])**2 ,axis=1))
            yerrUpp= 1.0/binFac * np.sqrt(np.sum( spec['yerrUpp'].data.reshape([-1,binFac])**2 ,axis=1))
    
            spec.remove_rows(np.arange(nBins,len(spec)))
            spec['waveMin']=waveMin; spec['waveMin'].unit='um'
            spec['waveMax']=waveMax; spec['waveMin'].unit='um'
            spec['wave']   =0.5*(spec['waveMin']+spec['waveMax'])
            spec['resPower'] = 0.5*(spec['waveMax']+spec['waveMin']) / (spec['waveMax']-spec['waveMin']); spec['waveMin'].unit='um'               
            spec['yval']   =yval
            spec['yerrLow']=yerrLow
            spec['yerrUpp']=yerrUpp
            
            waveunit = self.meta["waveunit"]
            spec['wave'].unit=waveunit
            spec['waveMin'].unit=waveunit
            spec['waveMax'].unit=waveunit
            #Meta data
            spec.meta['waveunit']=waveunit


        elif resPower>0:
            
            specbefore=deepcopy(self)
            spec=deepcopy(self[:0])

            binFactors=np.round(self['resPower']/resPower).astype(int)
            
            ilow=0
            while ilow+binFactors[ilow]-1 < len(specbefore):
                
                iupp = ilow + binFactors[ilow]-1
                
                row=specbefore[ilow]
                row['waveMax']  = specbefore['waveMax'][iupp]
                row['wave']     = 0.5*(row['waveMax']+row['waveMin'])
                row['resPower'] = 0.5*(row['waveMax']+row['waveMin']) / (row['waveMax']-row['waveMin'])                 
                row['yerrLow']  = 1.0 / (iupp-ilow+1) * np.sqrt( np.sum(specbefore['yerrLow'][ilow:iupp+1]**2) )
                row['yerrUpp']  = 1.0 / (iupp-ilow+1) * np.sqrt( np.sum(specbefore['yerrUpp'][ilow:iupp+1]**2) )
                row['yval']     = np.mean(specbefore['yval'][ilow:iupp+1])

                spec.add_row(row)
                print(row)
                print(ilow)

                ilow=iupp+1                
                if ilow>len(specbefore)-1:
                    break

        else:
            warnings.warn('Binning not successful. Specify either binFac or resPower!')

        return spec            


        
    def plot(self,ax=None,title=None,label=None,
             xscale='linear',figsize=None,ylim=None,showxerr=True,xticks=None,xticklabels=None,
             color="k",ls="",
             marker="o",capsize=0, alpha=0.9,markeredgecolor="k", markerfacecolor="w",
             **kwargs):
        """
        Plot transmission spectrum

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes object to plot on. If None, a new figure is created. Default is None.
        title : str, optional
            Title of the plot. Default is None.
        label : str, optional
            Label for the data series (used in legend). Defaults to object's meta label if None.
        xscale : str, optional
            X-axis scale: 'linear' or 'log'. Default is 'linear'.
        figsize : tuple, optional
            Size of the figure if a new one is created. Default is None.
        ylim : tuple, optional
            Y-axis limits. Default is None.
        showxerr : bool, optional
            Whether to display horizontal error bars. Default is True.
        xticks : list, optional
            Custom x-tick locations. Default is None.
        xticklabels : list, optional
            Custom x-tick labels. Default is None.
        color : str, optional
            Color of the plot line and markers. Default is "k" (black).
        ls : str, optional
            Line style. Default is "" (no line).
        marker : str, optional
            Marker symbol. Default is "o" (circle).
        capsize : float, optional
            Size of the caps on error bars. Default is 0.
        alpha : float, optional
            Transparency level. Default is 0.9.
        markeredgecolor : str, optional
            Edge color of markers. Default is "k" (black).
        markerfacecolor : str, optional
            Fill color of markers. Default is "w" (white).
        **kwargs : dict
            Additional keyword arguments passed to `ax.errorbar`.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object (if created).
        ax : matplotlib.axes.Axes
            The axes object the spectrum is plotted on.

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

#%%

def centersToEdges(centers):
    """
    from auxbenneke/utilities.py: convert bin centers to bin edges
    Parameters
    ----------
    centers : array-like
        List of bin centers


    Returns
    -------
    edges : array
        List of bin edges corresponding to the input bin centers.
        The first and last edges are extrapolated by half the bin width.
    """
    edges=0.5*(centers[0:-1]+centers[1:])
    edges=np.r_[ edges[0]-(edges[1]-edges[0]), edges , edges[-1]+(edges[-1]-edges[-2])   ];
    return edges
