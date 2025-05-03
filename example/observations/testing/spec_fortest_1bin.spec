# %ECSV 1.0
# ---
# datatype:
# - {name: instrname, datatype: string}
# - {name: reference, datatype: string}
# - {name: bandpass, datatype: string}
# - {name: iwave, datatype: int64}
# - {name: wave, unit: um, datatype: float64}
# - {name: waveMin, unit: um, datatype: float64}
# - {name: waveMax, unit: um, datatype: float64}
# - {name: xMin, datatype: float64}
# - {name: xMax, datatype: float64}
# - {name: yval, unit: ppm, datatype: float64}
# - {name: yerrLow, unit: ppm, datatype: float64}
# - {name: yerrUpp, unit: ppm, datatype: float64}
# - {name: wlcLow, datatype: float64}
# - {name: wlcUpp, datatype: float64}
# - {name: ignore, datatype: int64}
# - {name: referenceLink, datatype: string}
# delimiter: ','
# meta: !!omap
# - {waveunit: um}
# - {name: null}
# - {label: null}
# - {spectype: dppm}
# - {ylabel: 'Transit Depth [ppm]'}
# - sysModel: {addParas: null, bestfit: null, guess: null, sysModel: null}
# schema: astropy-2.0
instrname,reference,bandpass,iwave,wave,waveMin,waveMax,xMin,xMax,yval,yerrLow,yerrUpp,wlcLow,wlcUpp,ignore,referenceLink
JWST_NIRISS_SOSS,,uniform,0,1.7047959392890335,0.6,2.809591878578067,0.0,0.0,6767.5935,32.2496,32.2496,0.0,0.0,0,none
