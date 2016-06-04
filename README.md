# vernal-windows
This repository includes algorithms used to identify thresholds that occur during the vernal window,  which is a period that describes the end of of winter and the start of the growing season. The vernal window starts when air temperatures rise above zero degrees celcius and ends when forests leaf out and pastures, croplands, and lawns green up.  A series of dramatic transitions occurs within this window as the system crosses rapid thermodynamic and biogeochemical thresholds that drive energy, water, solute, and carbon flows.  Contosta et al. (2016) describe the sequence of thresholds that occur within the vernal window, lags between thresholds, and shifts in the timing of these lags with climate change.  This paper is currently in review.  The research was supported by the New Hampshire EPSCoR Ecosystem and Society Project (NSF-EPS 1101245).

Files in this repository:

- <b>air_zero.R</b>: R script that identifies the date where the smoothed air temperature exceeds zero celsius
- <b>baseflow.R</b>: R script that identifies the date where the smoothed discharge reaches baseflow conditions
- <b>MODIS data extraction.pdf</b>: schematics of ArcGIS models used to extract data from the MODIS dataset for subsequent calculation of snowmelt and canopy leaf-out thresholds.
- <b>MODIS snowmelt & canopy.txt</b>: visual basic script that identifies the dates of snowmelt and canopy leaf-out, using data for New Hampshire extracted from the MODIS dataset
- <b>peak & trough.R</b>: R script that identifies the dates instream and in soils where:
  - smoothed conductivity values reach a maximum value, which occurs just prior to snowmelt
  - smoothed conductivity values reach a minimum value, which occurs at peak springtime water level
- <b>piecewise regression.R</b>: R script that identifies the date that marks the onset of warming in soils and streams, by identifying the date correspoding to the breakpoint in a piecewise regression of the two best-fit lines of the smoothed temperature time series
- <b>snowfree.R</b>: R script that identifies the date of completion of snowmelt, defined as the date where the smoothed snow water equivalent (SWE) values drop to zero.
