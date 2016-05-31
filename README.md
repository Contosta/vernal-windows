# vernal-windows
This repository includes the algorithms used to identify the various thresholds that mark the onset of spring within weather, snow, terrestrial, and aquatic systems. The full vernal window is defined as the span between the earliest to the latest threshold date. Developed as part of the research under the New Hampshire ESPCoR Ecosystems and Society grant. See companion paper by Contasta et al 2016.

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