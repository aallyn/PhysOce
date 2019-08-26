# PhysOce
This is a collection of functions for gathering and summarizing historical, current and future physical oceanographic data -- mainly sea surface temperature data. 

# Code
## obs_sst_functions.R
This R script contains functions for gathering sea surface temperature for the observation period. The core function within this script is the "env_data_extract.R" function. This function accesses OISST (plans for MURSST, ERSST) data from the THREDDS server for a region and time period of interst and saves the downloaded data as a raster stack. 

## clim_sst_functions.R
This R script containes functions for gathering sea surface temeprature projections from climate models. There are two core functions within this script. First is the "make_climate_threddsURL.R" function. This function takes in information from an input file and generates correct THREDDS URLs (see example of required information "./Input/ClimateModelParams.csv" file). Second is the "clim_data_extract.R" function. This function uses the THREDDS URLs to access the climate projection data, crop the full extent to a specific region of interest and save the data as a raster stack. 

# To do
- Post-processing of climate data to get overall ensemble mean and percentiles
- Time series plotting function. There are the beginings of a function in the "env_data_timeseries" function in the "obs_sst_functions.R" script, though needs to be edited and checked. 
