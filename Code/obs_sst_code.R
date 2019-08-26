######
## Environmental data access and visualization code
######

# Preliminaries: Shared paths, sourcing libraries and obs_sst_function ----------------------
source("https://raw.githubusercontent.com/GMRI-SEL/LabFunctionsandCode/master/GenerateSharedPathsScript.R")

# Source "library_check" helper function to install/load required libraries
source(paste(lab.func.path, "GeneralHelpers.R", sep = ""))
library_check(c("raster", "tidyverse", "sf", "here", "ncdf4"))

# Source "obs_sst_functions.R" function in this Project Folder
source(here("Code", "obs_sst_functions.R"))

# Extracting OISST data ---------------------------------------------------
oisst.dat<- obs_sst_extract(data.set = "OISST", dates = NULL, box = c(-90, -30, 30, 70), out.dir = paste(res.data.path, "OISST/", sep = ""), out.file = "NWAtl")

# What did we get?
print(oisst.dat) # Daily OISST data, 13707 days to be exact.
plot(oisst.dat[[1]])

# Visualizing OISST data
env_data_timeseries(oisst.dat, baseline = c("1982-01-01", "2011-01-01"), regions = c("NELME", "GoM", "SNE-MAB"), out.dir = "~/Dropbox/Andrew/Work/GMRI/Projects/AllData/")

