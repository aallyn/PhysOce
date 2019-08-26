#####
## Climate SST extraction code
#####

# Setting paths -- you will need to update these
code.path<- "./Code/"
input.path<- "./Input/"
output.path<- "./Output/"

# First, need to source "clim_sst_functions.R"
source(paste(code.path, "clim_sst_functions.R", sep = ""))

# Load libraries
library_check(c("tidyverse", "raster", "ncdf4", "lubridate"))

# Load in input file, providing the information for which climate models to access from ESGF-GFDL THREDDS server. The input file should be checked against the data catalog here: https://esgdata.gfdl.noaa.gov/thredds/catalog/esgcet/catalog.html. After navigating to the website, you can manually enter results for a given run of interest and then check the httpserver link. Moving forward, we might want to look at how to filter all files based on more minimal user input (i.e., "rcp85", "tos"). For now though, this seems to be okay. 
input.params<- read_csv(paste(input.path, "ClimatemodelParams.csv", sep = ""))

# Generate THREDDS URLS by mapping the "make_climate_threddsURL.R" function to our input parameter file. Might want to add a stop for URL not found?
thredds.use<- input.params %>%
  mutate(., "Thredds.URLS" = pmap(list(project = Project, model.short = Model.Short, experiment = Experiment, frequency = Time.Frequency, model.realm1 = Model.Realm1, model.realm2 = Model.Realm2, ensemble = Ensemble, version = Version, variable = Variable), make_climate_threddsURL))

# Extract climate data by mapping the "clim_data_extract.R" function to the thredds.use Thredds.URLS generated above
box.use<- c(-77, -60, 35, 46)
clim.results<- thredds.use %>%
  mutate(., "Clim.Data" = pmap(list(thredds.urls = Thredds.URLS, box = list(box.use), project = Project, out.path = list(output.path)), clim_data_extract))

# Not exactly sure how to check if this all actually worked...
par(mfrow = c(1, 3))
plot(clim.results$Clim.Data[[1]][[1140]])
plot(clim.results$Clim.Data[[3]][[1140]])
plot(clim.results$Clim.Data[[5]][[1]])

                        