#####
## Processing climate model output with cdo terminal functions and R
#####
library(tidyverse)
library(raster)
library(zoo)
library(lubridate)
library(viridis)
library(sf)
library(ncdf4)

land.shp<- st_read("~/Box/Data/Shapefiles/ne_50m_land/ne_50m_land.shp")

### Background
# For a lot of the climate model output (and satelitte temperature data), things are fairly straightforward using a combination of the netcdf4, netcdf4.helpers and raster libraries. With those, I have functions that can subset a global netcdf file to a region and time of interest, load the data, and save it as a nice raster stack. Unfortunately, the CMIP5 projections for sea surface temperature were not playing as nice. There were a lot of bizarre grid things (for example, some going 0-360, others -280-80, and different resolutions), which was a nightmare to handle using just the functions available with R. So, after about a week of trying to hack a solution, I ended up coming across the cdo commands and figured out how to run them through R/Rstudio. 

### Preliminaries.
# We are going to use terminal and the collection of CDO commands. To be able to do this, I first had to install "homebrew" and then used "brew install" to install cdo (see here: https://code.mpimet.mpg.de/projects/cdo/wiki/MacOS_Platform). I was prompted about permission for xcode build license and followed instructions. # CDO has a number of excellent functions for dealing with climate data. Most importantly for our purposes right now is the "remap" function as each of the climate models can be on a different grid. The process for dealing with this is really straight forward. First, we create a grid. We do this using the exact same language we would in the terminal, however, we wrap it in a "system" call to execute it from RStudio.

# First, the paths to folders where we are doing all our work
cd.to.path<- list("RCP85" = "~/Box/Data/CMIP5_SST/RCP85/", "RCP45" = "~/Box/Data/CMIP5_SST/RCP45/", "Historical" = "~/Box/Data/CMIP5_SST/Historical/")

# For each of these folders, we are going to want to create a grid, then remap each of the climate model files located within the folder. 

# Some house keeping
# Region of interest
ext.x<- c(-76, -60)
ext.y<- c(35, 46)
ext.use<- as.numeric(c(ext.x, ext.y))
lonlatbox<- paste(as.character(ext.use), collapse = ",")

# Pattern to use for files
pattern.use<- "*.[[:digit:]]\\.*nc"

# Now looping through folders and over files
for(i in seq_along(cd.to.path)){
  # CD to folder and create our grid within the system call
  cd.to.path.use<- cd.to.path[[i]]
  system(paste("cd", cd.to.path.use, "&& cat > Std1degGrid <<EOF\n
               gridtype = lonlat\n
               xsize = 360\n
               ysize = 180\n
               xfirst = -179.5\n
               yfirst = -89.5\n
               yinc = 1\n
               xinc = 1\n
               EOF", sep=" "))
  
  clim.mod.files<- list.files(cd.to.path.use, pattern = pattern.use)
  clim.mods.df<- data.frame("File.Long" = clim.mod.files)
  clim.mods.df$Group<- unlist(lapply(strsplit(as.character(clim.mods.df$File.Long), "_"), "[", 3))
  
  for(j in seq_along(clim.mods.df$File.Long)){
    regrid.use<- "Std1degGrid"
    in.file<- as.character(clim.mods.df$File.Long[j])
    out.file.tempA<- paste(gsub(".nc", "", in.file), "_regrid.nc", sep = "")
    out.file.tempB<- paste(gsub(".nc", "", in.file), "_regrid_masked.nc", sep = "")
    out.file<- paste(gsub(".nc", "", in.file), "_regrid_masked_focal.nc", sep = "")
    run.command<- paste(paste("cd ", cd.to.path.use, sep = ""), 
                        paste(paste(" ; cdo remapbil,", regrid.use, sep = ""), paste(" -selname,tos,", in.file, out.file.tempA, sep = " "), sep = ""),
                        paste(paste("; cdo div "), paste(out.file.tempA, "landlakesmaskStdGridUse.nc", out.file.tempB, sep = " "), sep = ""),
                        #paste(paste(" ; cdo mul ", out.file.tempA, sep = ""), paste(" landlakesmaskStdGrid.nc", out.file.tempB, sep = " "), sep = ""), 
                        paste(paste(" ; cdo sellonlatbox,", lonlatbox, sep = ""), paste(out.file.tempB, out.file, sep = " ")), sep = "")
    system(run.command)
    
    # Add the info on the new files, needed later, to the clim.mods.df data frame
    clim.mods.df$File.Regrid[j]<- out.file.tempA
    clim.mods.df$File.Regrid.Masked[j]<- out.file.tempB
    clim.mods.df$File.Regrid.Masked.Focal[j]<- out.file
    
    # All done, next iteration
    print(paste(in.file, " is done!", sep = ""))
  }
  
  # Save the dataframe
  write_csv(clim.mods.df, path = paste(cd.to.path[[i]], "ClimModelData.csv", sep = ""))
}

## Raster stacks -- we are going to want to have a raster stack for each of the RCP45 and RCP85 models, which combines not only their information, but also the information from the matching historical runs. Let's do the raster conversion for the historical files first...
clim.mods.hist<- read_csv(paste("~/Box/Data/CMIP5_SST/Historical/", "ClimModelData.csv", sep = ""))
clim.mods.hist.groups<- unique(clim.mods.hist$Group)

for(i in seq_along(clim.mods.hist.groups)){
  
  clim.mods.sub<- clim.mods.hist %>%
    filter(., Group == unique(clim.mods.hist$Group)[i])
  
  # Empty raster stack to output processed data
  rast.stack.out<- stack()
  
  for(j in seq_along(clim.mods.sub$File.Regrid.Masked.Focal)){
    if(FALSE){
      rast.temp<- raster::stack(paste(cd.to.path.use, 
                                      "tos_Omon_IPSL-CM5B-LR_rcp85_r1i1p1_200601-210012_regrid_masked_focal.nc", sep = ""))
    }
    
    rast.temp<- raster::stack(paste("~/Box/Data/CMIP5_SST/Historical/", 
                                    clim.mods.sub$File.Regrid.Masked.Focal[j], sep = ""))
    names.use<- names(rast.temp)
    
    # Adjust values
    # Sometimes, land is 0, also sometimes have weird low values, makes these NA below Kelvin freezing point of salt water?
    rast.temp[rast.temp < 273.15]<- NA
    
    # Kelvin to Celsius
    rast.temp<- raster::calc(rast.temp, fun = function(x) x - 273.15)
    names(rast.temp)<- names.use
    
    # Save it
    rast.stack.out<- raster::stack(rast.stack.out, rast.temp)
  }
  
  # Write it out fully processed climate output for each of the modeling groups
  writeRaster(rast.stack.out, filename = paste("~/Box/Data/CMIP5_SST/Historical/", unique(clim.mods.hist$Group)[i], ".grd", sep = ""), overwrite = TRUE)
  
  print(paste(unique(clim.mods.hist$Group)[i], " is done!", sep = ""))
  
  
}

## Alright, now the CMP45 and CMP85 projections
clim.mods.rcp.paths<- list("RCP85" = "~/Box/Data/CMIP5_SST/RCP85/", "RCP45" = "~/Box/Data/CMIP5_SST/RCP45/")
clim.hist.rasts<- list.files(path = "~/Box/Data/CMIP5_SST/Historical", pattern = "*.*grd", full.names = TRUE)

for(i in seq_along(clim.mods.rcp.paths)){
  
  clim.mod.df.use<- read_csv(paste(clim.mods.rcp.paths[[i]], "ClimModelData.csv", sep = ""))
  
  #or(j in seq_along(unique(clim.mod.df.use$Group))){
  for(j in 18:length(unique(clim.mod.df.use$Group))){
    clim.mods.sub<- clim.mod.df.use %>%
      filter(., Group == unique(clim.mod.df.use$Group)[j])
    
    # Read in correct historical data
    clim.hist.rast<- raster::stack(clim.hist.rasts[which(grepl(clim.mods.sub$Group, clim.hist.rasts))])
    names.hist<- names(clim.hist.rast)
    
    # Empty raster stack to output projections data, which starts with historical data
    rast.stack.out<- stack(clim.hist.rast)
    names(rast.stack.out)<- names.hist
    
    for(k in seq_along(clim.mods.sub$File.Regrid.Masked.Focal)){
      if(FALSE){
        rast.temp<- raster::stack(paste(cd.to.path.use, 
                                        "tos_Omon_IPSL-CM5B-LR_rcp85_r1i1p1_200601-210012_regrid_masked_focal.nc", sep = ""))
      }
      
      rast.temp<- raster::stack(paste(clim.mods.rcp.paths[[i]], 
                                      clim.mods.sub$File.Regrid.Masked.Focal[k], sep = ""))
      names.use<- names(rast.temp)
      
      # Adjust values
      # Sometimes, land is 0, also sometimes have weird low values, makes these NA below Kelvin freezing point of salt water?
      rast.temp[rast.temp < 273.15]<- NA
      
      # Kelvin to Celsius
      rast.temp<- raster::calc(rast.temp, fun = function(x) x - 273.15)
      names(rast.temp)<- names.use
      
      # Save it
      rast.stack.out<- raster::stack(rast.stack.out, rast.temp)
    }
    
    # Write it out fully processed climate output for each of the modeling groups
    writeRaster(rast.stack.out, filename = paste(clim.mods.rcp.paths[[i]], unique(clim.mod.df.use$Group)[j], "Full.grd", sep = ""), overwrite = TRUE)
    
    # Make a quick plot to check visually
    # Need a time dimension...
    rast.time<- gsub("[.]", "-", gsub("X", "", names(rast.stack.out)))
    rast.temp.ts<- setZ(rast.stack.out, rast.time)
    rast.temp.agg<- zApply(rast.temp.ts, by = month, fun = mean, na.rm = TRUE)
    rast.temp.agg.df<- as.data.frame(rast.temp.agg, xy = TRUE)
    colnames(rast.temp.agg.df)[3:ncol(rast.temp.agg.df)]<- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
    rast.df.plot<- rast.temp.agg.df %>%
      gather(., "Month", "Temp", -x, -y) %>%
      mutate(., "Month.Factor" = factor(Month, levels = c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")))
    
    rast.check.plot<- ggplot() +
      geom_raster(data = rast.df.plot, aes(x = x, y = y, fill = Temp), na.rm = TRUE) +
      scale_fill_viridis(name = "SST", na.value = NA) +
      geom_sf(data = land.shp, fill = NA) +
      coord_sf(xlim = ext.x, ylim = ext.y, datum = NA) +
      facet_wrap(~Month.Factor, nrow = 4) +
      theme_bw()
    ggsave(paste(clim.mods.rcp.paths[[i]], unique(clim.mod.df.use$Group)[j], "PlotCheck.png", sep = ""), rast.check.plot)
    
    print(paste(unique(clim.mod.df.use$Group)[j], " is done!", sep = ""))
  }
}

## Git check....
 
