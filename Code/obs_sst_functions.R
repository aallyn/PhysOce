############## 
### Set of functions to facilitate extracting and visulizing observed temperature data in R
#############

# Helper functions --------------------------------------------------------
make360 <- function(lon) {
  ## Details
  # This is a simple function to translate negative longitudes (measured on -180:180 scale) into 0-360, which is coordinate system used by some environmental datasets.
  
  # Args:
  # lon = Longitude in -180:180 degrees
  
  # Returns: 0 to 360 longitude
  
  ## Start function
  
  ind <- which(lon < 0)
  lon[ind] <- lon[ind] + 360
  return(lon)
  
  ## End function
}

fix_raster<- function(x, lons.use, lats.use, x.min.use, x.max.use, y.min.use, y.max.use) {
  ## Details
  # This function helps convert an array (x, y, z) into raster layers with the correct dimensions and orientation.
  
  # Args:
  # x = A matrix. As currently implemented, x is an individiaul time slice from an x, y, z array.
  # lons.use = Longitudes, extracted using ncvar_get
  # lats.use = Latitudes, extracted using ncvar_get
  # x.min.use = Bounding box x minimum
  # x.max.use = Bounding box x maximum
  # y.min.use = Bounding box x minimum
  # y.max.use = Bounding box x maximum
  
  # Returns: Correctly oriented raster layer
  
  ## Start function
  r.temp<- t(x)[ncol(x):1,]
  rast.out<- raster(r.temp, xmn = min(lons.use), xmx = max(lons.use), ymn = min(lats.use), ymx = max(lats.use))
  return(rast.out)
  
  ## End function
}

fix_raster_clim<- function(x, lons.use, lats.use) {
  ## Details
  # This function helps convert an array (x, y, z) into raster layers with the correct dimensions and orientation.
  
  # Args:
  # x = A matrix. As currently implemented, x is an individiaul time slice from an x, y, z array.
  # lons.use = Longitudes, extracted using ncvar_get
  # lats.use = Latitudes, extracted using ncvar_get
  # x.min.use = Bounding box x minimum
  # x.max.use = Bounding box x maximum
  # y.min.use = Bounding box x minimum
  # y.max.use = Bounding box x maximum
  
  # Returns: Correctly oriented raster layer
  
  ## Start function
  r.temp<- t(x)[ncol(x):1,]
  rast.out<- raster(r.temp, xmn = min(lons.use), xmx = max(lons.use), ymn = min(lats.use), ymx = max(lats.use))
  return(rast.out)
  
  ## End function
}

oisst_list_to_stack<- function(x, ...) {
  ## Details
  # This function helps loop over daily OISST files
  
  # Args:
  # x = 
  
  # Returns: Stack of daily OISST files
  
  ## Start function
  # Install libraries
  library_check(c("ncdf4", "raster"))
  
  # Set arguments for debugging -- this will NOT run when you call the function. Though, you can run each line inside the {} and then you will have everything you need to walk through the rest of the function.
  if(FALSE){
    x = oisst.files
  }
  
  # OISST data stem
  data.stem<- "https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/noaa.oisst.v2.highres/"
  
  # Connecting and extracting lat/lon/time variables from netcdf file
  my.nc<- nc_open(paste(data.stem, x, sep = ""))
  lats<- ncvar_get(my.nc, var = "lat")
  lons<- ncvar_get(my.nc, var = "lon")
  times<- ncvar_get(my.nc, var = "time")
  
  # Make times a little bit easier to handle
  dates.full <- as.Date(times, origin='1800-01-01', tz= "GMT")
  
  # Find indices and windows corresponding to spatial box of interest, which are then used in the "start" and "count" arguments to the ncvar_get call for the sst variable
  b.box<- c(make360(box[1]), make360(box[2]), box[3], box[4])
  x.window<- which(lons > b.box[1] & lons < b.box[2])
  x.min<- min(x.window)
  x.max<- max(x.window)
  x.count<- ifelse(x.max-x.min > 0, x.max-x.min, 1)
  
  y.window<- which(lats > b.box[3] & lats < b.box[4])
  y.min<- min(y.window)
  y.max<- max(y.window)
  y.count<- ifelse(y.max - y.min > 0, y.max - y.min, 1) 
  
  if(!is.null(dates)) {
    times.window<- which(dates.full > as.Date(dates[1]) & dates.full < as.Date(dates[2]))
  } else {
    times.window<- times
  }
  
  time.min<- which.min(times.window)
  time.max<- which.max(times.window)
  time.count<- ifelse(time.max - time.min > 0, time.max - time.min, 1)
  
  # Now we have the lon,lat,time indices and windows, but need to match up their order with how they are handled in the ncvar_get call 
  dim.order <- sapply(my.nc$var$sst$dim, function(x) x$name)
  
  # Set start and counts
  start.use<- c("lon" = x.min, "lat" = y.min, "time" = time.min)
  count.use<- c("lon" = x.count, "lat" = y.count, "time" = time.count)
  
  # Run ncvar_get, adjusting order of start and count as needed
  temp<- ncvar_get(my.nc, varid = "sst", start = start.use[dim.order], count = count.use[dim.order])
  
  # Moving from the array format of temp to a raster stack
  temp.list<- lapply(seq(dim(temp)[3]), function(x) fix_raster(temp[,,x], lons.use = lons, lats.use = lats, x.min.use = x.min, x.max.use = x.max, y.min.use = y.min, y.max.use = y.max))
  rast.temp<- raster::rotate(stack(temp.list))
  return(rast.temp)
  
  ## End function
}

# Main satellite SST data extraction function ---------------------------------
obs_sst_extract<- function(data.set = "OISST", dates = NULL, box = c(-77, -60, 35, 46), out.dir = "~/Dropbox/Andrew/Work/GMRI/EnvData/", out.file) {
  
  ## Details
  # This function accesses webhosted satellite data and then downloads a subset of the data based on dates and long/lat bounding box. After downloading the data, the function processes it and saves it as one raster stack file. 
  
  # Args:
  # data.set = Env data to extract (options = ERSST, OISST, ...)
  # dates = If !NULL, subset full time series to specific dates. Dates should be specified as dates = c("YYYY-MM-DD", "YYYY-MM-DD")
  # box = If !NULL, crop rasters to sepcific box faster downloading and processing. Box should be specified as box = c(xmin, xmax, ymin, ymax).
  # out.dir = Directory to store resulting raster stack. Note, this will overwrite any rasters with the existing name. 

  # Returns: Stack of daily OISST files and also saves the raster stack as a .grd file in the output directory specified by out.dir
  
  ## Start function
  # Install libraries
  library_check(c("raster", "ncdf4"))

  # Set arguments for debugging -- this will NOT run when you call the function. Though, you can run each line inside the {} and then you will have everything you need to walk through the rest of the function.
  if(FALSE){
    data.set<- "OISST"
    dates<- c("2001-01-01", "2008-01-01")
    box<- c(-77, -60, 34, 46)
    out.dir<- "~/Dropbox/Andrew/Work/GMRI/EnvData/"
    out.file<- "NWAtl"
  }
  
  # Spatial projection infomation
  proj.wgs84<- CRS("+init=epsg:4326") #WGS84 projection for all rasters
  
  # Access data from the web depending on "data.set"
  # MURSST-- Download of full timeseries takes ~ 20 minutes.
  if(data.set == "MURSST"){
    
    # Set path to THREDDS server
    data.path<- "http://thredds.jpl.nasa.gov/thredds/dodsC/OceanTemperature/MUR-JPL-L4-GLOB-v4.1.nc"
    
    # Some steps to make things cooperate with nc functions. MURSST is already in -180 to 180, so we are good there..
    b.box<- c(box[1], box[2], box[3], box[4])
    
    # Connecting and extracting lat/lon/time variables from netcdf file
    my.nc<- nc_open(data.path)
    lats<- ncvar_get(my.nc, var = "lat")
    lons<- ncvar_get(my.nc, var = "lon")
    times<- ncvar_get(my.nc, var = "time")
    
    # Make times a little bit easier to handle
    dates.full<- as.Date(as.POSIXct(times, origin="1981-01-01", tz = "GMT"), "GMT", "%Y-%m-%d")
    
    x.window<- which(lons > b.box[1] & lons < b.box[2])
    x.min<- min(x.window)
    x.max<- max(x.window)
    x.count<- ifelse(x.max - x.min > 0, x.max-x.min, 1)
    
    y.window<- which(lats > b.box[3] & lats < b.box[4])
    y.min<- min(y.window)
    y.max<- max(y.window)
    y.count<- ifelse(y.max - y.min > 0, y.max - y.min, 1) 
    
    if(!is.null(dates)) {
      times.window<- which(dates.full > as.Date(dates[1]) & dates.full < as.Date(dates[2]))
    } else {
      times.window<- dates.full
    }
    
    time.min<- which.min(times.window)
    time.max<- which.max(times.window)
    
    # MURSST files are taking a while to get, so going to try to split up the processing...
    splits<- 569
    breaks<- seq(1, length(times), splits)[-1]
    
    for(i in seq_along(breaks)){
      
      # Find indices and windows corresponding to spatial box of interest, which are then used in the "start" and "count" arguments to the ncvar_get call for the sst variable
      if(i == 1){
        time.count<- seq(from = 1, to = breaks[i], by = 1)
        dates.use<- dates.full[1:breaks[i]]
      } else {
        time.count<- seq(from = breaks[i-1]+1, to = breaks[i], by = 1)
        dates.use<- dates.full[(breaks[i-1]+1):breaks[i]]
      }
      
      # Now we have the lon,lat,time indices and windows, but need to match up their order with how they are handled in the ncvar_get call 
      dim.order <- sapply(my.nc$var$analysed_sst$dim, function(x) x$name)
      
      # Set start and counts
      start.use<- c("lon" = x.min, "lat" = y.min, "time" = time.min)
      count.use<- c("lon" = x.count, "lat" = y.count, "time" = length(time.count))
      
      # Run ncvar_get, adjusting order of start and count as needed
      temp<- ncvar_get(my.nc, varid = "analysed_sst", start = start.use[dim.order], count = count.use[dim.order])
      
      # Moving from the array format of temp to a raster stack
      temp.list<- lapply(seq(dim(temp)[3]), function(x) fix_raster(temp[,,x], lons.use = lons, lats.use = lats, x.min.use = x.min, x.max.use = x.max, y.min.use = y.min, y.max.use = y.max))
      rast.temp<- raster::stack(temp.list)
      rast.temp<- raster::calc(rast.temp, fun = function(x) x - 273.15)
      names(rast.temp)<- dates.use
      
      # Update min time
      time.min<- breaks[i] + 1
      
      # Write out chunk
      proj4string(rast.temp)<- proj.wgs84
      writeRaster(rast.temp, filename = paste(out.dir, data.set, ".chunk", i, ".grd", sep = ""), overwrite = TRUE)
      
      print(paste("MURSST chunk ", i, " out of ", length(breaks), " is done!"))
    }
  }
  
  # ERSST
  if(data.set == "ERSST") {
    # Set data path
    data.path<- "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/noaa.ersst/sst.mnmean.v4.nc"
    
    # Some steps to make things cooperate with nc functions, including getting variables from -180:180 to 0-360
    b.box<- c(make360(box[1]), make360(box[2]), box[3], box[4])
    
    # Connecting and extracting lat/lon/time variables from netcdf file
    my.nc<- nc_open(data.path)
    lats<- ncvar_get(my.nc, var = "lat")
    lons<- ncvar_get(my.nc, var = "lon")
    times<- ncvar_get(my.nc, var = "time")
    
    # Make times a little bit easier to handle
    dates.full <- as.Date(times, origin='1800-01-01', tz= "GMT")
    
    # Find indices and windows corresponding to spatial box of interest, which are then used in the "start" and "count" arguments to the ncvar_get call for the sst variable
    x.window<- which(lons > b.box[1] & lons < b.box[2])
    x.min<- min(x.window)
    x.max<- max(x.window)
    x.count<- ifelse(x.max - x.min > 0, x.max-x.min, 1)
    
    y.window<- which(lats > b.box[3] & lats < b.box[4])
    y.min<- min(y.window)
    y.max<- max(y.window)
    y.count<- ifelse(y.max - y.min > 0, y.max - y.min, 1) 
    
    if(!is.null(dates)) {
      times.window<- which(dates.full > as.Date(dates[1]) & dates.full < as.Date(dates[2]))
    } else {
      times.window<- dates.full
    }
    
    time.min<- which.min(times.window)
    time.max<- which.max(times.window)
    time.count<- ifelse(time.max - time.min > 0, time.max - time.min, 1)
    
    # Now we have the lon,lat,time indices and windows, but need to match up their order with how they are handled in the ncvar_get call 
    dim.order <- sapply(my.nc$var$sst$dim, function(x) x$name)
    
    # Set start and counts
    start.use<- c("lon" = x.min, "lat" = y.min, "time" = time.min)
    count.use<- c("lon" = x.count, "lat" = y.count, "time" = time.count)
    
    # Run ncvar_get, adjusting order of start and count as needed
    temp<- ncvar_get(my.nc, varid = "sst", start = start.use[dim.order], count = count.use[dim.order])
    
    # Moving from the array format of temp to a raster stack
    # Adjustment for when x and y count are 1...
    if(length(dim(temp)) == 1){
      temp.list<- lapply(seq(dim(temp)[1]), function(x) raster(as.matrix(temp[x]), xmn = lons[x.min], xmx = lons[x.max], ymn = lats[y.max], ymx = lats[y.min]))
    } else {
      temp.list<- lapply(seq(dim(temp)[3]), function(x) raster(temp[ , , x], xmn = lons[x.min], xmx = lons[x.max], ymn = lats[y.max], ymx = lats[y.min]))
    }
    
    stack.out<- raster::rotate(stack(temp.list))
    stack.out[stack.out == 327.16]<- NA
    names(stack.out)<- dates.full[-length(dates.full)]
    
    # Write out raster stack
    proj4string(stack.out)<- proj.wgs84
    writeRaster(stack.out, filename = paste(out.dir, data.set, ".grd", sep = ""), overwrite = TRUE)
    return(stack.out)
  }
  
  # OISST -- Download of full timeseries takes ~ 20 minutes.
  if(data.set == "OISST") {
    stem.path<- "http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/noaa.oisst.v2.highres/"
    
    b.box<- c(make360(box[1]), make360(box[2]), box[3], box[4])
    
    # Get file list
    if(!is.null(dates)){
      # Get time series subset
      date.min<- format(as.Date(dates[1]), "%Y")
      date.max<- format(as.Date(dates[2]), "%Y")
      files<- paste(stem.path, paste("sst.day.mean.", seq(from = date.min, to = date.max, by = 1), ".v2.nc", sep = ""), sep = "")
    } else {
      current.year<- format(Sys.Date(), "%Y")
      files<- paste(stem.path, paste("sst.day.mean.", seq(from = 1981, to = current.year, by = 1), ".v2.nc", sep = ""), sep = "")
    }
    
    stack.out<- stack()
    
    for(i in seq_along(files)){
      # Connecting and extracting lat/lon/time variables from netcdf file
      my.nc<- nc_open(files[i])
      lats<- ncvar_get(my.nc, var = "lat")
      lons<- ncvar_get(my.nc, var = "lon")
      times<- ncvar_get(my.nc, var = "time")
      
      # Make times a little bit easier to handle
      dates.full <- as.Date(times, origin='1800-01-01', tz= "GMT")
      
      # Find indices and windows corresponding to spatial box of interest, which are then used in the "start" and "count" arguments to the ncvar_get call for the sst variable
      x.window<- which(lons > b.box[1] & lons < b.box[2])
      x.min<- min(x.window)
      x.max<- max(x.window)
      x.count<- ifelse(x.max - x.min > 0, x.max-x.min, 1)
      
      y.window<- which(lats > b.box[3] & lats < b.box[4])
      y.min<- min(y.window)
      y.max<- max(y.window)
      y.count<- ifelse(y.max - y.min > 0, y.max - y.min, 1) 
      
      time.min<- which.min(dates.full)
      time.max<- which.max(dates.full)
      time.count<- ifelse(time.max - time.min > 0, (time.max - time.min)+1, 1)
      
      # Now we have the lon,lat,time indices and windows, but need to match up their order with how they are handled in the ncvar_get call 
      dim.order <- sapply(my.nc$var$sst$dim, function(x) x$name)
      
      # Set start and counts
      start.use<- c("lon" = x.min, "lat" = y.min, "time" = time.min)
      count.use<- c("lon" = x.count, "lat" = y.count, "time" = time.count)
      
      # Run ncvar_get, adjusting order of start and count as needed
      temp<- ncvar_get(my.nc, varid = "sst", start = start.use[dim.order], count = count.use[dim.order])
      
      # Moving from the array format of temp to a raster stack
      temp.list<- lapply(seq(dim(temp)[3]), function(x) fix_raster(temp[,,x], lons.use = lons, lats.use = lats, x.min.use = x.min, x.max.use = x.max, y.min.use = y.min, y.max.use = y.max))
      rast.temp<- raster::rotate(stack(temp.list))
      
      if(i == 1) {
        stack.out<- rast.temp
      } else {
        stack.out<- stack(stack.out, rast.temp)
      }
      print(paste(files[i], "is done", sep = " "))
    }
    
    # Names
    stack.layers<- nlayers(stack.out)
    
    if(!is.null(dates)){
      stack.names<- seq(from = as.Date(dates[1]), to = as.Date(dates[1]) + (nlayers(stack.out)-1), by = "day")
      names(stack.out)<- stack.names
    } else {
      start.name<- as.Date("1981-09-01")
      end.name<- (start.name + stack.layers)-1
      stack.names<- seq(from = start.name, to = end.name, by = "day")
      names(stack.out)<- stack.names
    }
  }
  
  # Write out raster stack
  date.out<- format(Sys.Date(), "%Y%M")
  proj4string(stack.out)<- proj.wgs84
  writeRaster(stack.out, filename = paste(out.dir, data.set, date.out, ".grd", sep = ""), overwrite = TRUE)
  return(stack.out)
}

# Env data time series plot function --------------------------------------
obs_sst_timeseries<- function(oisst.stack, baseline = c("1982-01-01", "2011-01-01"), alt.trend = c("2004-01-01", "2018-12-31"), regions = c("NELME", "GoM", "SNE-MAB"), out.dir = "/Volumes/Shared/Research/Mills Lab/SST/") {
  ## Details
  # This function plots time series (for potentially different regions) and fits trend lines to evaluate warming rates
  
  # Args:
  # oisst.stack = Path to the raster OISST stack created by env_data_extract function
  # baseline = Start and end dates for the baseline period. Dates should be specified as dates = c("YYYY-MM-DD", "YYYY-MM-DD"). Defaults to 1982-2011.
  # Alt trend = Start and end dates for an alternative trendline to add to the plot. Dates should be specified as alt.trend = c("YYYY-MM-DD", "YYYY-MM-DD").
  # regions = Character vector determining what regions should be analyzed. Can be any combination of "NELME", "GoM", "SNE-MAB"
  # out.dir = Path to output both SST flat csv files and time series plots 
  
  # Returns: NA. Just outputs files and figures
  
  ## Start function
  # Install libraries
  library_check(c("tidyverse", "sf", "zoo", "forecast", "raster"))
  
  # Set arguments for debugging -- this will NOT run when you call the function. Though, you can run each line inside the {} and then you will have everything you need to walk through the rest of the function.
  if(FALSE){
    oisst.stack = "~/Dropbox/Andrew/Work/GMRI/Projects/AllData/OISST.grd"
    baseline = c("1982-01-01", "2011-01-01")
    alt.trend = c("2004-01-01", "2018-12-31")
    regions = c("NELME", "GoM", "SNE-MAB")
    out.dir<- "~/Dropbox/Andrew/Work/GMRI/Projects/AllData/"
  }
  
  # Read in OISST data
  oisst.dat<- raster::stack(oisst.stack)
  oisst.dates<- gsub("[.]", "-", gsub("X", "", names(oisst.dat)))
  oisst.rts<- setZ(oisst.dat, oisst.dates)
  
  # Extracting temperatures for different regions, saving files and plotting time series
  for(i in seq_along(regions)){
    # Focus on only region of interest
    mask.use<- switch(regions[i],
                      "NELME" = st_read("/Volumes/Shared/Research/Mills Lab/SST/Shapefiles/NELME_sf.shp"),
                      "GoM" = st_read("/Volumes/Shared/Research/Mills Lab/SST/Shapefiles/GoM_sf.shp"),
                      "SNE-MAB" = st_read("/Volumes/Shared/Research/Mills Lab/SST/Shapefiles/SNEandMAB_sf.shp"))
    oisst.m<- mask(oisst.rts, mask.use)
    
    # Need to get climatology from the OISST data 
    # Baseline
    oisst.m.base<- oisst.m[[which(getZ(oisst.m) >= baseline[1] & getZ(oisst.m) <= baseline[2])]]
    oisst.m.base<- setZ(oisst.m.base, seq.Date(from = as.Date(baseline[1]), to = as.Date(baseline[2]), by = "day"))
    
    dates.unique<- unique(format(as.Date(getZ(oisst.m)), "%m-%d"))
    daily.means<- stack(lapply(seq(length(dates.unique)), function(x) calc(oisst.m.base[[which(format(getZ(oisst.m.base), "%m-%d") == dates.unique[x])]], fun = mean, na.rm = TRUE)))
    names(daily.means)<- dates.unique
    daily.sd<- stack(lapply(seq(length(dates.unique)), function(x) calc(oisst.m.base[[which(format(getZ(oisst.m.base), "%m-%d") == dates.unique[x])]], fun = sd, na.rm = TRUE)))
    names(daily.sd)<- dates.unique
    
    # Daily temps for KM
    daily.temps<- do.call("cbind", lapply(seq(1:nlayers(oisst.m)), function(x) as.data.frame(oisst.m[[x]], xy = TRUE)))
    daily.temps2<- daily.temps %>%
      subset(., select=which(!duplicated(names(.)))) %>%
      gather(., Year, SST, -x, -y) 
    
    daily.temps2$Year<- gsub("[.]", "-", gsub("X", "", daily.temps2$Year))
    daily.dat<- daily.temps2
    names(daily.dat)[3]<- "Date"
    daily.dat$Date<- as.Date(daily.dat$Date)
    daily.dat<- daily.dat %>%
      drop_na(SST)
    
    # Daily data
    write_csv(daily.dat, path = paste(out.dir, regions[i], ".SSTDailyDegF.csv", sep = ""))
    
    # Baseline data
    clim.temps<- do.call("cbind", lapply(seq(1:nlayers(daily.means)), function(x) as.data.frame(daily.means[[x]], xy = TRUE)))
    clim.temps2<- clim.temps %>%
      subset(., select=which(!duplicated(names(.)))) %>%
      gather(., Day, SST, -x, -y) 
    
    clim.temps2$Day<- paste("1944", gsub("[.]", "-", gsub("X", "", clim.temps2$Day)), sep = "-")
    clim.dat<- clim.temps2
    names(clim.dat)[3]<- "Date"
    clim.dat$Date<- as.Date(clim.dat$Date)
    clim.dat<- clim.dat %>%
      drop_na(SST) %>%
      group_by(Date) %>%
      summarize(., Mean.SST = mean(SST))
    
    # Clim data
    write.csv(clim.dat, file = paste(out.dir, regions[i], ".SSTDailyClimatologyDegF.csv", sep = ""))
    
    # Alright, now substract each daily OISST from the daily climatology to get the anomaly
    anom.type<- "Non.standardized"
    daily.anoms<- switch(anom.type, 
                         Standardized = stack(lapply(seq(1:nlayers(oisst.m)), function(x) (oisst.m[[x]] - daily.means[[match(format(as.Date(getZ(oisst.m)[x]), "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.means))))]])/daily.sd[[match(format(as.Date(getZ(oisst.m)[x]), "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.sd))))]])),
                         Non.standardized = stack(lapply(seq(1:nlayers(oisst.m)), function(x) (oisst.m[[x]] - daily.means[[match(format(as.Date(getZ(oisst.m)[x]), "%m-%d"), gsub("[.]", "-", gsub("X", "", names(daily.means))))]]))))
    names(daily.anoms)<- getZ(oisst.m)
    
    # Save it
    writeRaster(daily.anoms, filename = paste(out.dir, regions[i], ".grd", sep = ""), overwrite = TRUE)
    
    # Getting the anomalies
    ts.wide.daily<- do.call("cbind", lapply(seq(1:nlayers(daily.anoms)), function(x) as.data.frame(daily.anoms[[x]], xy = TRUE)))
    ts.df.daily<- ts.wide.daily %>%
      subset(., select=which(!duplicated(names(.)))) %>%
      gather(., Year, SST, -x, -y) 
    names(ts.df.daily)[3]<- "Date"
    ts.df.daily$Date<- gsub("X", "", gsub("[.]", "-", ts.df.daily$Date))
    ts.df.daily<- ts.df.daily %>%
      separate(., Date, c("YYYY", "MM", "DD")) %>%
      group_by(., YYYY, MM, DD) %>%
      summarize_at(., "SST", mean, na.rm = T)
    
    # Save em
    write.csv(ts.df.daily, paste(out.dir, regions[i], "DailyAnomalies.csv", sep = ""))
    
    # Keep going for individual region...
    ts.df.dailymu<- ts.df.daily %>%
      mutate(., "Plot.Date" = as.Date(paste(YYYY, MM, DD, sep = "-"))) %>%
      data.frame
    
    # Smooth daily values
    # Creat zoo series
    zoo.dates<- as.Date(ts.df.dailymu$Plot.Date)
    zoo.vals<- zoo(ts.df.dailymu$SST, zoo.dates)
    zoo.ts<- as.ts(zoo.vals)
    ts.df.dailymu$smoothed<- as.numeric(ma(zoo.vals, order = 15, centre = FALSE))
    
    ts.df.monthlymu<- ts.df.dailymu %>%
      group_by(YYYY, MM) %>%
      dplyr::summarize(., Mean.SST = mean(SST, na.rm =T)) %>%
      mutate(., Plot.Date = as.Date(paste(YYYY, MM, "15", sep = "-"))) %>%
      data.frame
    
    ts.df.yearlymu<- ts.df.dailymu %>%
      group_by(YYYY) %>%
      dplyr::summarize(., Mean.SST = mean(SST, na.rm = T)) %>%
      mutate(., Plot.Date = as.Date(paste(YYYY, "06", "15", sep = "-"))) %>%
      filter(., YYYY != format(Sys.time(), "%Y")) %>%
      data.frame()
    
    ts.df.yearlymu$Year.Model<- as.integer(ts.df.yearlymu$YYYY) - 1982
    
    # Plots
    d15<- TRUE
    if(d15 == TRUE){
      # Get current year, then plot all data before the current year
      curr.year<- format(as.Date(Sys.time()), "%Y")
      ts.end<- paste(as.numeric(curr.year)-1, "-12-31", sep = "")
      
      # Full trend line
      sst.anom.lm.full<- lm(Mean.SST ~ Year.Model, data = ts.df.yearlymu)
      summary(sst.anom.lm.full)
      adj.r2.full<- paste("Adj R2 = ", round(summary(sst.anom.lm.full)$adj.r.squared, 3), sep = "")
      
      # Fitted model formula
      my.formula.full<- paste("Mean.SST = ", signif(round(sst.anom.lm.full$coef[[1]], 3), 5), " + ", signif(round(sst.anom.lm.full$coef[[2]], 3), 5), "*Year", sep = "")
      base.ts<- ggplot(data = subset(ts.df.dailymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date(ts.end, format = "%Y-%m-%d")), aes(x = Plot.Date, y = smoothed)) + 
        geom_line(col = "#d9d9d9", lwd = 0.15) +
        geom_point(data = subset(ts.df.yearlymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date(ts.end, format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), col = "black") +
        geom_line(data = subset(ts.df.yearlymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date(ts.end, format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), col = "black", lwd = 0.25) +
        geom_smooth(data = subset(ts.df.yearlymu, Plot.Date >= as.Date("1982-01-01", format = "%Y-%m-%d") & Plot.Date <= as.Date(ts.end, format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), method = "lm", formula = y ~ x, col = "#969696", size = 0.75, se = FALSE) +
        ylim(c(-3, 4)) +
        ylab("SST Anomaly (1982-2011 baseline)") + 
        xlab("Year") +
        theme_bw() +
        theme(axis.text=element_text(size=14),
              axis.title=element_text(size=16, face="bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        geom_text(aes(x = as.Date("1987-06-15"), y = -2.5, label = my.formula.full), col = "#969696") +
        geom_text(aes(x = as.Date("1997-06-15"), y = -2.5, label = adj.r2.full), col = "#969696") 
      ggsave(paste(out.dir, "sstanomaly_BaselineOnly.jpg", sep = ""), width = 8, height = 6, units = "in")
      dev.off()
      
      # Add a more recent trend line
      sst.anom.lm.plot2<- lm(Mean.SST ~ Year.Model, data = subset(ts.df.yearlymu, YYYY >= format(as.Date(alt.trend[1]), "%Y") & YYYY <= format(as.Date(alt.trend[2]), "%Y")))
      summary(sst.anom.lm.plot2)
      adj.r2.plot2<- paste("Adj R2 = ", round(summary(sst.anom.lm.plot2)$adj.r.squared, 3), sep = "")
      
      # Fitted model formula
      my.formula.plot2<- paste("Mean.SST = ", signif(round(sst.anom.lm.plot2$coef[[1]], 3), 5), " + ", signif(round(sst.anom.lm.plot2$coef[[2]], 3), 5), "*Year", sep = "")
      base.ts.alttrend<- base.ts +
        geom_smooth(data = subset(ts.df.yearlymu, Plot.Date >= as.Date(alt.trend[1], format = "%Y-%m-%d") & Plot.Date <= as.Date(alt.trend[2], format = "%Y-%m-%d")), aes(x = Plot.Date, y = Mean.SST), method = "lm", formula = y ~ x, col = "red", size = 0.75, se = FALSE) +
        geom_text(aes(x = as.Date("1987-06-15"), y = -2.8, label = my.formula.plot2), col = "red") +
        geom_text(aes(x = as.Date("1997-06-15"), y = -2.8, label = adj.r2.plot2), col = "red") 
      ggsave(paste(out.dir, "sstanomaly_Baseand", format(as.Date(alt.trend[1]), "%Y"), "to", format(as.Date(alt.trend[2]), "%Y"), ".jpg", sep = ""), width = 8, height = 6, units = "in")
      dev.off()
    }
    print(paste(regions[i], " is done!", sep = ""))
  }
  # End function
}

# FVCOM data extraction function ------------------------------------------
fvcom_extract<- function(threddsURL = "http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean", lonW = -82, lonE = -61.75, latS = 24, latN = 46.25, spacing, date1 = NULL, date2 = NULL, variable, out.path) {
  # This function extracts FVCOM data from the threddsURL. After extracting the data (and subsetting if dates are specified), it saves a the raw data as a large matrix and also as a regular gridded raster stack. To retrieve the full time series takes about ~5 minutes.
  
  # Args:
  # threddsURL = Link to FVCOM THREDDS URL. Default is "http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean" A quick note here, this path is actually found at the page "http://www.smast.umassd.edu:8080/thredds/hindcasts.html?dataset=fvcom/hindcasts/30yr_gom3/mean" and under the "Access" options. I think if you wanted to get something else, say the FVCOM forecasts, you would pull the THREDDS URL from this page: http://www.smast.umassd.edu:8080/thredds/forecasts.html?dataset=gom3_nocache. You would then probably need to do some function modifications/checking to make sure everything still works as you expect. 
  # lonW = Left bounding box decimal degree coordinate for raster stack
  # lonE = Right bounding boc decimal degree coordinate for raster stack
  # latS = Bottom bounding box decimal degree coordinate for raster stack
  # latN = Upper bounding box decimal degree coordinate for raster stack
  # spacing = Grid spacing for output raster stack. The FVCOM data is not on a regular grid, which is a requirement of raster stacks. To get around this, the function interpolates the FVCOM data onto a raster layer, defined by lon, lat and spacing. 
  # date1 = Either a date string ("YYYY-MM-DD") to set beginning of time series subset, or NULL to retrieve full time series. Default is NULL.
  # date2 = Either a date string ("YYYY-MM-DD") to set end of time series subset, or NULL to retrieve full time series. Default is NULL.
  # variable = FVCOM variable data to extract. Default is temp, could also be "salinity" -- not sure on others.
  # out.path = Path to save large matrix and raster stack FVCOM data
  
  # Returns: Raster stack of FVCOM variable data. 
  
  # Preliminaries
  library_check(c("ncdf4", "raster", "fields"))
  
  # Debugging
  if(FALSE) {
    threddsURL<- "http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3/mean"
    lonW<- -82
    lonE<- -61.75
    latS<- 24
    latN<- 46.25
    spacing<- 0.15
    date1<- NULL
    date2<- NULL
    variable<- "temp"
    out.path<- "~/Desktop/"
  }
  
  ## Start
  # Open file connection, extract lon/lat info
  nc<- nc_open(threddsURL)
  lon.nc<- ncvar_get(nc, varid = "lon")
  nlons<- dim(lon.nc)
  lat.nc<- ncvar_get(nc, varid = "lat")
  nlats<- dim(lat.nc)
  
  # Extract available dates from netCDF file
  ncdates<- nc$dim$time$vals
  ncdates<- as.Date(ncdates,origin = '1858-11-17') #available time points in nc
  
  # Subset full time series?
  if(is.null(date1) && is.null(date2)) {
    date1<- min(ncdates)
    date2<- max(ncdates)
  }
  
  if (class(date1) == 'Date'){
    # Get index of nearest time point
    date1indx = which.min(abs(date1 - ncdates)) 	
  } else if (class(date1) == 'character'){
    # Convert to a Date object first
    date1 = as.Date(date1)
    date1indx = which.min(abs(date1 - ncdates)) 
  }
  if (missing(date2)) {
    # If date2 isn't specified, reuse date1
    date2indx = which.min(abs(date1 - ncdates)) 
    cat('Only 1 date specified\n')
  } else {
    if (class(date2) == 'Date'){
      # If date2 exists, get index of nearest time point to date2
      date2indx = which.min(abs(date2 - ncdates)) 		
    } else if (class(date2) == 'character'){
      date2 = as.Date(date2)
      date2indx = which.min(abs(date2 - ncdates))
    }
  }
  
  # Get number of time steps to extract
  ndates <- (date2indx - date1indx) + 1 #get number of time steps to extract
  
  # Get variable data
  var.fvcom<- ncvar_get(nc, varid = variable, start =  c(1, 45, 1))
  
  # Little clean up to save all data as a large 
  fvcom.df.out<- data.frame("lon" = lon.nc, "lat" = lat.nc, var.fvcom)
  colnames(fvcom.df.out)[3:ncol(fvcom.df.out)]<- paste(variable, ncdates, sep = ".")
  write_csv(fvcom.df.out, paste(out.path, "FVCOM", variable, ".csv", sep = ""))
  
  # Now need to convert to a raster. In the var.fvcom matrix, each row is a point and the columns are the different dates. Need to extract and put into raster stack. One issue with this, though, is that the FVCOM data are not on a perfect grid -- so some interpolation is going to be needed. To do this, set up an empty raster layer
  rast.lon<- seq(from = lonE, to = lonW, by = -spacing)
  rast.lat<- seq(from = latS, to = latN, by = spacing)
  e<- extent(c(lonW, lonE, latS, latN))
  rast.temp<- raster(e, nrows = length(rast.lat), ncols = length(rast.lon))
  
  stack0<- stack()
  
  for (i in 1:ncol(var.fvcom)) {
    t1<- data.frame("long" = lon.nc, "lat" = lat.nc, "z" = var.fvcom[,i])
    rast<- rasterize(t1[,1:2], rast.temp, t1[,3], fun = mean)
    stack0<- stack(stack0, rast)
    print(paste(ncdates[i], " is done", sep = ""))
  }
  
  # Raster stack spatial projection and names
  proj4string(stack0)<- CRS("+init=epsg:4326") #WGS84
  names(stack0)<- as.character(ncdates)
  
  # Save it and return it
  writeRaster(stack0, paste(out.path, "FVCOM", variable, ".grd", sep = ""))
  return(stack0)
}

