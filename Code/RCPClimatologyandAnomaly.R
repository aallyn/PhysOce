# Preliminaries -----------------------------------------------------------
source("https://raw.githubusercontent.com/GMRI-SEL/LabFunctionsandCode/master/GenerateSharedPathsScript.R")
source(paste(lab.func.path, "GeneralHelpers.R", sep = ""))
source(here::here("Code", "clim_sst_functions.R"))
library_check(c("tidyverse", "sf", "lubridate", "here", "raster", "rts"))

# List of each individual raster stack -- for each one we are going to want to calculate the climatology, and then the anomalies
res.data.path<- "/Users/aallyn/Box/Data/"
all.rast.ts<- tibble("Raster.TS" = c(paste(res.data.path, "OISST/OISSTNWAtl.grd", sep = ""), unlist(list.files(paste(res.data.path, "CMIP5_SST/RCP85", sep = ""), pattern = "Full.grd", full.names = TRUE)), list.files(paste(res.data.path, "CMIP5_SST/RCP45", sep = ""), pattern = "Full.grd", full.names = TRUE)))

# Output path
all.rast.ts<- all.rast.ts %>%
  separate("Raster.TS", c("Path", "File"), sep = "/(?!.*/)", remove = FALSE)
all.rast.ts$res.ts<- ifelse(grepl("OISST", all.rast.ts$Raster.TS), "daily", "monthly")
  
# First, climatology with calc_climatology
all.rast.ts<- all.rast.ts %>%
  mutate(., "Raster.Clim" = pmap(list(raster.ts = Raster.TS, res.ts = res.ts, res.clim = "monthly", baseline = list(c("1982-01-01", "2011-01-01")), out.path = Path, out.file = File), calc_climatology))

# Okay, got the climatology, now need to calculate the anomalies. To do that, we probably want to split out the Raster.Clim list column and have two new columns...one for the mean and one for the SD.
rast.clim.temp<- all.rast.ts %>%
  dplyr::select(., Raster.TS, Raster.Clim) %>%
  unnest(cols = c(Raster.Clim)) 

rast.clim.mu<- rast.clim.temp[seq(from = 1, to = nrow(rast.clim.temp), by = 2),]
colnames(rast.clim.mu)[2]<- "Clim.Mu"
rast.clim.sd<- rast.clim.temp[seq(from = 2, to = nrow(rast.clim.temp), by = 2),]
colnames(rast.clim.sd)[2]<- "Clim.SD"
rast.clim.merge<- rast.clim.mu %>%
  left_join(., rast.clim.sd)

all.rast.ts<- all.rast.ts %>%
  left_join(., rast.clim.merge)

# Only going to do this for the climate models, not the OISST (row 1)
clim.rast.ts<- all.rast.ts[-1, ]

clim.rast.ts<- clim.rast.ts %>%
  mutate(., "Clim.Anom" = pmap(list(raster.ts = Raster.TS, raster.clim.mu = Clim.Mu, raster.clim.sd = Clim.SD, anom.type = "Standardized", out.path = Path, out.file = File), calc_anomaly))

# Now, going to apply the anomalies to the observed OISST climatology to get a biased corrected projected temperature.
clim.rast.ts<- clim.rast.ts %>%
  mutate(., "Proj.SST" = pmap(list(raster.anom = Clim.Anom, raster.clim = list(all.rast.ts$Clim.Mu[[1]]), out.path = Path, out.file = File), calc_proj_sst))

# Alright, now need to go from a raster stack Proj.SST to a data frame then calculate the mean and 95% CI
stack_to_df<- function(rast.stack){
  df.out<- as.data.frame(rast.stack, xy = TRUE)
  df.out<- df.out %>%
    gather(., "Year", "SST", -x, -y)
  return(df.out)
}

clim.rast.ts<- clim.rast.ts %>%
  mutate(., "Proj.SST.DF" = map(Proj.SST, stack_to_df))
saveRDS(clim.rast.ts, "~/Box/Data/CMIP5_SST/ProcessedSSTProjectionsWithRasters.rds")


clim.out<- clim.rast.ts %>%
  dplyr::select(., Raster.TS, File, Proj.SST.DF) %>%
  unnest()

clim.out$Year<- gsub("X", "", gsub("[.]", "-", clim.out$Year))
clim.out$Scenario<- ifelse(grepl("85", clim.out$Raster.TS), "RCP85", "RCP45")

clim.out<- clim.out %>%
  dplyr::filter(., Year >= "1982-01-01") %>%
  dplyr::select(., File, Scenario, Year, x, y, SST)
saveRDS(clim.out, "~/Box/Data/CMIP5_SST/ProcessedSSTProjections.rds")


# Compare our pulls with NCAR data ----------------------------------------
clim.out<- readRDS("~/Box/Data/CMIP5_SST/ProcessedSSTProjections.rds")

# Mean, rcp45 and rcp85
clim.summs<- clim.out %>%
  separate(Year, into = c("Year", "Month", "Day")) %>%
  group_by(Scenario, Year, Month, x, y) %>%
  summarize("Mean" = mean(SST, na.rm = TRUE),
            "Pct5th" = quantile(SST, probs = c(0.05), na.rm = TRUE, names = FALSE),
            "Pct95th" = quantile(SST, probs = c(0.95), na.rm = TRUE, names = FALSE))

clim.summs<- clim.summs %>%
  group_by(Scenario, Year, Month) %>%
  nest()

df_to_rast<- function(df, stat, mask = nelme) {
  if(FALSE){
    df<- clim.summs$data[[1]]
    stat<- "Mean.SST"
  }
  
  df.temp<- df %>%
    dplyr::select(x, y, stat)
  
  rast.temp<- rasterFromXYZ(df.temp)
  rast.out<- mask(rast.temp, nelme)
  return(rast.out)
}

clim.summs<- clim.summs %>%
  mutate(., "RasterStack.Mean" = map2(data, "Mean", df_to_rast),
         "RasterStack.Pct05" = map2(data, "Pct5th", df_to_rast),
         "RasterStack.Pct95" = map2(data, "Pct95th", df_to_rast))

rcp85<- clim.summs %>%
  dplyr::filter(., Scenario == "RCP85")
rcp45<- clim.summs %>%
  dplyr::filter(., Scenario == "RCP45")

rcp85.mu<- raster::stack(rcp85$RasterStack.Mean)
names(rcp85.mu)<- paste(rcp85$Year, rcp85$Month, sep = ".")
rcp85.5th<- raster::stack(rcp85$RasterStack.Pct05)
names(rcp85.5th)<- paste(rcp85$Year, rcp85$Month, sep = ".")
rcp85.95th<- raster::stack(rcp85$RasterStack.Pct95)
names(rcp85.95th)<- paste(rcp85$Year, rcp85$Month, sep = ".")

rcp45.mu<- raster::stack(rcp45$RasterStack.Mean)
names(rcp45.mu)<- paste(rcp45$Year, rcp45$Month, sep = ".")
rcp45.5th<- raster::stack(rcp45$RasterStack.Pct05)
names(rcp45.5th)<- paste(rcp45$Year, rcp45$Month, sep = ".")
rcp45.95th<- raster::stack(rcp45$RasterStack.Pct95)
names(rcp45.95th)<- paste(rcp45$Year, rcp45$Month, sep = ".")

rcp85.mu.df<- as.data.frame(rcp85.mu, xy = TRUE) %>%
  gather(., "Year", "SST.Mean", -x, -y)
rcp85.5th.df<- as.data.frame(rcp85.5th, xy = TRUE) %>%
  gather(., "Year", "SST.pct5th", -x, -y)
rcp85.95th.df<- as.data.frame(rcp85.95th, xy = TRUE) %>%
  gather(., "Year", "SST.pct95th", -x, -y)

rcp85.df.all<- rcp85.mu.df
rcp85.df.all$SST.pct5th<- rcp85.5th.df$SST.pct5th
rcp85.df.all$SST.pct95th<- rcp85.95th.df$SST.pct95th

rcp85.df.all$Year<- gsub("X", "", gsub("[.]", "-", rcp85.df.all$Year))
rcp85.df.all<- rcp85.df.all %>%
  separate(Year, into = c("Year", "Month", "Day")) %>%
  group_by(Year, Month) %>%
  summarize(NELME.SST.Mean = mean(SST.Mean, na.rm = TRUE),
            NELME.SST.pct5th = mean(SST.pct5th, na.rm = TRUE),
            NELME.SST.pct95th = mean(SST.pct95th, na.rm = TRUE))
rcp85.df.all$Scenario<- rep("RCP85", nrow(rcp85.df.all))

rcp45.mu.df<- as.data.frame(rcp45.mu, xy = TRUE) %>%
  gather(., "Year", "SST.Mean", -x, -y)
rcp45.5th.df<- as.data.frame(rcp45.5th, xy = TRUE) %>%
  gather(., "Year", "SST.pct5th", -x, -y)
rcp45.95th.df<- as.data.frame(rcp45.95th, xy = TRUE) %>%
  gather(., "Year", "SST.pct95th", -x, -y)

rcp45.df.all<- rcp45.mu.df
rcp45.df.all$SST.pct5th<- rcp45.5th.df$SST.pct5th
rcp45.df.all$SST.pct95th<- rcp45.95th.df$SST.pct95th

rcp45.df.all$Year<- gsub("X", "", gsub("[.]", "-", rcp45.df.all$Year))
rcp45.df.all<- rcp45.df.all %>%
  separate(Year, into = c("Year", "Month", "Day")) %>%
  group_by(Year, Month) %>%
  summarize(NELME.SST.Mean = mean(SST.Mean, na.rm = TRUE),
            NELME.SST.pct5th = mean(SST.pct5th, na.rm = TRUE),
            NELME.SST.pct95th = mean(SST.pct95th, na.rm = TRUE))
rcp45.df.all$Scenario<- rep("RCP45", nrow(rcp45.df.all))

rcp.all<- bind_rows(rcp85.df.all, rcp45.df.all)
rcp.all$Scenario<- factor(rcp.all$Scenario, levels = c("RCP85", "RCP45"))

# Plot, raw temps
rcp.all$Plot.Date<- as.Date(paste(rcp.all$Year, rcp.all$Month, "15", sep = "-"))
rcp.all<- rcp.all %>%
  dplyr::filter(., Plot.Date < "2100-01-15")

rcp.plot<- rcp.all %>%
  group_by(Year, Scenario) %>%
  summarize(NELME.SST.YrMean = mean(NELME.SST.Mean, na.rm = TRUE),
            NELME.SST.Yrpct5th = mean(NELME.SST.pct5th, na.rm = TRUE),
            NELME.SST.Yrpct95th = mean(NELME.SST.pct95th, na.rm = TRUE))

rcp.plot.out<- ggplot() + 
  geom_line(data = rcp.plot, aes(x = as.numeric(Year), y = NELME.SST.YrMean, color = Scenario)) +
  geom_ribbon(data = rcp.plot, aes(x = as.numeric(Year), ymin = NELME.SST.Yrpct5th, ymax = NELME.SST.Yrpct95th, fill = Scenario), alpha=0.15) +
  scale_color_manual(name = "Scenario", values = c("#e31a1c", "#1f78b4")) +
  scale_fill_manual(name = "Scenario", values = c("#e31a1c", "#1f78b4")) +
  labs(x = "Year", y = "Ensemble Projected SST", fill = "Scenario") +
  ylim(c(10, 20)) +
  theme(strip.background = element_blank(), axis.text.x=element_text(angle=90, hjust=1)) 
rcp.plot.out

# NCAR
proj.wgs<- "+proj=longlat +ellps=WGS84 +datum=WGS84" 
ncar.mu<- raster::stack("~/Dropbox/Andrew/Work/GMRI/Projects/COCA/Data/ClimateData/climate.sst.proj.grd")
ncar.mu<- projectRaster(ncar.mu, crs = proj.wgs) 
ncar.mu<- mask(ncar.mu, nelme)
ncar.5th<- raster::stack("~/GitHub/COCA/Data/climate.sst.proj.pct05.grd")
ncar.5th<- projectRaster(ncar.5th, crs = proj.wgs) 
ncar.5th<- mask(ncar.5th, nelme)
ncar.95th<- raster::stack("~/GitHub/COCA/Data/climate.sst.proj.pct95.grd")
ncar.95th<- projectRaster(ncar.95th, crs = proj.wgs) 
ncar.95th<- mask(ncar.95th, nelme)

ncar.mu.df<- as.data.frame(ncar.mu, xy = TRUE) %>%
  gather(., "Year", "SST", -x, -y) 
ncar.5th.df<- as.data.frame(ncar.5th, xy = TRUE) %>%
  gather(., "Year", "SST.pct5th", -x, -y)
ncar.95th.df<- as.data.frame(ncar.95th, xy = TRUE) %>%
  gather(., "Year", "SST.pct95th", -x, -y)

ncar.all.df<- ncar.mu.df
ncar.all.df$SST.pct5th<- ncar.5th.df$SST.pct5th
ncar.all.df$SST.pct95th<- ncar.95th.df$SST.pct95th

ncar.all.df$Year<- gsub("X", "", gsub("[.]", "-", ncar.all.df$Year))
ncar.all.df<- ncar.all.df %>%
  separate(Year, into = c("Year", "Month", "Day")) %>%
  group_by(Year, Month) %>%
  dplyr::summarize(NELME.SST.Mean = mean(SST, na.rm = TRUE),
            NELME.SST.pct5th = mean(SST.pct5th, na.rm = TRUE),
            NELME.SST.pct95th = mean(SST.pct95th, na.rm = TRUE))

ncar.plot<- ncar.all.df %>%
  group_by(Year) %>%
  summarize(NELME.SST.YrMean = mean(NELME.SST.Mean, na.rm = TRUE),
            NELME.SST.Yrpct5th = mean(NELME.SST.pct5th, na.rm = TRUE),
            NELME.SST.Yrpct95th = mean(NELME.SST.pct95th, na.rm = TRUE))

ncar.plot.out<- ggplot() + 
  geom_line(data = ncar.plot, aes(x = as.numeric(Year), y = NELME.SST.YrMean, color = "#b30000")) +
  geom_ribbon(data = ncar.plot, aes(x = as.numeric(Year), ymin = NELME.SST.Yrpct5th, ymax = NELME.SST.Yrpct95th, fill = "#b30000"), alpha=0.15) +
  scale_color_manual(name = "Scenario", values = c("#b30000")) +
  scale_fill_manual(name = "Scenario", values = c("#b30000")) +
  labs(x = "Year", y = "Ensemble Projected SST", fill = "#b30000") +
  ylim(c(10, 20)) +
  theme(strip.background = element_blank(), axis.text.x=element_text(angle=90, hjust=1)) 
ncar.plot.out

both.out<- plot_grid(rcp.plot.out, ncar.plot.out, labels = c("Our Data", "NCAR Data"))

# Looks pretty good -- check differences?
rcp85.only<- rcp.plot %>%
  dplyr::filter(., Scenario == "RCP85")

rcp.vs.ncar<- data.frame("Year" = rcp85.only$Year, 
                         "Mean.Diff" = rcp85.only$NELME.SST.YrMean - ncar.plot$NELME.SST.YrMean,
                         "Pct5th.Diff" = rcp85.only$NELME.SST.Yrpct5th - ncar.plot$NELME.SST.Yrpct5th,
                         "Pct95th.Diff" = rcp85.only$NELME.SST.Yrpct95th - ncar.plot$NELME.SST.Yrpct95th)

diff.plot.out<- ggplot() + 
  geom_line(data = rcp.vs.ncar, aes(x = as.numeric(Year), y = Mean.Diff, color = "#636363")) +
  geom_ribbon(data = rcp.vs.ncar, aes(x = as.numeric(Year), ymin = Pct5th.Diff, ymax = Pct95th.Diff, fill = "#636363"), alpha=0.15) +
  scale_color_manual(name = "Scenario", values = c("#636363")) +
  scale_fill_manual(name = "Scenario", values = c("#636363")) +
  labs(x = "Year", y = "Ensemble Projected SST\nOur pull vs NCAR Data", fill = "#636363") +
  theme(strip.background = element_blank(), axis.text.x=element_text(angle=90, hjust=1)) 
diff.plot.out

