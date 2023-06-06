library(tidyverse)
library(stringr)
library(sf)

#need to select gages to include (700?)
#get rid of hysets (lots of USA gages, randomly sample after taht)
gages <- st_read("data/Caravan/shapefiles/combined.shp") %>%
  filter(!str_detect(gauge_id, "hysets")) #is this the only set that includes canada?

#choose gages that have the entire timeframe?
#2009-10-02 to 2014-09-30


#sample 700 random gages
sampled <- sample(gages$gauge_id, 700) 

gages_samp <- gages %>%
  filter(gauge_id %in% sampled)

st_write(gages_samp, "data/gages700.shp")

topull <- gages_samp$gauge_id
dataset_folder <- str_extract(topull, "^.*(?=(_))")

#trim to just the necessary date range and columns
#current app wants: DATE,Q_mmd,prcp_mmd,T_mean_C,swe_mm

for(i in 1:700){
  read_csv(paste0("data/Caravan/timeseries/csv/",
                       dataset_folder[i], "/",
                       topull[i], ".csv")) %>%
        select(DATE = date,
               Q_mmd = streamflow, 
               prcp_mmd = total_precipitation_sum,
               T_mean_C = temperature_2m_mean,
               swe_mm = snow_depth_water_equivalent_mean) %>%
        filter(DATE >= mdy("10-01-2009") &
                 DATE < mdy("10-01-2014")) %>%
        write_csv(paste0("data/Caravan700/",
                         topull[i], ".csv"))
}

#gage information file
#HUC_02,GAGE_ID,GAGE_NAME,LAT,LONG,DRAINAGE_AREA_KM2,Elevation_m
sets <- unique(folders)

info <- read_csv(paste0("data/Caravan/attributes/", sets[1],
                 "/attributes_other_", sets[1], ".csv"))

for(z in 2:length(sets)){
  info <- info %>%
    bind_rows(read_csv(paste0("data/Caravan/attributes/", sets[z],
                    "/attributes_other_", sets[z], ".csv")))
}

filter(info, gauge_id %in% gages_samp$gauge_id) %>%
  rename(GAGE_ID = gauge_id,
         GAGE_NAME = gauge_name,
         LAT = gauge_lat,
         LONG = gauge_lon, 
         DRAINAGE_AREA_KM2 = area) %>%
  write_csv("data/gage_information_caravan700.csv")
