library(tidyverse)
library(stringr)
library(sf)

#need to select gages to include (700?)
#get rid of hysets (lots of USA gages, randomly sample after taht)
gages <- st_read("data/Caravan/shapefiles/combined.shp") %>%
  filter(!str_detect(gauge_id, "camels_")) #hysets has USA and Canada

topull <- "placeholder"
num_of_gages <- 0

while(num_of_gages < 700){
  #sample 1 random gage
  sampled <- sample(gages$gauge_id, 1) 
  dataset_folder <- str_extract(sampled, "^.*(?=(_))")
  
  #does gage have all data?
  #2009-10-02 to 2014-09-30
  check <- read_csv(paste0("data/Caravan/timeseries/csv/",
                  dataset_folder, "/",
                  sampled, ".csv")) %>%
          filter(date >= mdy("10-01-2009") &
                  date < mdy("10-01-2014"))
  
  #1825 is all days in 5 years
  if(length(check$streamflow) >= 1825){
    topull <- append(topull, sampled)
  }
  
  num_of_gages <- length(topull) - 1
}

topull <- topull[-1]
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

#create and output sampled shp file
gages_samp <- gages %>%
  filter(gauge_id %in% topull)

st_write(gages_samp, "data/gages700.shp", append = FALSE)

#gage information file
#HUC_02,GAGE_ID,GAGE_NAME,LAT,LONG,DRAINAGE_AREA_KM2,Elevation_m
sets <- unique(dataset_folder)

info <- read_csv(paste0("data/Caravan/attributes/", sets[1],
                 "/attributes_other_", sets[1], ".csv"))

for(z in 2:length(sets)){
  info <- info %>%
    bind_rows(read_csv(paste0("data/Caravan/attributes/", sets[z],
                    "/attributes_other_", sets[z], ".csv")))
}

filter(info, gauge_id %in% topull) %>%
  rename(GAGE_ID = gauge_id,
         GAGE_NAME = gauge_name,
         LAT = gauge_lat,
         LONG = gauge_lon, 
         DRAINAGE_AREA_KM2 = area) %>%
  write_csv("data/gage_information_caravan700.csv")
