library(tidyverse)

source("app/HBV.R")

#Get data, relatively uniform distribution of precip, minimal snow
#camels_02017500,"JOHNS CREEK AT NEW CASTLE, VA",United States of America,37.50624,-80.10671,285.2584061790841

#37, -73 degrees

latitude <- 37
original <- read_csv("data/Caravan/timeseries/csv/camels/camels_02017500.csv") %>%
  filter(date >= mdy("10-01-2009") &
           date < mdy("10-01-2014")) %>%
  select(DATE = date,
         Q_mmd = streamflow, 
         prcp_mmd = total_precipitation_sum,
         T_mean_C = temperature_2m_mean,
         swe_mm = snow_depth_water_equivalent_mean) 

#calculate PET
latrad <- (latitude/360) * 2 * pi #convert to radians

PET1 <- original %>% select(DATE) %>%
  mutate(DOY = yday(DATE)) %>% #DOY for dates
  mutate(tempvar = (2 * pi / 365) * DOY) %>%
  mutate(delta_h = 0.4093 * sin(tempvar - 1.405)) %>% # declination of the sun above the celestial equator in radians on day JulDay of the year#
  mutate(daylen = (2 * acos(-tan(delta_h) * tan(latrad)) / 0.2618)) %>% #day length, h
  mutate(PET = 29.8 * daylen * 0.611 * exp(17.3 * original$T_mean_C / (original$T_mean_C + 237.3)) / (original$T_mean_C + 273.2))  #PET Hamon method

PET <- PET1$PET

#generate random parameters
rando <- function(){
  FC    <- runif(1, min = 40   , max = 400)  #Max soil moisture storage, field capacity
  beta  <- runif(1, min = 1    , max = 6)    #Shape coefficient governing fate of water input to soil moisture storage
  LP    <- runif(1, min = 0.3  , max = 1)    #Threshold for reduction of evap
  SFCF  <- runif(1, min = 0.4  , max = 1.2)  #Snowfall correction factor
  TT    <- runif(1, min = -1.5 , max = 1.2)  #Threshold temperature
  CFMAX <- runif(1, min = 1    , max = 8)    #Degree-day factor
  k0    <- runif(1, min = 0.05 , max = 0.5)  #Recession constant (upper storage, near surface)
  k1    <- runif(1, min = 0.01 , max = 0.3)  #Recession constant (upper storage)
  k2    <- runif(1, min = 0.001, max = 0.15) #Recession constant (lower storage)
  UZL   <- runif(1, min = 0    , max = 70)   #Threshold for shallow storage
  PERC  <- runif(1, min = 0    , max = 4)    #Percolation, max flow from upper to lower storage
  MAXBAS<- rep(1)
  
  c(FC, beta, LP, SFCF, TT, CFMAX, k0, k1, k2, UZL, PERC, MAXBAS)
}

#run HBV with random parameters and original data
output <- HBV(rando(), original$prcp_mmd, original$T_mean_C, PET, routing = 0)

#replace discharge and snow in original data with modeled discharge
HBVland <- original
HBVland$Q_mmd <- output$q
HBVland$swe_mm <- output$SWE

write_csv(HBVland, "app/Caravan700/HBVland1.csv")

attributes <- read_csv("app/gage_information_caravan700.csv")
attributes <- rbind(attributes,
  c("HBVland1", "HBV Land", "HBV Land", 37, -73, 100))

write_csv(attributes, "app/gage_information_caravan700.csv")
