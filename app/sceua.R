params <- c(93.56175
            ,3.648981
            ,0.4683411
            ,0.7189698
            ,0.348700873
            ,6.284367
            ,0.24336037
            ,0.10287941
            ,0.006582184
            ,24.7481881
            ,2.460605500,1)

dat <- read_csv("Camels_Data/gage01013500.csv")%>%
  rename(Precip = prcp_mmd, Qobs = Q_mmd, 
         Temp = T_mean_C, Snow = swe_mm)

latrad <- (45/360) * 2 * pi #convert to radians

PET1 <- dat %>% select(DATE) %>%
  mutate(DOY = yday(DATE)) %>% #DOY for dates
  mutate(tempvar = (2 * pi / 365) * DOY) %>%
  mutate(delta_h = 0.4093 * sin(tempvar - 1.405)) %>% # declination of the sun above the celestial equator in radians on day JulDay of the year#
  mutate(daylen = (2 * acos(-tan(delta_h) * tan(latrad)) / 0.2618)) %>% #day length, h
  mutate(PET = 29.8 * daylen * 0.611 * exp(17.3 * dat$Temp / (dat$Temp + 237.3)) / (dat$Temp + 273.2))  #PET Hamon method

dat <- bind_cols(dat,PET = PET1$PET)
obs <- dat$Qobs

source("HBV.R")
#function that returns just NSE from HBV run
hbvnse <- function(params, dat, obs, routing){
    
  results <- HBV(params, dat$Precip, dat$Temp, dat$PET, routing = 0)
    
    results <- bind_cols(results, Qobs = obs)
    
    EvalStart <- floor(length(results$Qobs) * 0.4)
    EvalEnd <- length(results$Qobs)
    
    #trim the first 40% of the record so it isn't included in the NSE calculation
    results <- results[EvalStart:EvalEnd,]
    
    NSE <- 1 - ((sum((results$Qobs - results$q)^2))/
         sum((results$Qobs - mean(results$Qobs)) ^ 2))
      
    return(-NSE)
}


hbvnse(params, dat, obs, routing = 0)

library(rtop)

source("hbvnse.R")

lowpar <- c(40,1,0.3,0.4,-1.5 ,1  ,0.05 ,0.01 ,0.001,0,0,1)   
highpar <- c(400 ,6   ,1   ,1.2 ,1.2 ,8   ,0.5 ,0.3 ,0.15,70, 4,1  )

sceuaout <- sceua(hbvnse, pars = params, lower = lowpar, upper = highpar, dat = dat, obs = obs, routing = 0)
