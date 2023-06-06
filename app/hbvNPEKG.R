#Run HBV, return just negative NSE for optimization routine
#non-parametric version of kling-gupta from Pool et al 2018
source("NPE-KlingGupta.R")

hbvnpe <- function(params, dat, obs, routing){
  
  results <- HBV(params, dat$Precip, dat$Temp, dat$PET, routing = 0)
  
  results <- bind_cols(results, Qobs = obs)
  
  EvalStart <- floor(length(results$Qobs) * 0.4)
  EvalEnd <- length(results$Qobs)
  
  #trim the first 40% of the record so it isn't included in the NSE calculation
  results <- results[EvalStart:EvalEnd,]
  
  NPE <- RNP(sim = results$q,
            obs = results$Qobs)
  
  return(NPE)
}
