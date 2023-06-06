#Run HBV, return just negative NSE for optimization routine
hbvnse <- function(params, dat, obs, routing){
  
  results <- HBV(params, dat$Precip, dat$Temp, dat$PET, routing = 0)
  
  results <- bind_cols(results, Qobs = obs)
  
  EvalStart <- floor(length(results$Qobs) * 0.4)
  EvalEnd <- length(results$Qobs)
  
  #trim the first 40% of the record so it isn't included in the NSE calculation
  results <- results[EvalStart:EvalEnd,]
  
  NSE <- 1 - ((sum((results$Qobs - results$q)^2, na.rm = TRUE))/
                sum((results$Qobs - mean(results$Qobs)) ^ 2, na.rm = TRUE))
  
  return(-NSE)
}
