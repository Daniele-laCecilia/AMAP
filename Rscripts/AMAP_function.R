# Author: Daniele, la Cecilia (The University of Sydney @ daniele.lacecilia@sydney.edu.eu, now at eawag)
# Author: Giovanni M., Porta (Politecnico di Milano @ giovanni.porta@polimi.it)

# Should you use the AMAP function in your analysis, please cite the following publication
# Data and script were used in the journal article:
# Daniele la Cecilia, Giovanni M. Porta, Fiona H.M. Tang, Monica Riva, Federico Maggi
# Probabilistic indicators for soil and groundwater contamination risk assessment
# Ecological Indicators
# Volume 115
# 2020
# 106424
# https://doi.org/10.1016/j.ecolind.2020.106424
# http://www.sciencedirect.com/science/article/pii/S1470160X20303617

# Should you use the AMA family (AMAE, AMAV, AMAskew, AMAkurt) in your analysis, please cite the following publication
# The AMA family of indicators was introduced in:
# Aronne Dell'Oca, Monica Riva, Alberto Guadagnini
# Moment-based metrics for global sensitivity analysis of hydrological systems
# Hydrology and Earth System Sciences
# Volume 21
# 2017
# 6219-6234
# https://doi.org/10.5194/hess-21-6219-2017
# https://www.hydrol-earth-syst-sci.net/21/6219/2017/

###############################################################################################
###############################################################################################
###############################################################################################

# This function computes AMA sensitivity indices and principal Sobol Indices (Si) and AMAP
# starting from a set of Monte Carlo realizations. It also gives
# conditional moments with respect to one parameter at a time up to order 4
# inputs:

# The function needs a matrix called: Sampling_points , containing MC realizations of the parameters [n simulations x n parameters]
#                                     Output_Mat: Full model evaluations [n simulations x n time points]
# The user defines: nclass , the number of classes used to compute conditional statistics
#                   thr: threshold value
# The script automatically get: name_var , the name of each output variable (text string)

AMAP_func <- function(sampling_points, Output_Mat, nclass,thr,name_var) {
  
  ntime <- ncol(Output_Mat)
  Npar <- ncol(sampling_points)
  Nsim <- nrow(sampling_points)
  
  # initialize
  lower <- matrix(0,nrow=Npar, ncol=nclass)
  upper <- matrix(0,nrow=Npar, ncol=nclass)
  cond_mean <- array(0,dim=c(Npar, nclass))
  cond_var <- array(0,dim=c(Npar, nclass))
  cond_skew <- array(0,dim=c(Npar, nclass))
  cond_kurt <- array(0,dim=c(Npar, nclass))
  
  number <- array(0,dim=c(Npar, nclass))
  cond_PThr <- array(0,dim=c(Npar, nclass))
  
  AMAP <- array(0,dim=c(Npar, ntime))
  AMAE <- array(0,dim=c(Npar, ntime))
  AMAV <- array(0,dim=c(Npar, ntime))
  Si <- array(0,dim=c(Npar, ntime))
  AMAskew <- array(0,dim=c(Npar, ntime))
  AMAkurt <- array(0,dim=c(Npar, ntime))
  
  for (itime in 1:ntime) { # we calculated AMAP over time, but can be done with respect ot other variables
    output <- Output_Mat[,itime]
    
    Pthr_unc <- length(which(output>thr))/Nsim # number of exceedances
    
    for (ipar in 1:Npar) {
      param_i <- sampling_points[,ipar]
      sort_par <- sort(param_i, decreasing = FALSE) # ascending order
      isort <- order(param_i, decreasing = FALSE) # ascending order
      
      for (iclass in 1:nclass) {

        lower[ipar,iclass] <- param_i[isort[round(Nsim*(iclass-1)/nclass)+1]]
        upper[ipar,iclass] <- param_i[isort[round(Nsim*(iclass)/nclass)]]

        j_sampling_points_class <- which(param_i>lower[ipar,iclass] & param_i<upper[ipar,iclass])
        output_class <- output[j_sampling_points_class]
        
        cond_mean[ipar, iclass] <- mean(output_class)
        cond_var[ipar, iclass] <- var(output_class)
        cond_skew[ipar, iclass] <- skewness(output_class)
        cond_kurt[ipar, iclass] <- kurtosis(output_class)
  
        number[ipar,iclass] = length(j_sampling_points_class)
        cond_PThr[ipar, iclass] <- length(which(output_class>thr))/number[ipar,iclass]

      }
      
      AMAP[ipar,itime] <- mean(abs(cond_PThr[ipar, ] - Pthr_unc))
      AMAE[ipar,itime] <- mean(abs(cond_mean[ipar, ] - mean(output))) / abs(mean(output))
      AMAV[ipar,itime] <- mean(abs(cond_var[ipar, ] - var(output))) / var(output)
      Si[ipar,itime] <- var(cond_mean[ipar, ]) / var(output)
      AMAskew[ipar,itime] <- mean(abs(cond_skew[ipar, ] - skewness(output))) / abs(skewness(output))
      AMAkurt[ipar,itime] <- mean(abs(cond_kurt[ipar, ] - kurtosis(output))) / kurtosis(output)

    } # ipar
    
  } # itime
  
  # save outputs in the environment
  AMAP <<- AMAP
  AMAE <<- AMAE
  AMAV <<- AMAV
  Si <<- Si
  AMAskew <<- AMAskew
  AMAkurt <<- AMAkurt
  
}
