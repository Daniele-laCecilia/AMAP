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

library(ggplot2)
rm(list=ls()) # clear environment

##############################################################################################################
### INPUT

# which indicator to be plotted? replace "AMAP" with any indicator you wish to plot
# indicators: AMAP, AMAE, AMAV, AMAskew, AMAkurt, Si
var.out <- "AMAP" # chosen by the user

dotime <- 1 # chosen by the user
# dotime <- 0: calculate indicators at a fixed time
# dotime <- 1: calculate indicators over time

## THIS IS YOUR MATRIX Output_mat (for the example without analysis over time)
# For simplicity, instead of having the whole time series for all simulations, we reported 2 cases:
file <- "Matrice_REF_C_GLP_15years.csv" # predicted concentrations in aquifer in microg/l at a fixed time = 15 years
# organized as: each row is a different simulation, each column is a different time
mat_fixed <- as.matrix(read.csv(file, header = TRUE, sep = ";", dec = ".", check.names = FALSE))
labout_fixed <- colnames(mat_fixed)

## OR THIS IS YOUR MATRIX Output_mat (for the example with analysis over time)
file <- "Matrice_REF_C_GLP_by5years.csv" # predicted concentrations in aquifer in microg/l over time
# In this example, predicted concentrations every 5 years
# organized as: each row is a different simulation, each column is a different time
mat_overtime <- as.matrix(read.csv(file, header = TRUE, sep = ";", dec = ".", check.names = FALSE))
labout_overtime <- colnames(mat_overtime)

## THIS IS YOUR MATRIX sampling_points
file <- "Matrice_parametri.csv" # values of uncertain parameters sampled from a uniform distribution using a Quasi Monte-Carlo technique
sampling_points <- as.matrix(read.csv(file, header = TRUE, sep=";", dec=".", check.names = FALSE))
# organized as: each row is a different simulation, each column is a different parameter
# names(sampling_points) <- c("k","b","psi_s","phi","Sr_l","Sr_g")
name_var <- colnames(sampling_points)

nclass <- 20 # to be chosen by the user, in theory, the more the better

thr <- 0.1 # to be chosen by the user
# In this example, thr corresponds to 0.1 microg/l, the safety threshold concentration for a concerning plant protection product in Europe
# For a mixture of concerning plant protection products and transformation products, thr <- 0.5
# For fluxes of concerning PPP, we used thr <- 0.2*1000 * 0.01/100; % mg/m2/year, that is 0.01% of a gross application rate of 2 kg/ha/year
# For fluxes of mixtures of concerning PPP and transformation product, we used thr <- 0.2*1000 * 0.05/100; % mg/m2/year, that is 0.05% of a gross application rate of 2 kg/ha/year given that only 1 active ingredient and 1 transformation product were accounted for in this analysis

### END INPUT
##############################################################################################################

if (dotime == 1) {
  Output_Mat <- mat_overtime # Output_Mat goes in the function to calculate the indicators
} else {
  Output_Mat <- mat_fixed
}



# Function not implemented in R. This implementation is identical to the formula in Matlab with flag=0, to give consistent results
source("skewness_function.R")

# Function not implemented in R. This implementation is identical to the formula in Matlab with flag=0, to give consistent results
source("kurtosis_function.R")

# Function developed to calculate the indicators
source("AMAP_function.R")


##############################################################################################################
### RUN FUNCTION
AMAP_func(sampling_points, Output_Mat, nclass,thr,name_var)

##############################################################################################################
### PLOT
Npar <- ncol(sampling_points)

if (dotime == 1) {
  tmp <- eval(parse(text = var.out))
  x_pos <- seq(1,ncol(tmp),1)
  tmp <- matrix(tmp, nrow=nrow(tmp)*ncol(tmp), ncol=1)
  df_plot <- as.data.frame(tmp)
  # for plot
  break.out <- seq(1,max(x_pos),1)
  lab.out <- labout_overtime
} else {
  tmp <- eval(parse(text = var.out))
  x_pos <- 1
  df_plot <- as.data.frame(tmp)
  # for plot
  break.out <- seq(1,max(x_pos),1)
  lab.out <- labout_fixed
}

df_plot <- cbind(as.data.frame(name_var) ,
                 rep(x_pos, each=Npar),
                 df_plot)
names(df_plot) <- c("Name_var","x_pos","Var1")

ggplot(data=df_plot, aes(x=x_pos, y=Var1, fill=Name_var)) +
  geom_bar(stat="identity") +
  theme(
    axis.text.x = element_text(size=18, hjust = 0.5, color="black", angle = 0),
    axis.text.y = element_text(size=18, hjust = 0.5, color="black", angle = 0),
    axis.title.x = element_text(size=18, hjust = 0.5, color="black", angle = 0),
    axis.title.y = element_text(size=18, hjust = 0.5, color="black", angle = 90),
    legend.title = element_text(color = "black", size = 18),
    legend.text = element_text(color = "black", size = 18)
  ) +
  scale_x_continuous(breaks=break.out,
                   labels=lab.out) +
  xlab("Year") +
  ylab(var.out) +
  labs(fill = "Parameters")
  
# namefolder <- ""
# namefile <- ""
# ggsave(paste0(namefolder, "/", namefile, ".jpg"), width = 10, height = 8, units = 'in')
  

