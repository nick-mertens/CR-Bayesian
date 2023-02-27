# Setting work Directory
model_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian-LogisticProb/models"
setwd(model_directory)

# Import necessary packages
library(bayestestR)
library(rstanarm)
library(comprehenr)
library(ggplot2)

conv_test = function(model_name, show_params=FALSE, plot_chains=FALSE) {
  # Convergence test with Gelman-Rubin Statistics
  # Parameter 1: Specify name of the STAN fit model  
  # Parameter 2: Whether to print out divergent parameters
  # Parameter 3: Whether to plot the divergent chains
  
  # Load the saved fit model by Makes
  loaded_fit <- readRDS(model_name)
  
  # Gelman-Rubin Statistics
  # Get all the parameter names 
  params <- to_list(for (p in dimnames(as.array(loaded_fit))$parameters) if (!grepl("lp_", p, fixed=TRUE)) p)
  div_params <- list()
  
  for (p in params) {
    diag_ = gelman.diag(As.mcmc.list(loaded_fit, pars=p, include=TRUE), confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                        multivariate=TRUE)
    # If PSRF (Potential Scale Reduction Factor) is beyond normal range
    if ((diag_$psrf[1] >= 1.1)) {
      div_params <- append(div_params, p) 
    }
  }
  
  # print the number of parameters with divergence
  print(paste(length(div_params), "parameters failed in convergence test for", model_name))
  
  # if show_params is TRUE, list all parameters with divergence
  if (show_params == TRUE) {
    for (p in div_params) {
      print(p)
    }
  }
  
  # if plot_chains is TRUE, plot MCMC chains for divergent parameters
  if (plot_chains == TRUE) {
    if (length(div_params) > 0) {
      stan_trace(loaded_fit, pars=div_params, include=TRUE, inc_warmup=TRUE)
    }
  }
}

# Example Run
conv_test("fit_M8_Mini.rds", show_params=TRUE, plot_chains=TRUE)