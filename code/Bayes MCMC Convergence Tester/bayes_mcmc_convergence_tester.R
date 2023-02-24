# Setting work Directory
model_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian-LogisticProb"
setwd(model_directory)

# Import necessary packages
library(bayestestR)
library(rstanarm)
library(comprehenr)
library(ggplot2)

conv_test = function(file_path="../", show_params=FALSE) {
  # Convergence test with Gelman-Rubin Statistics
  # Specify file path where the STAN fit models are located
  
  # Collect all fit models by makes
  make_list <- list.files(path=file_path, pattern=".rds")
  
  # Dictionary to save all divergent parameters by makes 
  make_params = c()
  
  # Go through each make model to run diagnostics
  for (m in make_list) {
    # Load the saved fit model by Makes
    loaded_fit <- readRDS(paste(file_path, m, sep=""))
    
    # Monte Carlo Standard Errors 
    mc_res <- mcse(fit)
    
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
    
    # Store divergent parameters for each make 
    make_params[[m]] <- div_params
    
    # print the number of parameters with divergence
    print(paste(length(div_params), "parameters failed in convergence test for", m))
    
    # if show_params is TRUE, list all parameters with divergence
    if (show_params == TRUE) {
      for (p in div_params) {
        print(p)
      }
    }
  }
}


plot_mcmc = function(make_model) {
  # Plot MCMC Chains for divergent parameters 
  # Param 1: make_model: Stan Fit model name (Specify the stan file directory)
  loaded_fit <- readRDS(make_model)
  
  # Gelman-Rubin Statistics
  params <- to_list(for (p in dimnames(as.array(loaded_fit))$parameters) if (!grepl("lp_", p, fixed=TRUE)) p)
  div_params <- list()
  
  for (p in params) {
    diag_ = gelman.diag(As.mcmc.list(loaded_fit, pars=p, include=TRUE), confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                        multivariate=TRUE)
    if ((diag_$psrf[1] >= 1.1)) {
      div_params <- append(div_params,p) 
    }
  }

  if (length(div_params) < 1) {
    return(NULL)
  }
  stan_trace(loaded_fit, pars=div_params, include=TRUE, inc_warmup=TRUE)
}

# Example Run
conv_test("models/", show_params=TRUE)
plot_mcmc("models/fit_Nissan.rds")
