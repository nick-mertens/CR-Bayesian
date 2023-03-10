# Setting work Directory
model_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code/8 - Hierarchical Bayesian logistic probability model with Age and Miles"
setwd(model_directory)

source("bayes_hier_LP_M8.R")

# Import necessary packages
library(bayestestR)
library(rstanarm)
library(comprehenr)
library(ggplot2)
library("gridExtra")
library(ggpubr)

conv_test = function(fit_model_name, show_params=FALSE, plot_chains=FALSE) {
  # Convergence test with Gelman-Rubin Statistics
  # Parameter 1: Specify name of the STAN fit model
  # Parameter 2: Whether to print out divergent parameters
  # Parameter 3: Whether to plot the divergent chains
  
  # Load the saved fit model
  loaded_fit <- readRDS(fit_model_name)
  
  # Create Stan data
  # Extract make name from fit name 
  make = str_split(fit_model_name, "_")[[1]][2]
  temp_df = data_filter_func(df, make)
  stan_data = stan_data_func(temp_df, years)
  
  # Get all the parameter names 
  params <- to_list(for (p in dimnames(as.array(loaded_fit))$parameters) if (!grepl("lp_", p, fixed=TRUE)) p)
  div_params <- list()
  
  for (p in params) {
    # Check whether parameters can be skipped from checking
    if (grepl("alpha", p) | grepl("Beta", p)) {
      indices = str_split(str_extract_all(p, pattern="(\\d+),(\\d+)")[[1]], ",")[[1]]
      if (grepl("alpha", p)) {
        i = strtoi(indices[1])
        j = strtoi(indices[2])
      } else {
        j = strtoi(indices[1])
        i = strtoi(indices[2])
      }
      # If that specific MMT-MY pair is non-existent, no need to test
      if (stan_data$N_MY[j,i] == 0) {
        next
      }
    }
    
    # Gelman-Rubin Statistics
    diag_ = gelman.diag(As.mcmc.list(loaded_fit, pars=p, include=TRUE), confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                        multivariate=TRUE)
    
    # If PSRF (Potential Scale Reduction Factor) is beyond normal range
    if ((diag_$psrf[1] >= 1.1)) {
      div_params <- append(div_params, p) 
    }
  }
  
  # print the number of parameters with divergence
  print(paste(length(div_params), "parameters failed in convergence test for", fit_model_name))
  
  # if show_params is TRUE, list all divergent parameters 
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
  return(loaded_fit)
}

multiplot_posterior = function(model_path="models/", MakeName, i, j, iter=NULL, chains=NULL, downsampled=TRUE) {
  # Plot posterior distributions of parameters for Stanfit objects for comparison
  # Parameter 1: path where Stanfit objects are located
  # Parameter 2: Name of the Make
  # Parameter 3: index of MMT
  # Parameter 4: index of MY
  # Parameter 5: Filter models by number of iterations (List)
  # Parameter 6: Filter models by number of chains (List)
  # Parameter 7: Whether models are downsampled
  
  # Collect all Stanfit objects of the input Make
  rds_list = list.files("models/")
  
  # Whether to look at downsampled models 
  if (downsampled) {
    rds_list = to_list(for (m in rds_list) if (grepl(MakeName, m, fixed=TRUE) & (grepl("_ds_", m, fixed=TRUE))) m)
  } else {
    rds_list = to_list(for (m in rds_list) if (grepl(MakeName, m, fixed=TRUE) & !(grepl("_ds_", m, fixed=TRUE))) m)
  }
    
  if (!is.null(iter)) {
    rds_list = to_list(for (m in rds_list) if (str_split(m, "_|\\.")[[1]][4] %in% iter) m)
  }
  
  if (!is.null(chains)) {
    rds_list = to_list(for (m in rds_list) if (str_split(m, "_|\\.")[[1]][5] %in% chains) m)
  }
  
  # Parameters of interests: mu, kappa, rho, beta_0, beta_1, beta_2
  mu_df = data.frame()
  kappa_df = data.frame()
  rho_df = data.frame()
  beta0_df <- data.frame()
  beta1_df <- data.frame()
  beta2_df <- data.frame()
  
  for (r in rds_list) {
    loaded_fit <- readRDS(paste("models/", r, sep=""))
    iter_num = str_split(r, "_|\\.")[[1]][4]
    chains_num = str_split(r, "_|\\.")[[1]][5]
    ds_num = str_split(r, "_|\\.")[[1]][8]
    
    # Save parameter name as well as iterations & chains information 
    mu_df <- rbind(mu_df, data.frame(iter_chains=sub('NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), mu=as.data.frame(loaded_fit)$mu))
    kappa_df <- rbind(kappa_df, data.frame(iter_chains=sub('NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), kappa=as.data.frame(loaded_fit)$kappa))
    rho_df <- rbind(rho_df, data.frame(iter_chains=sub('NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), rho=extract(loaded_fit)$rho[,1]))
    beta0_df <- rbind(beta0_df, data.frame(iter_chains=sub('NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), beta0=extract(loaded_fit)$Beta_0[,j,i]))
    beta1_df <- rbind(beta1_df, data.frame(iter_chains=sub('NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), beta1=extract(loaded_fit)$Beta[,j,i,1]))
    beta2_df <- rbind(beta2_df, data.frame(iter_chains=sub('NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), beta2=extract(loaded_fit)$Beta[,j,i,2]))
  }
  
  # Plot posteriors for each parameter. Each overlay represents respective iteration & chains input for training
  mu <- ggplot(mu_df, aes(x=mu, color=iter_chains)) + geom_density()
  kappa <- ggplot(kappa_df, aes(x=kappa, color=iter_chains)) + geom_density()
  rho <- ggplot(rho_df, aes(x=rho, color=iter_chains)) + geom_density()
  beta0 <- ggplot(beta0_df, aes(x=beta0, color=iter_chains)) + geom_density()
  beta1 <- ggplot(beta1_df, aes(x=beta1, color=iter_chains)) + geom_density()
  beta2 <- ggplot(beta2_df, aes(x=beta2, color=iter_chains)) + geom_density()
  ggarrange(mu, kappa, rho, beta0, beta1, beta2, ncol=3, nrow=2, common.legend=TRUE)
}

plot_posterior = function(model_path="models/", model_name, i, j) {
  # Plot posterior distribution of specified Stanfit objects
  # Parameter 1: path where Stanfit objects are located
  # Parameter 2: RDS file name
  # Parameter 3: index of MMT
  # Parameter 4: index of MY
  
  # Parameters of interests: beta_1, beta_2
  beta1_df <- data.frame()
  beta2_df <- data.frame()
  
  loaded_fit <- readRDS(paste("models/", model_name, sep=""))
  iter_num = str_split(model_name, "_|\\.")[[1]][4]
  chains_num = str_split(model_name, "_|\\.")[[1]][5]
  beta1_df <- rbind(beta1_df, data.frame(iter_chains=sub('_NA', '', paste(iter_num,"_",chains_num,sep="")), beta1=extract(loaded_fit)$Beta[,j,i,1]))
  beta2_df <- rbind(beta2_df, data.frame(iter_chains=sub('_NA', '', paste(iter_num,"_",chains_num,sep="")), beta2=extract(loaded_fit)$Beta[,j,i,2]))
  # Save parameter name as well as iterations & chains information 
  beta1_mode <- max_density_func(beta1_df$beta1)
  beta2_mode <- max_density_func(beta2_df$beta2)
  
  size = 1
  axis_theme = theme(axis.title.x = element_text(size = 40), axis.text.x = element_text(size = 28))
  # Plot posteriors for each parameter. Each overlay represents respective iteration & chains input for training
  beta1 <- ggplot(beta1_df, aes(x=beta1, color=iter_chains)) + geom_density(size=size) + ylab(NULL) + axis_theme + 
    geom_vline(xintercept=beta1_mode, size=1.5, color="red") + geom_text(aes(x=beta1_mode, label=paste("Mode:", round(beta1_mode, 3)), y=0.3, vjust=0.8, hjust=-0.2), colour='black', size=14)
  beta2 <- ggplot(beta2_df, aes(x=beta2, color=iter_chains)) + geom_density(size=size) + ylab(NULL) + axis_theme + 
    geom_vline(xintercept=beta2_mode, size=1.5, color="red") + geom_text(aes(x=beta2_mode, label=paste("Mode:", round(beta2_mode, 3)), y=0.3, vjust=0.8, hjust=-0.2), colour='black', size=14)
  figure <- ggarrange(beta1, beta2, ncol=2, nrow=1, common.legend=TRUE) 
  annotate_figure(figure, fig.lab.size=16)
}

#plot_posterior("models/","fit_Acura_M8_5000_12.rds", 1, 1)

