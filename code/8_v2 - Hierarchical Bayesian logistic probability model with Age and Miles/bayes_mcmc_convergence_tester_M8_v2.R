# Setting work Directory
model_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code/8_v2 - Hierarchical Bayesian logistic probability model with Age and Miles"
setwd(model_directory)

source("bayes_hier_LP_M8_v2.R")

# Import necessary packages
library(bayestestR)
library(rstanarm)
library(comprehenr)
library(ggplot2)
library("gridExtra")
library(ggpubr)

conv_test = function(fit_model_name, problem_area, show_params=FALSE, plot_chains=FALSE) {
  # Convergence test with Gelman-Rubin Statistics
  # Parameter 1: Specify name of the STAN fit model
  # Parameter 2: Whether to print out divergent parameters
  # Parameter 3: Whether to plot the divergent chains
  
  # Load the saved fit model
  loaded_fit <- readRDS(fit_model_name)
  
  # Create Stan data
  # Extract make name from fit name 
  make = str_split(fit_model_name, "_")[[1]][2]
  temp_df = data_filter_func(df, make, problem_area=problem_area)
  stan_data = stan_data_func(temp_df, years, problem_area=problem_area)
  
  # Get all the parameter names 
  params <- to_list(for (p in dimnames(as.array(loaded_fit))$parameters) if ((!grepl("lp_", p, fixed=TRUE)) & (!grepl("log_lik", p, fixed=TRUE))) p)
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
  
  rds_list = to_list(for (m in rds_list) if (grepl(MakeName, m, fixed=TRUE)) m)
  
  # if # of iterations is specified
  if (!is.null(iter)) {
    rds_list = to_list(for (m in rds_list) if (str_split(m, "_|\\.")[[1]][5] %in% iter) m)
  }
  
  # if # of chains is specified
  if (!is.null(chains)) {
    rds_list = to_list(for (m in rds_list) if (str_split(m, "_|\\.")[[1]][6] %in% chains) m)
  }
  
  # Parameters of interests: mu, kappa, rho, beta_0, beta_1, beta_2
  mu1_df = data.frame()
  mu2_df = data.frame()
  kappa1_df = data.frame()
  kappa2_df = data.frame()
  rho1_df = data.frame()
  rho2_df = data.frame()
  beta0_df <- data.frame()
  beta1_df <- data.frame()
  beta2_df <- data.frame()
  
  for (r in rds_list) {
    loaded_fit <- readRDS(paste("models/", r, sep=""))
    iter_num = str_split(r, "_|\\.")[[1]][5]
    chains_num = str_split(r, "_|\\.")[[1]][6]
    ds_num = str_split(r, "_|\\.")[[1]][9]
    
    #legend = sub('_NA', '_fs', paste(iter_num,"_",chains_num,"_",ds_num,sep=""))
    legend_ds = sub('NA', 'Full', paste(ds_num, "samples"))
    
    # Save parameter name as well as iterations & chains information 
    mu1_df <- rbind(mu1_df, data.frame(iter_chains=legend_ds, mu1=extract(loaded_fit)$mu[,1]))
    mu2_df <- rbind(mu2_df, data.frame(iter_chains=legend_ds, mu2=extract(loaded_fit)$mu[,2]))
    # kappa1_df <- rbind(kappa1_df, data.frame(iter_chains=legend_ds, kappa1=extract(loaded_fit)$kappa[,1]))
    # kappa2_df <- rbind(kappa2_df, data.frame(iter_chains=legend_ds, kappa2=extract(loaded_fit)$kappa[,2]))
    # rho1_df <- rbind(rho1_df, data.frame(iter_chains=legend_ds, rho1=extract(loaded_fit)$rho[,1,i]))
    # rho2_df <- rbind(rho2_df, data.frame(iter_chains=legend_ds, rho2=extract(loaded_fit)$rho[,2,i]))
    # beta0_df <- rbind(beta0_df, data.frame(iter_chains=legend_ds, beta0=extract(loaded_fit)$Beta_0[,j,i]))
    beta1_df <- rbind(beta1_df, data.frame(iter_chains=legend_ds, beta1=extract(loaded_fit)$Beta[,j,i,1]))
    beta2_df <- rbind(beta2_df, data.frame(iter_chains=legend_ds, beta2=extract(loaded_fit)$Beta[,j,i,2]))
  }
  
  size = 1
  axis_theme = theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 16), legend.text = element_text(size=14))
  
  # Plot posteriors for each parameter. Each overlay represents respective iteration & chains input for training
  mu1 <- ggplot(mu1_df, aes(x=mu1, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme 
  mu2 <- ggplot(mu2_df, aes(x=mu2, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  # kappa1 <- ggplot(kappa1_df, aes(x=kappa1, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  # kappa2 <- ggplot(kappa2_df, aes(x=kappa2, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  # rho1 <- ggplot(rho1_df, aes(x=rho1, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  # rho2 <- ggplot(rho2_df, aes(x=rho2, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  # beta0 <- ggplot(beta0_df, aes(x=beta0, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  beta1 <- ggplot(beta1_df, aes(x=beta1, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  beta2 <- ggplot(beta2_df, aes(x=beta2, color=iter_chains)) + geom_density(linewidth=size) + ylab(NULL) + axis_theme
  figure <- ggarrange(mu1, mu2, beta1, beta2, ncol=2, nrow=2, common.legend=TRUE) 
  annotate_figure(figure, fig.lab.size=16)
}

plot_posterior = function(model_path="models/", model_name, i, j) {
  # Plot posterior distributions of parameters for a single Stanfit object
  # Parameter 1: path where Stanfit objects are located
  # Parameter 2: RDS file name
  # Parameter 3: index of MMT
  # Parameter 4: index of MY
  
  # Parameters of interests: mu, kappa, rho, beta_0, beta_1, beta_2
  mu1_df = data.frame()
  mu2_df = data.frame()
  beta1_df <- data.frame()
  beta2_df <- data.frame()

  loaded_fit <- readRDS(paste("models/", model_name, sep=""))
  iter_num = str_split(model_name, "_|\\.")[[1]][5]
  chains_num = str_split(model_name, "_|\\.")[[1]][6]
  ds_num = str_split(model_name, "_|\\.")[[1]][8]
  
  # Save parameter name as well as iterations & chains information 
  mu1_df <- rbind(mu1_df, data.frame(iter_chains=sub('_NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), mu1=extract(loaded_fit)$mu[,1]))
  mu2_df <- rbind(mu2_df, data.frame(iter_chains=sub('_NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), mu2=extract(loaded_fit)$mu[,2]))
  beta1_df <- rbind(beta1_df, data.frame(iter_chains=sub('_NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), beta1=extract(loaded_fit)$Beta[,j,i,1]))
  beta2_df <- rbind(beta2_df, data.frame(iter_chains=sub('_NA', '', paste(iter_num,"_",chains_num,"_",ds_num,sep="")), beta2=extract(loaded_fit)$Beta[,j,i,2]))
  mu1_mode <- max_density_func(mu1_df$mu1)
  mu2_mode <- max_density_func(mu2_df$mu2)
  beta1_mode <- max_density_func(beta1_df$beta1)
  beta2_mode <- max_density_func(beta2_df$beta2)
  
  size = 1
  axis_theme = theme(axis.title.x = element_text(size = 32), axis.text.x = element_text(size = 24))
  
  # Plot posteriors for each parameter. Each overlay represents respective iteration & chains input for training
  mu1 <- ggplot(mu1_df, aes(x=mu1, color=iter_chains)) + geom_density(size=size) + ylab(NULL) + axis_theme + 
    geom_vline(xintercept=mu1_mode, size=1.5, color="red") + geom_text(aes(x=mu1_mode, label=paste("Mode:", round(mu1_mode, 3)), y=0.5, vjust=0.8, hjust=-0.2), colour='black', size=11)
  mu2 <- ggplot(mu2_df, aes(x=mu2, color=iter_chains)) + geom_density(size=size) + ylab(NULL) + axis_theme + 
    geom_vline(xintercept=mu2_mode, size=1.5, color="red") + geom_text(aes(x=mu2_mode, label=paste("Mode:", round(mu2_mode, 3)), y=0.5, vjust=0.8, hjust=-0.2), colour='black', size=11)
  beta1 <- ggplot(beta1_df, aes(x=beta1, color=iter_chains)) + geom_density(size=size) + ylab(NULL) + axis_theme + 
    geom_vline(xintercept=beta1_mode, size=1.5, color="red") + geom_text(aes(x=beta1_mode, label=paste("Mode:", round(beta1_mode, 3)), y=0.3, vjust=0.8, hjust=-0.2), colour='black', size=11)
  beta2 <- ggplot(beta2_df, aes(x=beta2, color=iter_chains)) + geom_density(size=size) + ylab(NULL) + axis_theme + 
    geom_vline(xintercept=beta2_mode, size=1.5, color="red") + geom_text(aes(x=beta2_mode, label=paste("Mode:", round(beta2_mode, 3)), y=0.3, vjust=0.8, hjust=-0.2), colour='black', size=11)
  figure <- ggarrange(mu1, mu2, beta1, beta2, ncol=2, nrow=2, common.legend=TRUE) 
  annotate_figure(figure, fig.lab.size=16)
}

#conv_test("models/fit_Acura_M8_v2_5000_12.rds", "q19_2")
#plot_posterior("models/","fit_Acura_M8_v2_5000_12.rds", 1, 1)
# Posteriors for Nissan Rogue, MY 2019
multiplot_posterior(MakeName="Nissan", i=9, j=2, downsampled=TRUE)
