#Setting work Directory
working_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code/8_v2 - Hierarchical Bayesian logistic probability model with Age and Miles"
#working_directory = "/Users/nick.mertens/Library/CloudStorage/OneDrive-Blend360/Consumer Reports/Documents/2023 Bayesian Modeling - Phase II/CR-Bayesian/code/8 - Hierarchical Bayesian logistic probability model with Age and Miles"
setwd(working_directory)

#Importing required packages
library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(data.table)
library(coda)
library("lmtest")
library("sandwich") 
library("tidyverse")
library(lme4)
library(broom)
library(stringr)
library(brms)
library(cmdstanr)
library(rstan)
library(bayesplot)
library(tidybayes)
library(matrixStats)
library(TeachingDemos)
library(tidyr)
library(here)
library(parallel)
library(bruceR)
library(loo)

# Read original file into df
df = read.csv("mixeddata052722.csv")

# There are 3 unique years - 2018, 2019, 2020
years = length(unique(df$MY))
min_year = min(df$MY)

## data preprocessing
#### convert categorical features to factors
df$MY = as.factor(df$MY)

data_filter_func = function(df, make, problem_area){ # adding variable to allow for choice of problem area - NM 02/28/23
  ## filter data for a single Make
  df = df[df$MakeName == make, ]
  df = df[!is.na(df$Age) & !is.na(df[problem_area]) & !is.na(df$Miles), ]
  
  ## Standard Scaling Features
  df <- cbind(df, Age_scaled = scaler(df$Age, min = 0, max = 1))
  df <- cbind(df, Miles_scaled = scaler(df$Miles, min = 0, max = 1))
  
  return(df)
}

stan_data_func = function(df, years, problem_area){ # adding variable to allow for choice of problem area - NM 02/28/23
  ## create stan data
  #### observations, their length, and a matrix that counts the number of
  #### observations per MMT:MY. For incomplete MMT, the matrix is
  #### filled with zeros
  y = df[[problem_area]]
  age = df$Age_scaled  # assign scaled ages
  miles = df$Miles_scaled  # assign scaled mileages 
  N = length(y)
  n_mmt = length(unique(df$MMT))
  
  N_MY = df %>%
    group_by(MMT, MY) %>%
    summarize(count = n())
  
  # Add dummy column to avoid bug in complete() 
  N_MY = cbind(0, N_MY) 
  
  N_MY = lapply(unique(N_MY$MMT), 
                function(x) {
                  matrix(
                    complete(
                      N_MY[N_MY$MMT == x, ], 
                      MY, 
                      fill = list(count = 0)
                    )$count,
                    ncol = years,
                    byrow = TRUE
                  )
                })
  
  N_MY = do.call(cbind, N_MY)
  
  dim(N_MY) = c(3, n_mmt)
  
  return(list(N = N,
              n_years = years,
              n_mmt = n_mmt,
              N_MY = N_MY,
              y = y,
              age = age,
              miles = miles))
}

max_density_func <- function(x) {
  ## function to calculate the mode of a distribution
  max_density_x = which.max(density(x)$y)
  return(density(x)$x[max_density_x])
}

extract_coef_func <- function(fit, coef_mode=c("mode","mean")){
  ## extract Bayesian logistic coefficients from the object fit
  mcmc_df = as.data.frame(fit)
  
  Beta0_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_0[')]
  Beta0_df = mcmc_df[, names(mcmc_df) %in% Beta0_cols]
  
  Beta1_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta[') & endsWith(colnames(mcmc_df), '1]')]
  Beta1_df = mcmc_df[, names(mcmc_df) %in% Beta1_cols]
  
  Beta2_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta[') & endsWith(colnames(mcmc_df), '2]')]
  Beta2_df = mcmc_df[, names(mcmc_df) %in% Beta2_cols]
  
  # Coefficient selection mode
  if (coef_mode == "mode") {
    Beta0_values = apply(Beta0_df, MARGIN=2, FUN=max_density_func)
    Beta1_values = apply(Beta1_df, MARGIN=2, FUN=max_density_func)
    Beta2_values = apply(Beta2_df, MARGIN=2, FUN=max_density_func)
  } else {
    Beta0_values = colMeans(Beta0_df)
    Beta1_values = colMeans(Beta1_df)
    Beta2_values = colMeans(Beta2_df)
  }
  
  return(list(b2 = Beta2_values, b1 = Beta1_values, b0 = Beta0_values))
}

predict_func <- function(stan_data, coef, df){
  ## predict new outcomes based on each observation's age
  y_new = c()
  beta_df = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(beta_df) = c('MMT', 'MY', 'b0', 'b1', 'b2')
  
  n = 0
  for (i in 1:stan_data$n_mmt){
    for (j in 1:stan_data$n_years){
      if (stan_data$N_MY[j, i] > 0){
        for (t in 1:stan_data$N_MY[j, i]){
          b0_name = paste('Beta_0[', j, ',', i, ']', sep = "")
          b1_name = paste('Beta[', j, ',', i, ',1]', sep = "")
          b2_name = paste('Beta[', j, ',', i, ',2]', sep = "")
          Beta0 = coef$b0[b0_name]
          Beta1 = coef$b1[b1_name]
          Beta2 = coef$b2[b2_name]
          predict = 1 / (1 + exp(-(Beta0 + Beta1 * stan_data$age[t + n] + Beta2 * stan_data$miles[t + n])))
          y_new = c(y_new, predict)
        }
        mmt = unique(df$MMT)[i]
        my = unique(df[df$MMT == mmt, ]$MY)[j]
        beta_df[nrow(beta_df) + 1, ] = c(mmt, my, Beta0, Beta1, Beta2)
        n = n + stan_data$N_MY[j, i]
        
      }
    }
  }
  return(list(y_new = y_new, beta_df = beta_df))
  ## TO DO:
  #### 1- This process is not the most accurate
  #### 2- Use aposteriori samples for prediction
}

write(
  "//Bernoulli hierarchical model
data {
  int <lower=0> N;                                                //number of sample size per MMT
  int <lower=0> n_years;                                          //number of unique MY
  int <lower=0> n_mmt;                                            //number of unique MMT per MakeName
  int <lower=0> N_MY[n_years, n_mmt];          //number of observations per MY
  int <lower=0, upper=1> y[N];                                    //outcome
  real <lower=0, upper=1> age[N];                                           //age
  real <lower=0, upper=1> miles[N];                                         //Mileage
}

parameters {
  real <lower=0> kappa[2];                      //std for all levels of hierarchy
  real mu[2];                                   //Make-level mean parameter for normal priors
  real rho[2, n_mmt];                           //MMT-level mean parameter for normal priors
  real alpha[2, n_mmt, n_years];                //My-level mean parameter for normal priors
  real Beta_0[n_years, n_mmt];       //Intercept for the logistic regression at Make-MMT-MY level
  real Beta[n_years, n_mmt, 2];      //Coefficients for Age and Mileage at Make-MMT-MY level
}

model {
  int n = 0;                                                                                //** The int n initialization was moved to the first line. 
  kappa[1] ~ gamma(1, 1);
  kappa[2] ~ gamma(1, 1);
  mu[1] ~ normal(0, 1);
  mu[2] ~ normal(0, 1);
  for (i in 1:n_mmt){                                                                             
    rho[1,i] ~ normal(mu[1], 1 / kappa[1]);
    rho[2,i] ~ normal(mu[2], 1 / kappa[2]);
    for (j in 1:n_years){                                                                                                             
      if (N_MY[j,i] > 0){                                                      
        alpha[1,i,j] ~ normal(rho[1,i], 1/ kappa[1]);
        alpha[2,i,j] ~ normal(rho[2,i], 1/ kappa[2]);
        Beta_0[j, i] ~ normal(0, 1 / kappa[1]);
        Beta[j,i,1] ~ normal(alpha[1,i,j], 1 / kappa[1]);
        Beta[j,i,2] ~ normal(alpha[2,i,j], 1 / kappa[2]);
        for (t in 1:N_MY[j, i]){                                                  
          y[t + n] ~ bernoulli_logit(Beta_0[j, i] + Beta[j, i, 1] * age[t + n] + Beta[j, i, 2] * miles[t + n]);   
        }
        n = n + N_MY[j, i];                                                       
      }
    }
  }
}

generated quantities {
  // Declare a vector of length N to store the log likelihood for each observation
  vector[N] log_lik;
  
  {
    // Initialize a counter to keep track of the current observation number
    int n = 0;
    
    // Loop over each unique combination of i, j, and t
    for (i in 1:n_mmt) {
      for (j in 1:n_years) {
        if (N_MY[j,i] > 0) {
          for (t in 1:N_MY[j,i]) {
            // Calculate the log likelihood for the current observation using the
            // bernoulli_logit_lpmf function, which takes the observation y[t+n], and
            // the predicted log odds ratio Beta_0[j,i] + Beta[j,i,1]*age[t+n] + 
            // Beta[j,i,2]*miles[t+n].
            log_lik[t + n] = bernoulli_logit_lpmf(y[t + n] | Beta_0[j,i] + Beta[j,i,1]*age[t + n] + Beta[j,i, 2]*miles[t + n]);
          }
          // Increment the counter by the number of observations for the current 
          // combination of i and j
          n = n + N_MY[j,i];
        }
      }
    }
  }
}
", 
  "hier_LG_M8_v2.stan")

## set Stan options for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

################################################################################
######## main loop starts here

run_model <- function(df, make, iter=5000, chains=4, problem_area, save_fit=TRUE) { # adding variable to allow for choice of problem area - NM 02/28/23
  # df: Original DataFrame
  # iter: # of iterations to run MCMC
  # chains: # of MCMC chains
  # make_list: List of makes 
  # problem_area: Specific problem area
  # save_fit: whether to save the stanfit object
  
  ## generate an empty data frame to store Bayesian results
  #bayes_coef = data.frame()
  #bayes_pred = data.frame()
  
  ## Bayesian hierarchical models are run for each MMT
  #### 1- data is filtered for each MMT
  #### 2- stan data for feeding to the stan model is created
  #### 3- MCMC chains are run
  #### 4- results are filtered and written to file
  
  res_df = data.frame()
  temp_df = data_filter_func(df, make, problem_area)
  stan_data = stan_data_func(temp_df, years, problem_area)
  
  print(paste("Starting MCMC training for", make, "with", iter, "iterations with", chains, "chains"))
  
  fit <- stan(
    file = "hier_LG_M8_v2.stan",
    data = stan_data,
    iter = iter,
    warmup = 1000,
    chains = chains,
    cores = detectCores(),
    thin = 10,
    init_r = 1,
    #control = list(max_treedepth=10),
    seed = 1231,
    verbose = FALSE
  )
  
  # Save the Model if save_fit = TRUE
  if (save_fit) {
    saveRDS(fit, paste("models/fit_", make, "_M8_v2_", as.character(iter), "_", as.character(chains),".rds", sep=""))
  }
  
  return(fit)
  
}

pred_prob <- function(df, fit_model_name=NULL, problem_area, coef_mode=c("mode","mean")) {
  ## Compute naive and predicted probabilities by MMT-MY
  # Parameter 1: Original Dataframe
  # Parameter 2: Trained Stanfit model file - .rds file 
  # Parameter 3: Problem area of interest
  # Parameter 4: Coefficient selection method. "Mode" finds the mode of the posterior distribution. "Mean" finds the mean of the distribution.
  
  # Extract Make Name from fit model
  make = str_split(fit_model_name, "_")[[1]][2]
  
  # Filter data
  temp_df = data_filter_func(df, make, problem_area)
  years = length(unique(temp_df$MY))
  
  # create Stan Data 
  stan_data = stan_data_func(temp_df, years, problem_area)
  
  # Call Stan model if model name is given
  loaded_fit <- readRDS(fit_model_name)

  # Predict
  coef = extract_coef_func(loaded_fit, coef_mode)
  pred_res = predict_func(stan_data, coef, temp_df)
  
  temp_df['y_pred'] = pred_res$y_new
  
  temp_df = temp_df %>%
    select(MakeName, MMT, MY, problem_area, y_pred)
  
  # Compute mean cell probability & predicted probability 
  res_df = temp_df %>%
    group_by(MakeName, MMT, MY) %>% 
    dplyr::summarize(cnt=n(), across(everything(), list(mean=mean)))
  
  res_df = subset(res_df, select = -c(cnt_mean))
  res_df['%_deviation'] <- round((res_df$y_pred_mean - res_df$q19_2_mean) / res_df$q19_2_mean * 100, 1)
  res_df[c('q19_2_mean','y_pred_mean')] = round(res_df[c('q19_2_mean','y_pred_mean')], 4)
  
  return(res_df)
}

calculate_diagnostics <- function(filename){
  # Function to calculate diagnostics for model comparison
  # parameter 1: Trained Stanfit model file - .rds file
  
  # load in the model
  fit <- readRDS(filename)
  
  # extract the log likelihood values
  fit_log <- extract_log_lik(fit, "log_lik", merge_chains = FALSE)
  
  # calculate relative sample size
  fit_eff <- relative_eff(fit_log)
  
  # present WAIC
  print("WAIC")
  print(waic(fit_log))
  
  # present LOO
  print("Leave One Out CV")
  print(loo(fit_log, r_eff = fit_eff))
}

## Example Train
#for (iter in c(20000)) {
#  for (chain in c(12))
#  run_model(df, iter=iter, chains=chain, "Acura", "q19_2")
#}


#run_model(df, make="Nissan", iter=5000, chains=12,"q19_2", save_fit=TRUE)
# Example Calculate Diagnostics
#calculate_diagnostics("./models/fit_Acura_M8_v2_5000_12.rds")

# Example Compute Probability
#M8_v2_res_df = pred_prob(df, fit_model_name="models/fit_Nissan_M8_v2_5000_12.rds", problem_area = "q19_2", coef_mode="mode")

# View resulting table
#view(M8_v2_res_df)
