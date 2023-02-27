#Setting work Directory
working_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian-LogisticProb"
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

# Read original file into df
df = read.csv("../mixeddata052722.csv")

# There are 3 unique years - 2018, 2019, 2020
years = length(unique(df$MY))
min_year = min(df$MY)

## data preprocessing
#### convert categorical features to factors
df$MY = as.factor(df$MY)

data_filter_func = function(df, make){
  ## filter data for a single Make
  df = df[df$MakeName == make, ]
  df = df[!is.na(df$Age) & !is.na(df$q19_2) & !is.na(df$Miles), ]
  
  ## Standard Scaling Features
  df <- cbind(df, Age_scaled = scaler(df$Age, min = 0, max = 1))
  df <- cbind(df, Miles_scaled = scaler(df$Miles, min = 0, max = 1))
  
  return(df)
}

stan_data_func = function(df, years){
  ## Create stan data
  #### observations, their length, and a matrix that counts the number of
  #### observations per MMT:MY. For incomplete MMT, the matrix is
  #### filled with zeros
  y = df$q19_2
  age = df$Age_scaled
  miles = df$Miles_scaled
  N = length(y)
  n_mmt = length(unique(df$MMT))
  
  N_MY = df %>%
    group_by(MMT, MY) %>%
    summarise(count = n())
  
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

extract_coef_func <- function(fit){
  ## extract Bayesian logistic coefficients from the object fit
  mcmc_df = as.data.frame(fit)
  
  Beta0_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_0[')]
  Beta0_df = mcmc_df[, names(mcmc_df) %in% Beta0_cols]
  Beta0_values = apply(Beta0_df, MARGIN=2, FUN=max_density_func)
  
  Beta1_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_1[')]
  Beta1_df = mcmc_df[, names(mcmc_df) %in% Beta1_cols]
  Beta1_values = apply(Beta1_df, MARGIN=2, FUN=max_density_func)
  
  Beta2_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_2[')]
  Beta2_df = mcmc_df[, names(mcmc_df) %in% Beta2_cols]
  Beta2_values = apply(Beta2_df, MARGIN=2, FUN=max_density_func)
  
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
          b1_name = paste('Beta_1[', j, ',', i, ']', sep = "")
          b2_name = paste('Beta_2[', j, ',', i, ']', sep = "")
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
  real <lower=0> kappa;                      //std for all levels of hierarchy
  real mu;                                   //Make-level mean parameter for normal priors
  real rho[n_mmt];                           //MMT-level mean parameter for normal priors
  real alpha[n_mmt, n_years];                //My-level mean parameter for normal priors
  real Beta_0[n_years, n_mmt];       //Intercept for the logistic regression at Make-MMT-MY level
  real Beta_1[n_years, n_mmt];         //Logistic regression coefficient for age at Make-MMT-MY level
  real Beta_2[n_years, n_mmt];         //Logistic regression coefficient for miles at Make-MMT-MY level
}

model {
  int n = 0;    //** The int n initialization was moved to the first line. 
  kappa ~ gamma(1, 1);                                                        
  mu ~ normal(0, 1);                                                                                
  for (i in 1:n_mmt){                                                                             
    rho[i] ~ normal(mu, 1 / kappa);                                                                
    for (j in 1:n_years){                                                                                                             
      if (N_MY[j,i] > 0){                                                      
        alpha[i, j] ~ normal(rho[i], 1/ kappa);                                                     
        Beta_0[j, i] ~ normal(0, 1 / kappa);
        Beta_1[j, i] ~ normal(alpha[i, j], 1 / kappa);                                      //MMT:MY-level prior problem rate
        Beta_2[j, i] ~ normal(alpha[i, j], 1 / kappa);                                      //MMT:MY-level prior problem rate
  
        for (t in 1:N_MY[j, i]){                                                  
          y[t + n] ~ bernoulli_logit(Beta_0[j, i] + Beta_1[j, i] * age[t + n] + Beta_2[j, i] * miles[t + n]);            
        }
        n = n + N_MY[j, i];                                                       
      }
    }
  }
}

", 
  "hier_LG.stan")

## set Stan options for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

################################################################################
######## main loop starts here

run_model <- function(df, make_list) {
  ## generate an empty data frame to store Bayesian results
  #bayes_coef = data.frame()
  #bayes_pred = data.frame()
  
  ## initialize an iterator to show progression
  i = 0
  
  ## Bayesian hierarchical models are run for each MMT
  #### 1- data is filtered for each MMT
  #### 2- stan data for feeding to the stan model is created
  #### 3- MCMC chains are run
  #### 4- results are filtered and written to file
  
  # Number of iterations
  iter = 20000
  res_df = data.frame()
  
  for (make in (make_list)){
    temp_df = data_filter_func(df, make)
    stan_data = stan_data_func(temp_df, years)
    
    fit <- stan(
      file = "hier_LG.stan",
      data = stan_data,
      iter = iter,
      warmup = 1000,
      #chains = detectCores(), # number of chains
      chains = 8,
      cores = detectCores(),
      thin = 10,
      init_r = 0,
      #control = list(adapt_delta = 0.99),
      seed = 1231,
      verbose = FALSE
    )
    # Save the Model
    saveRDS(fit, paste("models/fit_M8_", as.character(iter), "_", make, ".rds", sep=""))
    
    coef = extract_coef_func(fit)
    pred_res = predict_func(stan_data, coef, temp_df)
    
    temp_df['y_pred'] = pred_res$y_new
    temp_df = temp_df %>%
      select(MakeName, MMT, MY, q19_2, y_pred)
    res_df = rbind(res_df, temp_df)
  }
}

# Example Run
run_model(df, c("Nissan"))