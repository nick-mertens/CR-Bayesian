#Setting work Directory
working_directory = "set it to data file destination"
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

################################################################################
#### SET INITIAL VARIABLES #####################################################
################################################################################
df = read.csv("mixeddata052722.csv")

## NOTE: 
#### to run for any problem area, please change q19_2 throughout the code

################################################################################
#### BAYESIAN HIERARCHICAL BERNOULLI-BETA: FUNCTIONALIZED ######################
################################################################################
years = length(unique(df$MY))
mile_groups = length(unique(df$MileGrpATC))
min_year = min(df$MY)

## data preprocessing
#### convert categorical features to factors
df$MY = as.factor(df$MY)
df$MileGrpATC = as.factor(df$MileGrpATC)

data_filter_func = function(df, make){
  ## filter data for a single Make
  df = df[df$MakeName == make, ]
  df = df[!is.na(df$Age) & !is.na(df$q19_2), ]
  return(df)
}

stan_data_func = function(df, years, mile_groups){
  ## function to create stan data
  #### observations, their length, and a matrix that counts the number of
  #### observations per MMT:MY:MileGrpATC. For incomplete MMT, the matrix is
  #### filled with zeros
  y = df$q19_2
  age = df$Age
  N = length(y)
  n_mmt = length(unique(df$MMT))
  N_MY_MileGrpATC = df %>%
    group_by(MMT, MY, MileGrpATC) %>%
    summarise(count = n())
  N_MY_MileGrpATC <- cbind(list(1:nrow(N_MY_MileGrpATC)), N_MY_MileGrpATC) # ADDED TO AVOID ERROR WITH COMPLETE FUNCTION - NM 02/23/23
  N_MY_MileGrpATC = lapply(unique(N_MY_MileGrpATC$MMT), 
                           function(x) {
                             matrix(
                               complete(
                                 N_MY_MileGrpATC[N_MY_MileGrpATC$MMT == x, ], 
                                 MY, 
                                 MileGrpATC, 
                                 fill = list(count = 0)
                               )$count,
                               nrow = years, 
                               ncol = mile_groups, 
                               byrow = TRUE
                             )
                           })
  N_MY_MileGrpATC = do.call(cbind, N_MY_MileGrpATC)
  dim(N_MY_MileGrpATC) = c(3, 3, n_mmt)
  return(list(N = N,
              n_years = years,
              n_grps = mile_groups,
              n_mmt = n_mmt,
              N_MY_MileGrpATC = N_MY_MileGrpATC, 
              y = y,
              age = age))
}

max_density_func <- function(x) {
  ## function to calculate the mode of a distribution
  max_density_x = which.max(density(x)$y)
  return(density(x)$x[max_density_x])
}

extract_coef_func <- function(fit){
  ## extract Bayesian logistic coefficients from the object fit
  mcmc_df = as.data.frame(fit)
  Beta_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta[')]
  Beta_df = mcmc_df[, names(mcmc_df) %in% Beta_cols]
  Beta_values = colMeans(Beta_df)
  
  Beta0_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_0[')]
  Beta0_df = mcmc_df[, names(mcmc_df) %in% Beta0_cols]
  Beta0_values = colMeans(Beta0_df)
  
  return(list(b = Beta_values, b0 = Beta0_values))
}

predict_func <- function(stan_data, coef, df){
  ## predict new outcomes based on each observation's age
  y_new = c()
  beta_df = data.frame(matrix(ncol = 5, nrow = 0))
  colnames(beta_df) = c('MMT', 'MY', 'MileGrpATC', 'b0', 'b')
  
  n = 0
  for (i in 1:stan_data$n_mmt){
    for (j in 1:stan_data$n_years){
      for (k in 1:stan_data$n_grps){
        if (stan_data$N_MY_MileGrpATC[j, k, i] > 0){
          for (t in 1:stan_data$N_MY_MileGrpATC[j, k, i]){
            b_name = paste('Beta[', j, ',', k, ',', i, ']', sep = "")
            b0_name = paste('Beta_0[', j, ',', k, ',', i, ']', sep = "")
            Beta = coef$b[b_name]
            Beta0 = coef$b0[b0_name]
            predict = 1 / (1 + exp(-(Beta0 + Beta * stan_data$age[t + n])))
            y_new = c(y_new, predict)
          }
          mmt = unique(df$MMT)[i]
          my = unique(df[df$MMT == mmt, ]$MY)[j]
          mpg = unique(df[df$MMT == mmt & df$MY == my, ]$MileGrpATC)[k]
          beta_df[nrow(beta_df) + 1, ] = c(mmt, my, mpg, Beta0, Beta)
          n = n + stan_data$N_MY_MileGrpATC[j, k, i]
        }
      }
    }
  }
  return(list(y_new = y_new, beta_df = beta_df))
  ## TO DO:
  #### 1- This process is not the most accurate
  #### 2- Use aposteriori samples for prediction
}

## write the stan hierarchical model
#### the following stan model assumes a certain hierarchy where
#### y^t_ijk | Beta_0_ijk, Beta_ijk ~ Bernoulli(Beta_0_ijk + Beta_ijk * age^t_ijk),
#### Beta_0_ijk ~ Normal(0, 1 / kappa),
#### Beta_ijk | omega_ijk ~ Normal(omega_ijk, 1 / kappa)
#### omega_ijk | alpha_ij ~ Normal(alpha_ij, 1 / kappa),
#### alpha_ij | rho_i ~ Beta(rho_i, 1 / kappa),
#### rho_i | mu ~ Normal (mu, 1 / kappa)
#### mu ~ Normal (0, 1)
#### kappa ~ Gamma (1, 1)
#### i in MMT, j in MMT:MY, k in MMT:MY:MileGrpATC

## NOTE:
#### the model is built for a single Make but will work dynamically with all Makes
write(
  "//Bernoulli hierarchical model
data {
  int <lower=0> N;                                                //number of sample size per MMT
  int <lower=0> n_years;                                          //number of unique MY
  int <lower=0> n_grps;                                           //number of unique MileGrpATC
  int <lower=0> n_mmt;                                            //number of unique MMT per MakeName
  int <lower=0> N_MY_MileGrpATC[n_years, n_grps, n_mmt];          //number of observations per MY, MileGrpATC
  int <lower=0, upper=1> y[N];                                    //outcome
  int <lower=0> age[N];                                           //age
}

parameters {
  real <lower=0> kappa;                      //std for all levels of hierarchy
  real mu;                                   //Make-level mean parametr for normal priors
  real rho[n_mmt];                           //MMT-level mean parameter for normal priors
  real alpha[n_mmt, n_years];                //My-level mean parameter for normal priors
  real omega[n_years, n_grps, n_mmt];        //MileGrpATC-level mean parameter for normal prior
  real Beta[n_years, n_grps, n_mmt];         //Logistic regression coefficient for age at Make-MMT-MY-MileGrpATC level
  real Beta_0[n_years, n_grps, n_mmt];       //Intercept for the logistic regression at Make-MMT-MY-MileGrpATC level
}

model {
  int n = 0;                                                                                        //to keep track of number of observation per MMT:MY:MileGrpATC
  kappa ~ gamma(1, 1);                                                                              //Make-level prior uncertainty parameter
  mu ~ normal(0, 1);                                                                                //Make-level prior shape paramet
    for (i in 1:n_mmt){                                                                             //MMT-level prior shape parameter
    rho[i] ~ normal(mu, 1 / kappa);                                                                 //MMT-level prior shape parameter
    for (j in 1:n_years){                                                                           //loop over MY for MY- and MY:MileGrpATC-level priors                                                             
      if (sum(N_MY_MileGrpATC[j, , i]) > 0){                                                        //assign priors for MY- and MY:MileGrpATC-level parameters only if there is at least one observation in that MY
        alpha[i, j] ~ normal(rho[i], 1/ kappa);                                                     //MY-level prior shape parameter
        for (k in 1:n_grps){                                                                        //loop over MY:MileGrpATC-level priors
          if (N_MY_MileGrpATC[j, k, i] > 0){                                                        //assign priors for over MY:MileGrpATC-level parameters only if there is at least one observation per MY:MileGrpATC
            omega[j, k, i] ~ normal(alpha[i, j], 1 / kappa);                                        //MMT:MY:MileGrpATC-level prior shape parameter
            Beta[j, k, i] ~ normal(omega[j, k, i], 1 / kappa);                                      //MMT:MY:MileGrpATC-level prior problem rate
            Beta_0[j, k, i] ~ normal(0, 1 / kappa);
            for (t in 1:N_MY_MileGrpATC[j, k, i]){                                                  //loop over vehicle observations per MMT:MY:MileGrpATC for model specification
              y[t + n] ~ bernoulli_logit(Beta_0[j, k, i] + Beta[j, k, i] * age[t + n]);             //vehicle observation is a function of problem rate per MMT:MY:MileGrpATC
            }
            n = n + N_MY_MileGrpATC[j, k, i];                                                       //to keep track of number of observation per MMT:MY:MileGrpATC 
          }
        }
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

## generate an empty data frame to store Bayesian results
bayes_coef = data.frame()
bayes_pred = data.frame()

## initialize an iterator to show progression
i = 0

## Bayesian hierarchical models are run for each MMT
#### 1- data is filtered for each MMT
#### 2- stan data for feeding to the stan model is created
#### 3- MCMC chains are run
#### 4- results are filtered and written to file

for (make in (unique(df$MakeName))){
  col.names = ifelse(i == 0, TRUE, FALSE) #to write col names in the final csv
  i = i + 1
  print(i)
  temp_df = data_filter_func(df, make)
  stan_data = stan_data_func(temp_df, years, mile_groups)
  fit <- stan(
    file = here::here("hier_LG.stan"),
    data = stan_data,
    iter = 20000,
    warmup = 1000,
    chains = detectCores(), # number of chains
    cores = detectCores(),
    thin = 10,
    init_r = 0,
    #control = list(adapt_delta = 0.99),
    seed = 1231
  )
  coef = extract_coef_func(fit)
  pred_res = predict_func(stan_data, coef, temp_df)
  
  temp_df['y_pred'] = pred_res$y_new
  temp_df = temp_df %>%
    select(MakeName, MMT, MY, MileGrpATC, q19_2, y_pred)
  
  ## write prediction results (MileGrpATC level probabilities) to a csv file
  #### csv and workspace save
  write.table(pred_res$beta_df, 
              file = "q19_2-BAYES_hier_LG_COEF.csv", 
              sep = ",",
              col.names = col.names,
              row.names = FALSE,
              append = TRUE)
  bayes_coef = rbind(bayes_coef, pred_res$beta_df)
  write.table(temp_df, 
              file = "q19_2-BAYES_hier_LG_pred.csv", 
              sep = ",",
              col.names = col.names,
              row.names = FALSE,
              append = TRUE)
  bayes_pred = rbind(bayes_pred, temp_df)
  gc() #garbage removal to free memory
}

## write prediction results (MY level probabilities) to a csv file
#### evaluate cell probabilities for each MY
#### assuming all model years have equal weight`  `
prob_df = bayes_pred %>% 
  group_by(MMT, MY, MileGrpATC) %>%
  summarise(prob_rate = mean(as.numeric(y_pred)))
write.csv(prob_df, 
          paste("q19_2-BAYES_hier_LG",
                "_prob_rate.csv",
                sep = ""))

mean_df = bayes_pred %>% 
  group_by(MMT, MY) %>%
  summarise(MY_prob_rate = mean(as.numeric(y_pred)))
write.csv(prob_df, 
          paste("q19_2-BAYES_hier_LG",
                "_mean_rate.csv",
                sep = ""))