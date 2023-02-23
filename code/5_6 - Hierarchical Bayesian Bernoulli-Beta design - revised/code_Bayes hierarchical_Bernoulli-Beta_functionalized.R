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

data_filter_func = function(df, mmt){
  ## filter data for a single MMT
  df = df[df$MMT == mmt, ]
  return(df)
}

stan_data_func = function(df, years, mile_groups){
  ## function to create stan data
  #### observations, their length, and a matrix that counts the number of
  #### observations per MMT:MY:MileGrpATC. For incomplete MMT, the matrix is
  #### filled with zeros
  y = df$q19_2[!is.na(df$q19_2)]
  N = length(y)
  N_MY_MileGrpATC = df %>% 
    group_by(MY, MileGrpATC) %>% 
    summarise(count = n())
  N_MY_MileGrpATC = matrix(complete(N_MY_MileGrpATC, MY, MileGrpATC, fill = list(count = 0))$count, 
                           nrow = years, 
                           ncol = mile_groups, 
                           byrow = TRUE)
  return(list(mode = mean(df$q19_2, na.rm = TRUE),
              N = N, 
              n_years = years,
              n_grps = mile_groups,
              N_MY_MileGrpATC = N_MY_MileGrpATC, 
              y = y))
}

max_density_func <- function(x) {
  ## function to calculate the mode of a distribution
  max_density_x = which.max(density(x)$y)
  return(density(x)$x[max_density_x])
}

theta_values_func <- function(df, fit, mmt, stan_data){
  ## extract mode of theta; or problem rates
  #### return theta values for valid MY:MileGrpATC
  #### valid My:MileGrpATC are those for which number of observation > 0
  theta_values = tibble(
    j = sort(rep(1:years, mile_groups)),
    k = rep(1:mile_groups, years),
    .variable = str_c(expression("theta"), "_", k, "_", j),
    values = as.vector(matrix(apply(matrix(extract(fit)$theta, 
                                           ncol = years * mile_groups, 
                                           byrow = FALSE), 
                                    2, 
                                    max_density_func), 
                              ncol = years * mile_groups, 
                              byrow = TRUE))
  ) %>% dplyr::select(.variable, values)
  temp_res = data.frame()
  for (year in 1:years){
    for (grp in 1:mile_groups){
      n = stan_data$N_MY_MileGrpATC[year, grp]
      if (n > 0){
        theta = theta_values$values[theta_values$.variable == paste('theta_', year, "_", grp, sep = "")]
        temp_res = rbind.data.frame(temp_res,
                                    cbind.data.frame(MMT = mmt,
                                                     MY = min_year + year - 1,
                                                     MileGrpATC = grp,
                                                     prob_rate = theta))
      }
    }
  }
  return(temp_res)
}
  
## write the stan hierarchical model
#### the following stan model assumes a certain hierarchy where
#### y^t_ijk | theta_ijk ~ Bernoulli(theta_ijk),
#### theta_ijk | omega_ijk ~ Beta(omega_ijk, kappa),
#### omega_ijk | alpha_ij ~ Beta(alpha_ij, kappa),
#### alpha_ij | rho_i ~ Beta(rho_i, kappa),
#### rho_i ~ Beta (mode, 100)
#### kappa ~ Gamma (1, 1)
#### i in MMT, j in MMT:MY, k in MMT:MY:MileGrpATC

## NOTE:
#### Beta distribution are parametrized with mode and concentration
#### the model is built for a single MMT but will work dynamically with all MMTs
write(
  "//Bernoulli hierarchical model
data {
  real <lower=0> mode;                                   //mode parameter
  int <lower=0> N;                                       //number of sample size per MMT
  int <lower=0> n_years;                                 //number of unique MY
  int <lower=0> n_grps;                                  //number of unique MileGrpATC
  int <lower=0> N_MY_MileGrpATC[n_years, n_grps];        //number of observations per MY, MileGrpATC
  int <lower=0, upper=1> y[N];                           //outcome
}

parameters {
  real <lower=2> kappa;                                  //MMT level concentration parameter
  real <lower=2> nu[n_years];                            //MMT-MY level concentration parameter
  real <lower=0, upper=1> rho;                           //MMT level mode parameter for Beta priors
  real <lower=0, upper=1> alpha[n_years];                //MMT-MY level mode parameter for Beta priors
  real <lower=0, upper=1> theta[n_years, n_grps];        //problem rate for MMT-MY-MileGrpATC level Bernoulli response
}

model {
  kappa ~ gamma(0.01, 0.01);                                                                        //MMT-level prior concentration parameter
  rho ~ beta(mode * (100 - 2) + 1, (1 - mode) * (100 - 2) + 1);                                     //MMT-level prior informed by MMT-level mode
  int n = 0;                                                                                        //to keep track of number of observation per MMT:MY:MileGrpATC
  for (j in 1:n_years){                                                                             //loop over MY for MY- and MY:MileGrpATC-level priors                                                             
    if (sum(N_MY_MileGrpATC[j, ]) > 0){                                                             //assign priors for MY- and MY:MileGrpATC-level parameters only if there is at least one observation in that MY
      alpha[j] ~ beta(rho * (kappa - 2) + 1, (1 - rho) * (kappa - 2) + 1);                          //MY-level prior mode parameter
      nu[j] ~ gamma(0.01, 0.01);                                                                    //MY-level prior concentration parameter
      for (k in 1:n_grps){                                                                          //loop over MY:MileGrpATC-level priors
        if (N_MY_MileGrpATC[j, k] > 0){                                                             //assign priors for over MY:MileGrpATC-level parameters only if there is at least one observation per MY:MileGrpATC
          theta[j, k] ~ beta(alpha[j] * (nu[j] - 2) + 1, (1 - alpha[j]) * (nu[j] - 2) + 1);         //MMT:MY:MileGrpATC-level prior problem rate
          for (t in 1:N_MY_MileGrpATC[j, k]){                                                       //loop over vehicle observations per MMT:MY:MileGrpATC for model specification
            y[t + n] ~ bernoulli(theta[j, k]);                                                      //vehicle observation is a function of problem rate per MMT:MY:MileGrpATC
          }
          n = n + N_MY_MileGrpATC[j, k];                                                            //to keep track of number of observation per MMT:MY:MileGrpATC 
        }
      }
    }
  }
}

", 
"hier_BB.stan")

## set Stan options for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

################################################################################
######## main loop starts here

## generate an empty data frame to store Bayesian results
bayes_res = data.frame()

## initialize an iterator to show progression
i = 0

## Bayesian hierarchical models are run for each MMT
#### 1- data is filtered for each MMT
#### 2- stan data for feeding to the stan model is created
#### 3- MCMC chains are run
#### 4- results are filtered and written to file

for (mmt in (unique(df$MMT))){
  col.names = ifelse(i == 0, TRUE, FALSE) #to write col names in the final csv
  i = i + 1
  print(i)
  temp_df = data_filter_func(df, mmt)
  stan_data = stan_data_func(temp_df, years, mile_groups)
  fit <- stan(
    file = here::here("hier_BB.stan"),
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
  temp_bayes_res = theta_values_func(temp_df, fit, mmt, stan_data)
  
  ## write prediction results (MileGrpATC level probabilities) to a csv file
  #### csv and workspace save
  write.table(temp_bayes_res, 
              file = "q19_2-BAYES_hier_BB_prob_rate.csv", 
              sep = ",",
              col.names = col.names,
              row.names = FALSE,
              append = TRUE)
  bayes_res = rbind(bayes_res, temp_bayes_res)
  gc() #garbage removal to free memory
}

## write prediction results (MY level probabilities) to a csv file
#### evaluate cell probabilities for each MY
#### assuming all model years have equal weight
mean_df = bayes_res %>% 
  group_by(MMT, MY) %>%
  summarise(MY_prob_rate = mean(as.numeric(prob_rate)))
write.csv(mean_df, 
          paste("q19_2-BAYES_hier_BB",
                "_mean_rate.csv",
                sep = ""))