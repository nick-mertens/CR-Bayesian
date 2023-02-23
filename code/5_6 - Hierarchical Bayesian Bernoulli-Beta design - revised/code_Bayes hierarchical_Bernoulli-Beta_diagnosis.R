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
#### to diagnose a problem for a different MMT, change line 47

################################################################################
#### BAYESIAN HIERARCHICAL BERNOULLI-BETA: SINGLE MMT TRAINING #################
################################################################################
years = length(unique(df$MY))
mile_groups = length(unique(df$MileGrpATC))

## data preprocessing
#### convert categorical features to factors
df$MY = as.factor(df$MY)
df$MileGrpATC = as.factor(df$MileGrpATC)

#### filter the data for a single MMT
df = df[df$MMT == 'Acura MDX', ]

##create stan data to feed to the stan hierarchical model
#### observations, their length, and a matrix that counts the number of
#### observations per MMT:MY:MileGrpATC. For incomplete MMT, the matrix is
#### filled with zeros
y = df$q19_2[!is.na(df$q19_2)]
N = length(y)
N_MY_MileGrpATC = df %>%
  group_by(MY, MileGrpATC) %>%
  summarise(count = n())
N_MY_MileGrpATC = matrix(complete(N_MY_MileGrpATC, 
                                  MY, 
                                  MileGrpATC, 
                                  fill = list(count = 0))$count,
                         nrow = years, 
                         ncol = mile_groups, 
                         byrow = TRUE)

#### creating a list of data objects to feed into the stan model
stan_data = list(mode = mean(df$q19_2, na.rm = TRUE),
                 N = N,
                 n_years = years,
                 n_grps = mile_groups,
                 N_MY_MileGrpATC = N_MY_MileGrpATC, 
                 y = y)

## set Stan options for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## write the stan hierarchical model
## NOTE:
#### Beta distribution are parametrized with mode and concentration
#### the model is built for a single MMT but will work dynamically with all MMTs

## model with non-informative prior
#### the following stan model assumes a certain hierarchy where
#### y^t_ijk | theta_ijk ~ Bernoulli(theta_ijk),
#### theta_ijk | alpha_ij ~ Beta(alpha_ij, nu_ij),
#### alpha_ij | rho_i ~ Beta(rho_i, kappa),
#### nu_ij ~ Gamma(0.01, 0.01)
#### rho_i ~ Beta (1, 1)
#### kappa ~ Gamma (0.01, 0.01)
#### i in MMT, j in MMT:MY, k in MMT:MY:MileGrpATC

# write(
#   "//Bernoulli hierarchical model
# data {
#   int <lower=0> N;                                       //number of sample size per MMT
#   int <lower=0> n_years;                                 //number of unique MY
#   int <lower=0> n_grps;                                  //number of unique MileGrpATC
#   int <lower=0> N_MY_MileGrpATC[n_years, n_grps];        //number of observations per MY, MileGrpATC
#   int <lower=0, upper=1> y[N];                           //outcome
# }
#
# parameters {
#   real <lower=2> kappa;                                  //MMT level concentration parameter
#   real <lower=2> nu[n_years];                            //MMT-MY level concentration parameter
#   real <lower=0, upper=1> rho;                           //MMT level mode parameter for Beta priors
#   real <lower=0, upper=1> alpha[n_years];                //MMT-MY level mode parameter for Beta priors
#   real <lower=0, upper=1> theta[n_years, n_grps];        //problem rate for MMT-MY-MileGrpATC level Bernoulli response
# }
#
# model {
#   kappa ~ gamma(0.01, 0.01);                                                                        //MMT-level prior concentration parameter
#   rho ~ beta(1, 1);                                                                                 //MMT-level prior mode parameter
#   int n = 0;                                                                                        //to keep track of number of observation per MMT:MY:MileGrpATC
#   for (j in 1:n_years){                                                                             //loop over MY for MY- and MY:MileGrpATC-level priors                                                             
#     if (sum(N_MY_MileGrpATC[j, ]) > 0){                                                             //assign priors for MY- and MY:MileGrpATC-level parameters only if there is at least one observation in that MY
#       alpha[j] ~ beta(rho * (kappa - 2) + 1, (1 - rho) * (kappa - 2) + 1);                          //MY-level prior mode parameter
#       nu[j] ~ gamma(0.01, 0.01);                                                                    //MY-level prior concentration parameter
#       for (k in 1:n_grps){                                                                          //loop over MY:MileGrpATC-level priors
#         if (N_MY_MileGrpATC[j, k] > 0){                                                             //assign priors for over MY:MileGrpATC-level parameters only if there is at least one observation per MY:MileGrpATC
#           theta[j, k] ~ beta(alpha[j] * (nu[j] - 2) + 1, (1 - alpha[j]) * (nu[j] - 2) + 1);         //MMT:MY:MileGrpATC-level prior problem rate
#           for (t in 1:N_MY_MileGrpATC[j, k]){                                                       //loop over vehicle observations per MMT:MY:MileGrpATC for model specification
#             y[t + n] ~ bernoulli_logit(theta[j, k]);                                                //vehicle observation is a function of problem rate per MMT:MY:MileGrpATC
#           }
#           n = n + N_MY_MileGrpATC[j, k];                                                            //to keep track of number of observation per MMT:MY:MileGrpATC 
#         }
#       }
#     }
#   }
# }
#
# ", 
# "hier_mmt_gen.stan")

## model with informative prior
#### the following stan model assumes a certain hierarchy where
#### y^t_ijk | theta_ijk ~ Bernoulli(theta_ijk),
#### theta_ijk | alpha_ij ~ Beta(alpha_ij, nu_ij),
#### alpha_ij | rho_i ~ Beta(rho_i, kappa),
#### nu_ij ~ Gamma(0.01, 0.01)
#### rho_i ~ Beta (mode, 100)
#### kappa ~ Gamma (0.01, 0.01)
#### i in MMT, j in MMT:MY, k in MMT:MY:MileGrpATC
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
"hier_mmt_BB.stan")

## fitting the Bayesian model with MCMC
fit <- stan(
  file = here::here("hier_mmt_gen.stan"), # stan hierarchical model built above
  data = stan_data, # feeding stan data
  iter = 20000, # number of MCMC iterations
  warmup = 1000, # number of warm-up iterations
  chains = detectCores(), # number of chains
  cores = detectCores(),
  thin = 10, # thinning factor
  init_r = 0, # initialization factor
  #control = list(adapt_delta = 0.99), # learning rate
  seed = 1231
)

## plotting the MCMC chain results for diagnosis
#### divergence plot of MCMC chains, chain mix
fit %>%
  mcmc_trace()

## posterior estimates
#### the functions below are helpers to evaluate posterior estimates

max_density_func <- function(x) {
  #### function to calculate the mode of a distribution
  max_density_x = which.max(density(x)$y)
  return(density(x)$x[max_density_x])
}

hpd_func <- function(x, limit = 'l') {
  ## DEPRECATED
  #### function to evaluate highest posterior density (hpd) limits
  #### given 'limit' it will return lower or upper hpd limits
  hpd_interval = emp.hpd(x, conf = 0.95)
  if (limit == 'l'){
    return(hpd_interval[1])
  } else {
    return(hpd_interval[2])
  }
}

## extract mode of theta; or problem rates
#### for MMTs with less MileGrpATC or MY, the model still estimates theta for 
#### all possible MY and MileGrpATC. Those additional combinations will always
#### have a median of 0.5 and are uniformly distributed
theta_values <- tibble(
  j = sort(rep(1:years, mile_groups)),
  k = rep(1:mile_groups, years),
  .variable = str_c(expression("theta"), "_", k, "_", j),
  values = as.vector(matrix(apply(matrix(extract(fit)$theta, 
                                         ncol = years * mile_groups, 
                                         byrow = FALSE), 
                                  2, 
                                  max_density_func), 
                            ncol = years * mile_groups, 
                            byrow = TRUE))) %>% dplyr::select(.variable, values)

## DEPRECATED
#### extract mode and highest posterior density (hpd) limits
#### only uncomment if interested in knowing the hpd limits
# theta_values <- tibble(
#  j = sort(rep(1:years, mile_groups)),
#  k = rep(1:mile_groups, years),
#  .variable = str_c(expression("theta"), "_", k, "_", j),
#  values = as.vector(matrix(apply(matrix(extract(fit)$theta, 
#                                         ncol = years * mile_groups, 
#                                         byrow = FALSE), 
#                                  2, 
#                                  max_density_func), 
#                            ncol = years * mile_groups, 
#                            byrow = TRUE)),
#  l_hpd = as.vector(matrix(apply(matrix(extract(fit)$theta,
#                                        ncol = years * mile_groups,
#                                        byrow = FALSE),
#                                 2,
#                                 function(x)
#                                 hpd_func(x, limit = 'l')),
#                           ncol = years * mile_groups,
#                           byrow = TRUE)),
#  u_hpd = as.vector(matrix(apply(matrix(extract(fit)$theta,
#                                        ncol = years * mile_groups,
#                                        byrow = FALSE),
#                                 2,
#                                 function(x)
#                                 hpd_func(x, limit = 'u')),
#                           ncol = years * mile_groups,
#                           byrow = TRUE))) %>% dplyr::select(.variable, values, l_hpd, u_hpd)

## draw the distribution of problem rate
#### modes are highlighted with a red line
#### medians are highlighted with a black diamond
fit %>%
  gather_draws(theta[j, k]) %>%
  unite(.variable, .variable, j, k) %>%
  ggplot(aes(x = .value, y = .variable)) +
  geom_halfeyeh(.width = .95) +
  geom_vline(aes(xintercept = values), theta_values, color = "red") +
  geom_text(aes(x = values, label = round(values, 3), y = 1.5), theta_values) + 
  facet_wrap(
    ~ .variable,
    nrow = years,
    ncol = mile_groups,
    scales = "free"
  )

## uncomment below to save working space
# save.image("data-Bayes_Hierarchical_Bernoulli-Beta.RData")