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
#df = df[df$MMT == 'Ford Focus', ]
#### filter the data for a single MakeName
df = df[df$MakeName == 'Acura', ]
n_mmt = length(unique(df$MMT))

##create stan data to feed to the stan hierarchical model
#### observations, their length, and a matrix that counts the number of
#### observations per MMT:MY:MileGrpATC. For incomplete MMT, the matrix is
#### filled with zeros
df = df[!is.na(df$Age) & !is.na(df$q19_2), ]
y = df$q19_2
age = df$Age
N = length(y)

N_MY_MileGrpATC = df %>%
  group_by(MMT, MY, MileGrpATC) %>%
  summarise(count = n())

# ADDED THIS CODE TO AVOID ERROR IN lapply BELOW - NM 02/22/23
###############################################################################
# N_MY_MileGrpATC <- N_MY_MileGrpATC %>%
#   mutate(row_num = row_number()) %>%
#   select(row_num, everything())

N_MY_MileGrpATC <- cbind(list(1:nrow(N_MY_MileGrpATC)), N_MY_MileGrpATC)
###############################################################################

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
dim(N_MY_MileGrpATC) = c(3, 3, 3)
#N_MY_MileGrpATC = list(N_MY_MileGrpATC = N_MY_MileGrpATC)

#### creating a list of data objects to feed into the stan model
stan_data = list(N = N,
                 n_years = years,
                 n_grps = mile_groups,
                 n_mmt = n_mmt,
                 N_MY_MileGrpATC = N_MY_MileGrpATC, 
                 y = y,
                 age = age)

## set Stan options for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## write the stan hierarchical model
#### the following stan model assumes a certain hierarchy where
#### y^t_ijk | theta_ijk ~ Bernoulli(theta_ijk),
#### theta_ijk | omega_ijk ~ Beta(omega_ijk, kappa),
#### omega_ijk | alpha_ij ~ Beta(alpha_ij, kappa),
#### alpha_ij | rho_i ~ Beta(rho_i, kappa),
#### rho_i ~ Beta (1, 1)
#### kappa ~ Gamma (1, 1)
#### i in MMT, j in MMT:MY, k in MMT:MY:MileGrpATC

## NOTE:
#### Beta distribution are parametrized with mode and concentration
#### the model is built for a single MMT but will work dynamically with all MMTs
##### CHANGED THIS CODE TO PUT int n = 0; AS FIRST LINE OF model ELEMENT - NM 02/22/23
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
  real alpha[n_mmt, n_years];                //MY-level mean parameter for normal priors
  real omega[n_mmt, n_years, n_grps];        //MileGrpATC-level mean parameter for normal prior
  real Beta[n_mmt, n_years, n_grps];         //Logistic regression coefficient for age at Make-MMT-MY-MileGrpATC level
  real Beta_0[n_mmt, n_years, n_grps];       //Intercept for the logistic regression at Make-MMT-MY-MileGrpATC level
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
"hier_reg_gen_quant.stan")

## fitting the Bayesian model with MCMC
fit <- stan(
  file = here::here("hier_reg_gen_quant.stan"), # stan hierarchical model built above
  data = stan_data, # feeding stan data
  iter = 8000, # number of MCMC iterations
  warmup = 1000, # number of warm-up iterations
  chains = 8, # number of chains
  thin = 10, # thinning factor
  init_r = 0, # initialization factor
  control = list(adapt_delta = 0.99), # learning rate
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


mcmc_df = as.data.frame(fit)
Beta_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta[')]
Beta_df = mcmc_df[, names(mcmc_df) %in% Beta_cols]
Beta_values = colMeans(Beta_df)

Beta0_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_0[')]
Beta0_df = mcmc_df[, names(mcmc_df) %in% Beta0_cols]
Beta0_values = colMeans(Beta0_df)

y_new = c()

for (i in 1:n_mmt){
  for (j in 1:years){
    for (k in 1:mile_groups){
      n = 0
      for (t in 1:N_MY_MileGrpATC[j, k, i]){
        Beta_name = paste('Beta[', j, ',', k, ',', i, ']', sep = "")
        Beta0_name = paste('Beta_0[', j, ',', k, ',', i, ']', sep = "")
        Beta = Beta_values[Beta_name]
        Beta0 = Beta0_values[Beta0_name]
        predict = 1 / (1 + exp(-(Beta0 + Beta * age[t + n])))
        y_new = c(y_new, predict)
      }
      n = n + N_MY_MileGrpATC[j, k, i]
    }
  }
}

df['y_new'] = y_new
prob_df = df %>% 
  group_by(MMT, MY, MileGrpATC) %>%
  summarise(MY_MileGrpATC_prob_rate = mean(y_new))

mean_df = prob_rate_df %>% 
  group_by(MMT, MY) %>%
  summarise(MY_prob_rate = mean(MY_MileGrpATC_prob_rate))

model_type = 'hier_log'
write.csv(prob_df, 
          paste("q19_2-BAYES_",
                model_type,
                "_prob_rate.csv",
                sep = ""))

write.csv(mean_df, 
          paste("q19_2-BAYES_",
                model_type,
                "_mean_rate.csv",
                sep = ""))

ggplot(aes(x = value, y = variable), data = reshape2::melt(Beta_df)) + # ADJUSTED CODE TO AVOID WARNING MESSAGE - NM 02/22/23
  stat_halfeye() + 
  facet_wrap(
    ~ variable,
    nrow = 9,
    ncol = 3,
    scales = "free"
  )

## extract mode of theta; or problem rates
#### for MMTs with less MileGrpATC or MY, the model still estimates theta for 
#### all possible MY and MileGrpATC. Those additional combinations will always
#### have a median of 0.5 and are uniformly distributed
Beta_values <- tibble(
  j = sort(rep(1:years, mile_groups * n_mmt)),
  k = rep(1:mile_groups, years * n_mmt),
  i = rep(1:mile_groups, years * mile_groups),
  # .variable = str_c(expression("Beta"), "_", k, "_", j, "_", i), # CHANGED BY NM 02/22/23
  .variable = paste0(expression("Beta"), "_", k, "_", j, "_", i),
  values = as.vector(matrix(apply(matrix(extract(fit)$Beta, 
                                         ncol = years * mile_groups * n_mmt, 
                                         byrow = FALSE), 
                                  2, 
                                  max_density_func), 
                            ncol = years * mile_groups * n_mmt, 
                            byrow = TRUE))) %>% dplyr::select(.variable, values)


fit %>%
  gather_draws(Beta[j, k, i]) %>%
  unite(.variable, .variable, .variable, j, k, i) %>%
  ggplot(aes(x = .value, y = .variable)) +
  # geom_halfeyeh(.width = .95) + # CHANGED BY NM 02/22/23 - geom_halfeyeh IS DEPRICATED
  stat_halfeye(.width = .95) +
  #geom_vline(aes(xintercept = values), Beta_values, color = "red") +
  #geom_text(aes(x = values, label = round(values, 3), y = 1.5), Beta_values) + 
  facet_wrap(
    ~ .variable,
    nrow = 9,
    ncol = 3,
    scales = "free"
  )

## uncomment below to save working space
# save.image("data-Bayes_Hierarchical_Bernoulli-Beta.RData")