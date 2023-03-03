#Setting work Directory
# working_directory = "set it to data file destination"
working_directory = "/Users/nick.mertens/Library/CloudStorage/OneDrive-Blend360/Consumer Reports/Documents/2023 Bayesian Modeling - Phase II/CR-Bayesian/code/8 - Hierarchical Bayesian logistic probability model with Age and Miles"

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
library(bayestestR)
library(rstanarm)
library(comprehenr)
library(ggplot2)
library(loo)

# load in the data
df = read.csv("mixeddata052722.csv")

# data preprocessing
years = length(unique(df$MY))
min_year = min(df$MY)
df$MY = as.factor(df$MY)

# filter dataset to one make
df = df[df$MakeName == 'Acura', ]
# df = df[df$MakeName == 'Nissan', ]
n_mmt = length(unique(df$MMT))
df = df[!is.na(df$Age) & !is.na(df$Miles) & !is.na(df$q19_2), ]

# standardize features
df <- cbind(df, Age_scaled = scaler(df$Age, min = 0, max = 1))
df <- cbind(df, Miles_scaled = scaler(df$Miles, min = 0, max = 1))


##create stan data to feed to the stan hierarchical model
#### observations, their length, and a matrix that counts the number of
#### observations per MMT:MY. For incomplete MMT, the matrix is
#### filled with zeros
y = df$q19_2
age = df$Age_scaled
miles = df$Miles_scaled
N = length(y)

N_MY = df %>%
  group_by(MMT, MY) %>%
  summarise(count = n())

# Add dummy column to avoid bug in complete() 
N_MY = cbind(list(1:nrow(N_MY)), N_MY) 

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

stan_data = list(N = N,
                 n_years = years,
                 n_mmt = n_mmt,
                 N_MY = N_MY, 
                 y = y,
                 age = age,
                 miles = miles)

## set Stan options for parallel computing
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## write the stan hierarchical model
write(
  "//Bernoulli hierarchical model
data {
  int <lower=0> N;                                                //number of sample size per MMT
  int <lower=0> n_years;                                          //number of unique MY
  int <lower=0> n_mmt;                                            //number of unique MMT per MakeName
  int <lower=0> N_MY[n_years, n_mmt];                             //number of observations per MY
  int <lower=0, upper=1> y[N];                                    //outcome
  real <lower=0> age[N];                                          //age
  real <lower=0> miles[N];                                        //Mileage
}
parameters {
  real <lower=0> kappa;                      //std for all levels of hierarchy
  real mu;                                   //Make-level mean parameter for normal priors
  real rho[n_mmt];                           //MMT-level mean parameter for normal priors
  real alpha[n_mmt, n_years];                //My-level mean parameter for normal priors
  real Beta_0[n_years, n_mmt];               //Intercept for the logistic regression at Make-MMT-MY level
  //real Beta_1[n_years, n_mmt];             //Logistic regression coefficient for age at Make-MMT-MY level
  //real Beta_2[n_years, n_mmt];             //Logistic regression coefficient for miles at Make-MMT-MY level
  real Beta[n_years, n_mmt, 2];              //Coefficients for Age & Mileage at Make-MMT-MY level - Combined Beta_1 & Beta_2 into Beta array
}

model {
  int n = 0;                                                                                //** The int n initialization was moved to the first line. 
  kappa ~ gamma(1, 1);                                                                      //Make-level prior uncertainty parameter
  mu ~ normal(0, 1);                                                                        //Make-level prior shape parameter        
  for (i in 1:n_mmt){                                                                             
    rho[i] ~ normal(mu, 1 / kappa);                                                         //MMT-level prior shape parameter       
    for (j in 1:n_years){                                                                                                             
      if (N_MY[j,i] > 0){                                                                   //assign priors for MY-level parameters only if at least one obs in MY
        alpha[i, j] ~ normal(rho[i], 1/ kappa);                                             //MY-level prior shape parameter        
        Beta_0[j, i] ~ normal(0, 1 / kappa);
        //Beta_1[j, i] ~ normal(alpha[i, j], 1 / kappa);                                      //MMT:MY-level prior problem rate
        //Beta_2[j, i] ~ normal(alpha[i, j], 1 / kappa);                                      //MMT:MY-level prior problem rate
        Beta[j,i,:] ~ normal(alpha[i, j], 1 / kappa);
        
        for (t in 1:N_MY[j, i]){                                                  
          //y[t + n] ~ bernoulli_logit(Beta_0[j, i] + Beta_1[j, i] * age[t + n] + Beta_2[j, i] * miles[t + n]);
          y[t + n] ~ bernoulli_logit(Beta_0[j, i] + Beta[j, i, 1] * age[t + n] + Beta[j, i, 2] * miles[t + n]);   
        }
        n = n + N_MY[j, i];                                                                 //to keep track of number of observation per MMT:MY
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
"hier_LG_model8.stan")

## fitting the Bayesian model with MCMC
fit <- stan(
  file = here::here("code", "8 - Hierarchical Bayesian logistic probability model with Age and Miles","hier_LG_model8.stan"), # stan hierarchical model built above
  data = stan_data, # feeding stan data
  iter = 8000, # number of MCMC iterations
  warmup = 1000, # number of warm-up iterations
  chains = 8, # number of chains
  thin = 10, # thinning factor
  init_r = 0, # initialization factor
  control = list(adapt_delta = 0.99), # learning rate
  seed = 1231,
  verbose = FALSE
)

# Save the model data
saveRDS(fit, "fit_model8_nissan_standardized.rds")

# If model already run, load it in
fit <- readRDS("fit_model8_acura_standardized.rds")

## checking the loo CV
fit_log <- extract_log_lik(fit, "log_lik", merge_chains = F)
fit_eff <- relative_eff(fit_log)
waic(fit_log)
loo(fit_log, r_eff = fit_eff)

# plotting the MCMC chain results for diagnosis
plot(fit, pars = c("mu", "Beta_0", "Beta_1", "Beta_2"))
print(fit, pars = c("Beta_0", "Beta_1", "Beta_2"))
traceplot(fit, pars = c("Beta_0", "Beta_1", "Beta_2"))

## posterior estimates
max_density_func <- function(x) {
  ## function to calculate the mode of a distribution
  max_density_x = which.max(density(x)$y)
  return(density(x)$x[max_density_x])
}

# Extracting coefficients
## grabbing both mean and mode of distributions
mcmc_df = as.data.frame(fit)
Beta_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta')]
Beta_df = mcmc_df[, names(mcmc_df) %in% Beta_cols]
Beta_mean_values = colMeans(Beta_df) # Mean
Beta_mode_values = apply(Beta_df, MARGIN=2, FUN=max_density_func) # Mode

Beta0_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_0[')]
Beta0_df = mcmc_df[, names(mcmc_df) %in% Beta0_cols]
Beta0_mean_values = colMeans(Beta0_df)
Beta0_mode_values = apply(Beta0_df, MARGIN=2, FUN=max_density_func)

#Beta1_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_1[')]
Beta1_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta[') & endsWith(colnames(mcmc_df), '1]')]
Beta1_df = mcmc_df[, names(mcmc_df) %in% Beta1_cols]
Beta1_mean_values = colMeans(Beta1_df)
Beta1_mode_values = apply(Beta1_df, MARGIN=2, FUN=max_density_func)

#Beta2_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta_2[')]
Beta2_cols = colnames(mcmc_df)[startsWith(colnames(mcmc_df), 'Beta[') & endsWith(colnames(mcmc_df), '2]')]
Beta2_df = mcmc_df[, names(mcmc_df) %in% Beta2_cols]
Beta2_mean_values = colMeans(Beta2_df)
Beta2_mode_values = apply(Beta2_df, MARGIN=2, FUN=max_density_func)

y_new = c()
beta_df = data.frame(matrix(ncol = 5, nrow = 0))
colnames(beta_df) = c('MMT', 'MY', 'b0', 'b1', 'b2')

n = 0
for (i in 1:stan_data$n_mmt){
  for (j in 1:stan_data$n_years){
    if (stan_data$N_MY[j, i] > 0){
      for (t in 1:stan_data$N_MY[j, i]){
        b_name  = paste('Beta', j, ',', i, ']', sep = "")
        b0_name = paste('Beta_0[', j, ',', i, ']', sep = "")
        #b1_name = paste('Beta_1[', j, ',', i, ']', sep = "")
        #b2_name = paste('Beta_2[', j, ',', i, ']', sep = "")
        b1_name = paste('Beta[', j, ',', i, ',1]', sep = "")
        b2_name = paste('Beta[', j, ',', i, ',2]', sep = "")
        Beta = Beta_mode_values[b_name] # change mode to mean depending on what values you want to use
        Beta0 = Beta0_mode_values[b0_name]
        Beta1 = Beta1_mode_values[b1_name]
        Beta2 = Beta2_mode_values[b2_name]
        predict = 1 / (1 + exp(-(Beta0 + Beta1 * stan_data$age[t + n] + Beta2 * stan_data$miles[t + n]))) # should Beta2 be added or subtracted?
        y_new = c(y_new, predict)
      }
      mmt = unique(df$MMT)[i]
      my = unique(df[df$MMT == mmt, ]$MY)[j]
      beta_df[nrow(beta_df) + 1, ] = c(mmt, my, Beta0, Beta1, Beta2)
      n = n + stan_data$N_MY[j, i]
      
    }
  }
}

df['y_pred'] = y_new
df <- df %>%
  select(MakeName, MMT, MY, Age, Miles, q19_2, y_pred)

print(df %>% group_by(MMT, MY) %>% summarise(mean_pred_prob=mean(y_pred)))
print(df %>% group_by(MMT, MY) %>% summarise(prob=sum(q19_2) / n()))