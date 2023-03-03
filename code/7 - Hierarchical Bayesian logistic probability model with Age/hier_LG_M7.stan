//Bernoulli hierarchical model
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
  real Beta_1[n_years, n_grps, n_mmt];         //Logistic regression coefficient for age at Make-MMT-MY-MileGrpATC level
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
            Beta_1[j, k, i] ~ normal(omega[j, k, i], 1 / kappa);                                      //MMT:MY:MileGrpATC-level prior problem rate
            Beta_0[j, k, i] ~ normal(0, 1 / kappa);
            for (t in 1:N_MY_MileGrpATC[j, k, i]){                                                  //loop over vehicle observations per MMT:MY:MileGrpATC for model specification
              y[t + n] ~ bernoulli_logit(Beta_0[j, k, i] + Beta_1[j, k, i] * age[t + n]);             //vehicle observation is a function of problem rate per MMT:MY:MileGrpATC
            }
            n = n + N_MY_MileGrpATC[j, k, i];                                                       //to keep track of number of observation per MMT:MY:MileGrpATC 
          }
        }
      }
    }
  }
}

generated quantities {
  vector[N] log_lik;
  int n = 0;
  for (i in 1:n_mmt){
    for (j in 1:n_years){
      if (sum(N_MY_MileGrpATC[j, , i]) > 0){
        for (k in 1:n_grps){
          if (N_MY_MileGrpATC[j, k, i] > 0){
            for (t in 1:N_MY_MileGrpATC[j, k, i]){
              log_lik[t + n] = bernoulli_logit_lpmf(y[t + n] | Beta_0[j, k, i] + Beta_1[j, k, i] * age[t + n]);
            }
            n = n + N_MY_MileGrpATC[j, k, i];
          }
        }
      }    
    }
  }
}

