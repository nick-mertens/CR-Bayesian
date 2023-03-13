//Bernoulli hierarchical model
data {
  int <lower=0> N;                             //number of sample size per MMT
  int <lower=0> n_years;                       //number of unique MY
  int <lower=0> n_mmt;                         //number of unique MMT per MakeName
  int <lower=0> N_MY[n_years, n_mmt];          //number of observations per MY
  int <lower=0, upper=1> y[N];                 //outcome
  real <lower=0, upper=1> age[N];              //age
  real <lower=0, upper=1> miles[N];            //Mileage
}

parameters {
  real <lower=0> kappa[2];                      //std for all levels of hierarchy
  real mu[2];                                   //Make-level mean parameters for normal priors
  real rho[2, n_mmt];                           //MMT-level mean parameters for normal priors
  real alpha[2, n_mmt, n_years];                //My-level mean parameters for normal priors
  real Beta_0[n_years, n_mmt];                  //Intercept for the logistic regression at Make-MMT-MY level
  real Beta[n_years, n_mmt, 2];                 //Coefficients for Age and Mileage at Make-MMT-MY level
}

model {
  int n = 0;                                    //** The int n initialization was moved to the first line. 
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

