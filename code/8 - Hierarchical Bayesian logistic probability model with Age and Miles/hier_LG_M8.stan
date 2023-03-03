//Bernoulli hierarchical model
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
  //real Beta_0[n_years, n_mmt];       //Intercept for the logistic regression at Make-MMT-MY level
  //real Beta_1[n_years, n_mmt];         //Logistic regression coefficient for age at Make-MMT-MY level
  //real Beta_2[n_years, n_mmt];         //Logistic regression coefficient for miles at Make-MMT-MY level
  real Beta[n_years, n_mmt, 3];
}

model {
  int n = 0;                                                                                //** The int n initialization was moved to the first line. 
  kappa ~ gamma(1, 1);                                                        
  mu ~ normal(0, 1);                                                                                
  for (i in 1:n_mmt){                                                                             
    rho[i] ~ normal(mu, 1 / kappa);                                                                
    for (j in 1:n_years){                                                                                                             
      if (N_MY[j,i] > 0){                                                      
        alpha[i, j] ~ normal(rho[i], 1/ kappa);                                                     
        //Beta_0[j, i] ~ normal(alpha[i, j], 1 / kappa);
        //Beta_1[j, i] ~ normal(alpha[i, j], 1 / kappa);                                      //MMT:MY-level prior problem rate
        //Beta_2[j, i] ~ normal(alpha[i, j], 1 / kappa);                                      //MMT:MY-level prior problem rate
        Beta[j,i,:] ~ normal(alpha[i, j], 1 / kappa);
        
        for (t in 1:N_MY[j, i]){                                                  
          //y[t + n] ~ bernoulli_logit(Beta_0[j, i] + Beta_1[j, i] * age[t + n] + Beta_2[j, i] * miles[t + n]);  
          y[t + n] ~ bernoulli_logit(Beta[j, i, 1] + Beta[j, i, 2] * age[t + n] + Beta[j, i, 3] * miles[t + n]);   
        }
        n = n + N_MY[j, i];                                                       
      }
    }
  }
}


