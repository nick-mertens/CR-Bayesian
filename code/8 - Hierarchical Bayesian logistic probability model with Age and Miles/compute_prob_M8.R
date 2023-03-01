#### Compute Prediction Probability for Model 8

working_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code/8 - Hierarchical Bayesian logistic probability model with Age and Miles"
setwd(working_directory)

# Source from main script- should be in the same directory
source("bayes_hier_LP_M8.R")

# Read original file
df = read.csv("mixeddata052722.csv")

compute_prob <- function(df, fit_model_name, problem_area, coef_mode=c("mode","mean")) {
  ## Compute naive and predicted probabilities by MMT-MY
  # Extract Make Name from fit model
  make = str_split(fit_model_name, "_")[[1]][2]
  
  # Filter data
  temp_df = data_filter_func(df, make, problem_area)
  years = length(unique(temp_df$MY))
  min_year = min(temp_df$MY)
  
  # convert categorical features to factors
  temp_df$MY = as.factor(temp_df$MY)
  # create Stan Data 
  stan_data = stan_data_func(temp_df, years, problem_area)
  
  # Call Stan model
  loaded_fit <- readRDS(fit_model_name)
  
  # Predict
  coef = extract_coef_func(loaded_fit, coef_mode)
  pred_res = predict_func(stan_data, coef, temp_df)
  
  temp_df['y_pred'] = pred_res$y_new
  temp_df = temp_df %>%
    select(MakeName, MMT, MY, problem_area, y_pred)
  
  # Compute naive probability (# of total problem count / # of total row count) and predicted probability
  res_df = temp_df %>%
    group_by(MakeName, MMT, MY) %>% 
    summarise(cnt=n(), round(across(everything(), list(mean=mean)), 4))
  
  res_df = subset(res_df, select = -c(cnt_mean))
  
  print(res_df)
  
  return(res_df)
}

# Example Run
res_df = compute_prob(df, fit_model_name="models/fit_Mercedes-Benz_M8_20000_12.rds", 
             problem_area = "q19_2", coef_mode="mean")
  