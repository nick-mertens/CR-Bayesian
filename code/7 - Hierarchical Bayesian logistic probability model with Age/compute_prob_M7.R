#### Compute Prediction Probability for Model 7

working_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code/7 - Hierarchical Bayesian logistic probability model with Age"
setwd(working_directory)

# Source from main script- should be in the same directory
source("code_Bayes hierarchical LG_functionalized.R")

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
  stan_data = stan_data_func(temp_df, years, mile_groups, problem_area)
  
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
  
  return(res_df)
}

# Example Run
res_df = compute_prob(df, fit_model_name="models/fit_Nissan_M7_20000_12.rds", 
             problem_area = "q19_2", coef_mode="mode")
  