library(stringr)

#Setting work Directory
org_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code"
setwd(org_directory)

m7_dir = "7 - Hierarchical Bayesian logistic probability model with Age/"
m8_dir = "8 - Hierarchical Bayesian logistic probability model with Age and Miles/"
m7_stan = "models/fit_Nissan_M7_20000_12.rds"
m8_stan = "models/fit_Nissan_M8_20000_4.rds"
problem_area = "q19_2"

compare_prediction <- function() {
  # Make sure both models pertain to same Makes
  if (str_split(m7_stan, "_")[[1]][2] != str_split(m7_stan, "_")[[1]][2]) {
    stop("The Stanfit objects pertain to different Makes")
  }
  
  # Source from Model 7 
  source(paste(m7_dir, "code_Bayes hierarchical LG_functionalized.R", sep = ""))
  
  M7_res_df = pred_prob(df, fit_model_name = m7_stan, 
                           problem_area = problem_area, coef_mode="mode")
  
  setwd(org_directory)
  
  # Source from Model 8 
  source(paste(m8_dir, "bayes_hier_LP_M8.R", sep = ""))
  
  M8_res_df = pred_prob(df, fit_model_name = m8_stan, 
                           problem_area = problem_area, coef_mode="mode")
  
  merged_df = merge(M7_res_df, M8_res_df, by=c("MakeName","MMT","MY","cnt",paste(problem_area, "_mean", sep="")))
  
  merged_df <- setnames(merged_df, old=c("y_pred_mean.x","y_pred_mean.y"), new=c("pred_mean_M7","pred_mean_M8"))
  
  view(merged_df)
}

compare_prediction()