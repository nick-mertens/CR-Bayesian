#Setting work Directory
working_directory = "set it to data file destination"
setwd(working_directory)

#Importing required packages
library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(data.table)
library(coda)
library(lmtest)
library(sandwich) 
library(tidyverse)
library(lme4)
library(broom)
library(stringr)
library(brms)
library(cmdstanr)
library(parallel)

################################################################################
#### SET INITIAL VARIABLES #####################################################
################################################################################
df = read.csv("mixeddata052722.csv")

## NOTE: 
#### to run for any problem area, please change q19_2 throughout the code
#### to run Bayesian logistic or linear regression, uncomment line 315, 316

################################################################################
#### FEATURE ENGINEERING #######################################################
################################################################################

## data pre-processing
#### Reading data
df = read.csv("mixeddata052722.csv")

#### convert categorical features to factors
df$MY = as.factor(df$MY)
df$MileGrpATC = as.factor(df$MileGrpATC)
df$MMT = as.factor(df$MMT)

#### create 0-1 nested features
#### model.matrix creates numerical matrix; for factors to binary
#### -1 removes the intercept in binary creation
mat1 = model.matrix(q19_2 ~ MMT - 1, data = df)
mat2 = model.matrix(q19_2 ~ MMT:MY - 1,data = df)
mat3 = model.matrix(q19_2 ~ MMT:MY:MileGrpATC - 1, data = df)
mat = cbind.data.frame(mat1, 
                       mat2, 
                       mat3, 
                       q19_2 = df$q19_2) # put features and response together
rm(list = c('mat1', 'mat2', 'mat3')) # remove unnecessary objects


## NOTE:
#### the nested feature engineering creates impossible combinations
#### impossible combinations are those without any observation in the data
#### the process below will identify those features and remove them
#### in doing so, we also change the naming convention to match CR's SAS output

## removing impossible combinations
#### determine combinations of MMT, MY, and MileGrpATC supported by the data
agg_MMT = aggregate(df$MMT, by = list(df$MMT), FUN = length)
agg_MMT_MY = aggregate(df$MMT, by = list(df$MMT, df$MY), FUN = length)
agg_MMT_MY_MileGrpATC = aggregate(df$MMT, 
                                  by = list(df$MMT, df$MY, df$MileGrpATC), 
                                  FUN = length)

#### feature name formatting
#### changing R feature names to match with that of CR's SAS
agg_MMT_names = paste("MMT", agg_MMT$Group.1, sep = "")
MMT_MY_temp1 = paste("MMT", agg_MMT_MY$Group.1, sep = "")
MMT_MY_temp2 = paste("MY", agg_MMT_MY$Group.2, sep = "")
agg_MMT_MY_names = paste(MMT_MY_temp1, MMT_MY_temp2, sep = ":")
MMT_MY_MileGrpATC_temp1 = paste("MMT", agg_MMT_MY_MileGrpATC$Group.1, sep ="")
MMT_MY_MileGrpATC_temp2 = paste("MY", agg_MMT_MY_MileGrpATC$Group.2, sep ="")
MMT_MY_MileGrpATC_temp3 = paste("MileGrpATC", agg_MMT_MY_MileGrpATC$Group.3, 
                                sep = "")
agg_MMT_MY_MileGrpATC_names = paste(MMT_MY_MileGrpATC_temp1, 
                                    MMT_MY_MileGrpATC_temp2, 
                                    MMT_MY_MileGrpATC_temp3, sep = ":")

#### merging all transformed feature names into a single list
positive_cnt_combinations = append(agg_MMT, agg_MMT_MY_combined)
positive_cnt_combinations = append(positive_cnt_combinations, 
                                   agg_MMT_MY_MileGrpATC_combined)
rm(list = c('agg_MMT', 'MMT_MY_temp1', 'MMT_MY_temp2',
            'agg_MMT_MY_combined', 'MMT_MY_MileGrpATC_temp1', 
            'MMT_MY_MileGrpATC_temp2', 'MMT_MY_MileGrpATC_temp3',
            'agg_MMT_MY_MileGrpATC_combined',
            'agg_MMT', 'agg_MMT_MY', 
            'agg_MMT_MY_MileGrpATC')) # remove unnecessary objects

#### extract all the combinations of MMT, MY, MileGrpATC; 
#### those with or without available data
all_combinations = c(colnames(mat))

#### identifying set difference of feature names between all potential ...
#### ...combinations and those supported by the data
#### removing impossible nesting combinations of MMT, MY, MileGrpATC;
#### those without any data available
set_diff = setdiff(all_combinations, positive_cnt_combinations)
set_diff = set_diff[1:length(set_diff) - 1] #-1 no to remove problem area term 
mat = mat[, !colnames(mat) %in% set_diff]
rm(list = c('all_combinations', 'positive_cnt_combinations', 'set_diff'))

## save R workspace to speed up process
# save.image("data-FE.RData") # uncomment to save the working space

################################################################################
#### BAYESIAN GLM MODEL: FUNCTIONALIZED ########################################
################################################################################

## load R workspace if not repeating the frequentist process
load("data-FE.RData") #uncomment if loading saved workspace

################################################################################
######## functionalized Bayes regression training for sequential runs of all MMT

data_filter_func = function(df, mat, mmt){
  ## filter data for a single MMT
  #### feature column that correspond to that MMT
  #### indices that correspond to that MMT
  ind = which(df$MMT == mmt)
  col = which(sapply(strsplit(colnames(mat), ":", fixed = TRUE),
                     '[[',
                     1) == paste('MMT', mmt, sep = ""))
  mat = mat[ind, c(col, length(colnames(mat)))]
  
  ## remove collinearity
  #### remove perfect collinearity of 0-1 data
  corr_cols = c()
  mys = as.character(unique(df$MY[df$MMT == mmt]))
  corr_cols = c(corr_cols, last(mys)) #remove latest model year
  for (my in mys){ #remove last mileage group in each model year
    grps = unique(df$MileGrpATC[df$MMT == mmt & df$MY == my])
    corr_cols = c(corr_cols, 
                  paste(my, ":MileGrpATC", last(grps), sep = ""))
  }
  corr_cols = paste(corr_cols, "$", sep = "")
  regex = paste(corr_cols, collapse = "|")
  rm_col = which(grepl(regex, colnames(mat), fixed = FALSE))
  mat = mat[, -c(rm_col)]
  
  return(mat)
}

col_format_func = function(mat_){
  ## replace special chars in column names because of BRMS package
  colnames(mat_) = gsub(":", "_", colnames(mat_))
  colnames(mat_) = gsub(" ", "_", colnames(mat_))
  colnames(mat_) = gsub("/", ".", colnames(mat_))
  colnames(mat_) = gsub("-", ".", colnames(mat_))
  colnames(mat_) = gsub("&", ".", colnames(mat_))
  
  return(mat_)
}

prior_prob_rate_func = function(df, mmt, y_name){
  ## evaluate average problem rate per MMT:MY
  avg_prob_rate = df[df$MMT == mmt, ] %>% 
    group_by(MMT, MY, MileGrpATC) %>%
    summarise(avg = mean(get(y_name), na.rm = TRUE),
              sd = sd(get(y_name), na.rm = TRUE))
  
  avg_mean_rate = df[df$MMT == mmt, ] %>%
    group_by(MMT, MY) %>%
    summarise(avg = mean(get(y_name), na.rm = TRUE),
              sd = sd(get(y_name), na.rm = TRUE))
  
  return(list(prob_rate = avg_prob_rate, mean_rate = avg_mean_rate))
}

set_prior_func = function(mat_, y_name, prior_prob_rate, prior_mean_rate){
  ## set prior for regression coefficient
  #### prior object is an rstan dataframe for each class of parameters
  #### mileage group level and model year level covariates have the same prior:
  ######## this is set to average problem rate of that model year
  ######## model year is extracted from covariate names for appropriate prior
  #### model make twin level covariates have a different prior:
  ######## this is set to average problem rate of that mmt
  prior_lst = data.frame()
  mmt_my_cols = which(grepl("MY2018|MY2019|MY2020", 
                            colnames(mat_), 
                            fixed = FALSE))
  mmt_my_cols = colnames(mat_)[mmt_my_cols]
  mmt_my_MileGrpATC_cols = which(grepl("MileGrpATC1|MileGrpATC2|MileGrpATC3",
                                       colnames(temp_mat),
                                       fixed = FALSE))
  mmt_my_MileGrpATC_cols = colnames(temp_mat)[mmt_my_MileGrpATC_cols]
  
  overall_rate = mean(mat_[[y_name]], na.rm = TRUE)
  # overall_sd = sd(mat_[[y_name]], na.rm = TRUE) #do not uncomment
  for (col in colnames(mat_)) {
    if (col %in% mmt_my_cols){
      if (col %in% mmt_my_MileGrpATC_cols){
        my = str_extract(col, "(?<=MY)\\d{4}")[1]
        grp = str_extract(col, "(?<=MileGrpATC)\\d{1}")[1]
        my_grp_avg = prior_prob_rate$avg[prior_prob_rate$MY == my & 
                                         prior_prob_rate$MileGrpATC == grp]
        if (my_grp_avg != 0){
          dist = paste("normal(", my_grp_avg, ",", 1, ")", sep = "")
        } else {
          dist = paste("normal(", 0, ",", 1, ")", sep = "")
        }
        prior = set_prior(dist, class = "b", coef = col)
        prior_lst = rbind(prior_lst, prior)
      } else {
        my = str_extract(col, "(?<=MY)\\d{4}")[1]
        my_avg = prior_mean_rate$avg[prior_mean_rate$MY == my]
        if (my_avg != 0){
          dist = paste("normal(", my_avg, ",", 1, ")", sep = "")
        } else {
          dist = paste("normal(", 0, ",", 1, ")", sep = "")
        }
        prior = set_prior(dist, class = "b", coef = col)
        prior_lst = rbind(prior_lst, prior)
      }
    } else if (col == y_name) {
      next
    } else {
      dist = paste("normal(", overall_rate, ",", 1, ")", sep = "")
      prior = set_prior(dist, class = "b", coef = col)
      prior_lst = rbind(prior_lst, prior)
    }
  }
  return(prior_lst)
  # To DO:
  # 1- find out why informative priors for sd don't produce results as expected
  # 2- find a better prior for prob_rate = 0 - DONE
}

bayes_fit_func = function(mat_, y_name, prior_lst = NULL, col.names = FALSE) {
  ## main function fitting the Bayesian model
  bayes_res = brm(formula = q19_2 ~ . - 1, # -1 means without intercept
                  data = mat_,
                  prior = prior_lst,
                  warmup = 1000, # warm up iterations
                  iter = 20000, # MCMC iterations
                  chains = detectCores(), # number of MCMC chains
                  init= "0", 
                  thin = 10, # thining factor, using ever other k samples
                  cores = detectCores(), # number of CPUs
                  seed = 1231,
                  backend = 'cmdstanr',
                  #control = list(adapt_delta = 0.99), # learning rate
                  silent = 2,
                  refresh = 0)
  return(bayes_res)
  # TO DO:
  # 1- automatically configure 'cmdstan' back-end engine
  # 2- fix the issue for GPU parallelization with Open Cl
}

bayes_prediction_func = function(bayes_res_obj, mat, temp_mat, mmt, y_name, col.names = FALSE){
  ## cell probabilities
  #### temporary save Bayesian predictions
  mmt_name = paste('MMT', mmt, sep = "")
  col = which(sapply(strsplit(colnames(mat), ":", fixed = TRUE),
                     '[[',
                     1) == paste('MMT', mmt, sep = ""))
  ind = which(df$MMT == mmt)
  col_data = mat[ind, col]
  col_data = col_data[, !names(col_data) %in% y_name]
  col_data = col_data[!duplicated(col_data), ]
  
  new_data = temp_mat
  new_data = new_data[, !names(new_data) %in% y_name]
  new_data = new_data[!duplicated(new_data), ]
  predict_bayes = brms:::predict.brmsfit(bayes_res, newdata = new_data)
  
  pred_df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(pred_df) = c('MMT', 'MY', 'MileGrpATC')
  for (i in 1:dim(col_data)[1]){
    cols = which(col_data[i, ] == 1)
    col_names = colnames(col_data)[cols]
    mmt = sub("MMT", "", col_names[1])
    my = sub("MY", "", str_extract(col_names[2], "MY\\d{4}"))
    mpg = sub("MileGrpATC", "", str_extract(col_names[3], "MileGrpATC\\d{1}"))
    pred_df[nrow(pred_df) + 1, ] = c(mmt, my, mpg)
  }
  pred_df = cbind.data.frame(pred_df, prob = predict_bayes)
  write.table(pred_df, 
              file = "q19_2-BAYES_temp_prob.csv", 
              sep = ",",
              col.names = col.names,
              row.names = FALSE,
              append = TRUE)
  return(pred_df)
}

coef_proccessing_func = function(bayes_res){
  ## DEPRECATED
  ## not needed if regressions are run without an intercept
  #### needed to sum intercept and MMT coefficients when both are in the model
  intercept_df = bayes_res %>% 
    filter(startsWith(row.names(bayes_res), "Intercept")) %>%
    select(Estimate, Est.Error)
  rownames(intercept_df) = gsub("Intercept_", "", rownames(intercept_df))
  
  intercept_filtered = bayes_res %>%
    filter(!startsWith(row.names(bayes_res), "Intercept")) %>%
    select(Estimate, Est.Error)
  
  for (row in rownames(intercept_df)){
    new_est = intercept_df$Estimate[rownames(intercept_df) == row] + 
      intercept_filtered$Estimate[rownames(intercept_filtered) == row]
    new_err = sqrt(intercept_df$Est.Error[rownames(intercept_df) == row] ^ 2 + 
      intercept_filtered$Est.Error[rownames(intercept_filtered) == row] ^ 2)
    intercept_filtered$Estimate[rownames(intercept_filtered) == row] = new_est
    intercept_filtered$Est.Error[rownames(intercept_filtered) == row] = new_err
  }
  return(intercept_filtered)
}

term_name_format_func = function(bayes_res, mat, mmt, y_name, col.names = col.names){
  ## coefficient extraction
  #### temporary save Bayesian coefficients
  bayes_coef = fixef(bayes_res)
  bayes_coef = setDT(as.data.frame(bayes_coef), keep.rownames = TRUE)
  
  ## name formatting function
  #### matching CR'snaming convention
  mmt_name = paste('MMT', mmt, sep = "")
  col = which(sapply(strsplit(colnames(mat), ":", fixed = TRUE),
                     '[[',
                     1) == paste('MMT', mmt, sep = ""))
  
  org_name = colnames(mat[, col])
  org_name_p_rm = gsub("[[:punct:]]", " ", org_name)
  coef_name = bayes_coef$rn
  coef_name_p_rm = gsub("[[:punct:]]", " ", coef_name)
  
  for (term in org_name_p_rm){
    if (term %in% coef_name_p_rm){
      bayes_coef$rn[which(coef_name_p_rm == term)] = org_name[org_name_p_rm == term]
    } else if (term != y_name) {
      row = data.frame(cbind(rn = org_name[org_name_p_rm == term],
                             Estimate = 0, 
                             Est.Error = '.', 
                             Q2.5 = '.',
                             Q97.5 = '.'))
      bayes_coef = rbind.data.frame(bayes_coef, row)
    }
  }
  bayes_coef = bayes_coef[order(bayes_coef$rn)]
  
  write.table(bayes_coef, 
              file = "q19_2-BAYES_temp_COEF.csv", 
              sep = ",",
              col.names = col.names,
              row.names = FALSE,
              append = TRUE)
  
  return(bayes_coef)
}

################################################################################
######## main loop starts here

## generate an empty data frame to store Bayesian results
bayes_coef_df = data.frame()
bayes_pred_df = data.frame()

## initialize an iterator to show progression
i = 0

## Bayesian linear (and logistic) models are run block by block for each MMT
#### 1- feature matrix is filtered for each MMT
#### 2- names are formatted because BRMS package doesn't recognize CR format
#### 3- average problem rate per model year is calculated for ech MMT
#### 4- priors are set for MMT, MMT:MY, and MMT:MY:MileGrpATC features
#### 5- model is run and parameter results + predictions are appended

for (mmt in (unique(df$MMT))){
  col.names = ifelse(i == 0, TRUE, FALSE) #to write col names in the final csv
  i = i + 1
  print(i)
  temp_mat = data_filter_func(df, mat, mmt)
  temp_mat = col_format_func(temp_mat)
  rate_res = prior_prob_rate_func(df, 
                                  mmt, 
                                  'q19_2')
  prior_lst = set_prior_func(temp_mat, 
                             'q19_2', 
                             rate_res$prob_rate, 
                             rate_res$mean_rate)
  bayes_res = bayes_fit_func(temp_mat, 
                             'q19_2', 
                             prior_lst = prior_lst,
                             col.names)
  bayes_pred = bayes_prediction_func(bayes_res,
                                     mat,
                                     temp_mat,
                                     mmt,
                                     'q19_2',
                                     col.names = col.names)
  bayes_coef = term_name_format_func(bayes_res, 
                                     mat,
                                     mmt,
                                     'q19_2',
                                     col.names = col.names)
  bayes_coef_df = rbind(bayes_coef_df, 
                        bayes_coef)
  bayes_pred_df = rbind(bayes_pred_df,
                        bayes_pred)
  gc() #garbage removal to free memory
}

## save R workspace to speed up process
# save.image("data-BAYES_LM.RData") # uncomment to save the working space

################################################################################
######## print results

## DEPRECATED
#### add intercept and MMT estimates
# bayes_coef = coef_proccessing_func(bayes_coef)

## write coefficient results
write.csv(bayes_coef_df, 
          "q19_2-BAYES_LM_COEF.csv")

## write prediction results (MileGrpATC level probabilities) to a csv file
#### csv and workspace save
write.csv(bayes_pred_df,
          "q19_2-BAYES_LM_prob_rate.csv")

## write prediction results (MY level probabilities) to a csv file
#### evaluate cell probabilities for each MY
#### assuming all model years have equal weight
mean_df = bayes_pred_df %>% 
  group_by(MMT, MY) %>%
  summarise(MY_prob_rate = mean(as.numeric(prob.Estimate)))
write.csv(mean_df, 
          "q19_2-BAYES_LM_mean_rate.csv")