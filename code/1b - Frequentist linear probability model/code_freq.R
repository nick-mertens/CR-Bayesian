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

################################################################################
#### SET INITIAL VARIABLES #####################################################
################################################################################
df = read.csv("mixeddata052722.csv")

## NOTE: 
#### to run for any problem area, please change q19_2 throughout the code

################################################################################
#### FEATURE ENGINEERING #######################################################
################################################################################

## data pre-processing
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

## Adding optional feature
#### Age can be an additional optional feature in the design matrix
# mat$age = df$Age #uncomment to add age as a feature

## save R workspace to speed up process
# save.image("FE.RData") # uncomment to save the working space

################################################################################
#### FREQUENTIST MODEL #########################################################
################################################################################

## load R workspace if not repeating the feature engineering process
# load("FE.RData") #uncomment if loading saved workspace

## create the regression formula
#### problem area ~ BX + e
#### ~. in formula means every feature in the data set
formula = formula(terms(as.formula("q19_2 ~ ."), data = mat))

## fit a linear model
lm_res = lm(q19_2 ~ ., data = mat)

## tidy-up the results
#### replace NA features with 0 for estimate and . for other statistics
#### features that made the desing matrix X "quasi-singluar" are set to NA by lm
tidy_lm = tidy(lm_res)
tidy_lm$estimate[is.na(tidy_lm$estimate)] = 0
tidy_lm$std.error[is.na(tidy_lm$std.error)] = "."
tidy_lm$statistic[is.na(tidy_lm$statistic)] = "."
tidy_lm$p.value[is.na(tidy_lm$p.value)] = "."

## write regression results (parameter estimates) to a csv file
#### csv and workspace save
write.csv(tidy_lm, "q19_2-FREQ_COEF.csv")
# save.image("LM.RData") # uncomment to save the working space

## TO DO:
#### explore biglm & speedglm to speed up process
#### explore sparse matrices to speed up process

################################################################################
#### PREDICT FREQUENTIST CELL PROBABILITIES ####################################
################################################################################

## load R workspace if not repeating the model fitting process
# load("LM.RData") #uncomment if loading saved workspace

## create an empty dataframe to store prediction probabilities
#### predictions are done for each MMT:MY:MileGrpATC
pred_df = data.frame(matrix(ncol = 4, nrow = 0))
colnames(pred_df) = c('MMT', 'MY', 'MileGrpATC', 'prob')

## To speed up process, prediction is looped over MMT
#### filter the design matrix for each MMT and remove duplicates
#### for each unique MMT, MY, MileGrpATC combination
######## extract specific MMT, MY, MileGrpATC and produce a prediction
for (mmt in unique(df$MMT)){
  
  row = which(df$MMT == mmt)
  pred_x = mat[row, !names(mat) %in% c('q19_2')]
  pred_x = setDF(pred_x)
  pred_x = pred_x[!duplicated(pred_x),]
  
  for (i in 1:dim(pred_x)[1]){
    cols = which(pred_x[i, ] == 1)
    col_names = colnames(pred_x)[cols]
    mmt = sub("MMT", "", col_names[1])
    my = sub("MY", "", str_extract(col_names[2], "MY\\d{4}"))
    mpg = sub("MileGrpATC", "", str_extract(col_names[3], "MileGrpATC\\d{1}"))
    pred = predict.lm(lm_res, pred_x[i, ])
    pred_df[nrow(pred_df) + 1, ] = c(mmt, my, mpg, pred)
  }
}

## write prediction results (MileGrpATC leve probabilities) to a csv file
#### csv and workspace save
pred_df = pred_df[!is.na(pred_df$MMT), ]
write.csv(pred_df, "q19_2-FREQ_prob_rate.csv")

## write prediction results (MY level probabilities) to a csv file
#### evaluate cell probabilities for each MY
#### assuming all model years have equal weight
mean_df = pred_df %>% 
  group_by(MMT, MY) %>%
  summarise(MY_prob_rate = mean(as.numeric(prob)))
write.csv(mean_df, "q19_2-FREQ_mean_rate.csv")