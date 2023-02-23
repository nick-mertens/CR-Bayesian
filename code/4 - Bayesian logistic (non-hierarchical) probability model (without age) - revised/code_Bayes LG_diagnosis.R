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
#### to run a Bayesian Logistic regression, uncomment line 194
#### to diagnose a problem for a different MMT, change line 117

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
# save.image("FE.RData") # uncomment to save the working space

################################################################################
#### BAYESIAN GLM MODEL: SINGLE MMT TRAINING ###################################
################################################################################

## load R workspace if not repeating the feature engineering process
load("data-FE.RData") #uncomment if loading saved workspace

## filter data for a single MMT
#### feature column that correspond to that MMT
#### indecis that correspond to that MMT
MMT = 'Acura MDX'
y_name = 'q19_2'
ind = which(df$MMT == MMT)
col = which(colnames(mat) %like% paste('MMT', MMT, sep = ""))
temp_mat = mat[ind, c(col, length(colnames(mat)))]

## remove collinearity
#### remove perfect collinearity of 0-1 data
corr_cols = c()
mys = as.character(unique(df$MY[df$MMT == MMT]))
corr_cols = c(corr_cols, last(mys)) #remove latest model year
for (my in mys){ #remove last mileage group in each model year
  grps = unique(df$MileGrpATC[df$MMT == MMT & df$MY == my])
  corr_cols = c(corr_cols, 
                paste(my, ":MileGrpATC", last(grps), sep = ""))
}
corr_cols = paste(corr_cols, "$", sep = "")
regex = paste(corr_cols, collapse = "|")
rm_col = which(grepl(regex, colnames(temp_mat), fixed = FALSE))
temp_mat = temp_mat[, -c(rm_col)]

## replace special chars in column names because of BRMS package
colnames(temp_mat) = gsub(":", "_", colnames(temp_mat))
colnames(temp_mat) = gsub(" ", "_", colnames(temp_mat))
colnames(temp_mat) = gsub("/", ".", colnames(temp_mat))
colnames(temp_mat) = gsub("-", ".", colnames(temp_mat))
colnames(temp_mat) = gsub("&", ".", colnames(temp_mat))

## evaluate average problem rate per MMT:MY
avg_prob_rate = df[df$MMT == MMT, ] %>% 
  group_by(MMT, MY, MileGrpATC) %>%
  summarise(avg = mean(get(y_name), na.rm = TRUE),
            sd = sd(get(y_name), na.rm = TRUE))

avg_mean_rate = df[df$MMT == MMT, ] %>%
  group_by(MMT, MY) %>%
  summarise(avg = mean(get(y_name), na.rm = TRUE),
            sd = sd(get(y_name), na.rm = TRUE))

## set prior for regression coefficient
#### prior object is an rstan dataframe for each class of parameters
#### mileage group level and model year level covariates have the same prior:
######## this is set to average problem rate of that model year
######## model year is extracted from covariate names for appropriate prior
#### model make twin level covariates have a different prior:
######## this is set to average problem rate of that mmt
prior_lst = data.frame()
mmt_my_cols = which(grepl("MY2018|MY2019|MY2020", 
                          colnames(temp_mat), 
                          fixed = FALSE))
mmt_my_cols = colnames(temp_mat)[mmt_my_cols]
mmt_my_MileGrpATC_cols = which(grepl("MileGrpATC1|MileGrpATC2|MileGrpATC3",
                                     colnames(temp_mat),
                                     fixed = FALSE))
mmt_my_MileGrpATC_cols = colnames(temp_mat)[mmt_my_MileGrpATC_cols]

## experimenting with prior specification when actual rate = 0
overall_loc = mean(temp_mat[[y_name]], na.rm = TRUE)

for (col in colnames(temp_mat)) {
  if (col %in% mmt_my_cols){
    if (col %in% mmt_my_MileGrpATC_cols){
      my = str_extract(col, "(?<=MY)\\d{4}")[1]
      grp = str_extract(col, "(?<=MileGrpATC)\\d{1}")[1]
      my_grp_avg = avg_prob_rate$avg[avg_prob_rate$MY == my & 
                                     avg_prob_rate$MileGrpATC == grp]
      if (my_grp_avg != 0){
        loc = log(my_grp_avg) - log (1 - my_grp_avg)
        dist = paste("normal(", loc, ",", 10, ")", sep = "")
      } else {
        loc = -1e+13
        dist = paste("normal(", loc, ",", 2, ")", sep = "")
      }
      prior = set_prior(dist, class = "b", coef = col)
      prior_lst = rbind(prior_lst, prior)
    } else {
      my = str_extract(col, "(?<=MY)\\d{4}")[1]
      my_avg = avg_mean_rate$avg[avg_mean_rate$MY == my]
      if (my_avg != 0){
        loc = log(my_avg) - log (1 - my_avg)
        dist = paste("normal(", loc, ",", 10, ")", sep = "")
      } else {
        loc = -1e+13
        dist = paste("normal(", loc, ",", 2, ")", sep = "")
      }
      prior = set_prior(dist, class = "b", coef = col)
      prior_lst = rbind(prior_lst, prior)
    }
  } else if (col == y_name) {
    next
  } else {
    if (overall_loc != 0){
      overall_loc = log(overall_loc) - log (1 - overall_loc)
      dist = paste("normal(", overall_loc, ",", 10, ")", sep = "")
    } else {
      overall_loc = -1e+13
      dist = paste("normal(", overall_loc, ",", 2, ")", sep = "")
    }
    prior = set_prior(dist, class = "b", coef = col)
    prior_lst = rbind(prior_lst, prior)
  }
}

## define and fit a Bayesian regression
#### family option defaults to Gaussian -> linear regression
#### for logistic regression, family = Bernoulli(link = "logit")
bayes_res = brm(formula = q19_2 ~ . - 1,
                data = temp_mat,
                family = bernoulli(link = "logit"), #uncomment for logistic
                prior = prior_lst,
                warmup = 1000,
                iter = 20000,
                chains = detectCores(), 
                init= "0",
                thin = 10,
                cores = detectCores(),
                seed = 1231,
                backend = 'cmdstanr',
                #control = list(adapt_delta = 0.99), # learning rate
                silent = 2,
                refresh = 0,
                #opencl = opencl(c(0, 0))
)

new_data = temp_mat
new_data = new_data[, !names(new_data) %in% y_name]
new_data = new_data[!duplicated(new_data), ]

predict_bayes = brms:::predict.brmsfit(bayes_res, newdata = new_data)
write.table(predict_bayes, 
            file = paste("q19_2-BAYES_diagnosis_",
                         MMT,
                         '_prob_rate.csv',
                         sep = ""),
            sep = ",",
            col.names = TRUE,
            row.names = FALSE,
            append = FALSE)


bayes_coef = fixef(bayes_res)
bayes_coef = setDT(as.data.frame(bayes_coef), keep.rownames = TRUE)

write.table(bayes_coef, 
            file = paste("q19_2-BAYES_diagnosis_",
                         MMT,
                         "_COEF.csv",
                         sep = ""),
            sep = ",",
            col.names = TRUE,
            row.names = FALSE,
            append = FALSE)

## save R workspace to speed up process
#### uncomment to save the working space
# save.image(paste("BAYES_diagnosis_", mmt, ".RData", sep = ""))

## plotting the MCMC chain results for diagnosis
#### divergence plot of MCMC chains, chain mix
mcmc_plot(bayes_res, 
          type = "trace")

#### autocorrelation plot of MCMC samples
mcmc_plot(bayes_res, 
          type = "acf_bar")

#### estimates of coefficients
mcmc_plot(bayes_res, 
          type = "areas",
          prob = 0.95)

#### model summary
summary(bayes_res)