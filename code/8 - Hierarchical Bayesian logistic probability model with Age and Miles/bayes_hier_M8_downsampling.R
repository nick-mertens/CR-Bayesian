#working_directory = "/Users/nick.mertens/Library/CloudStorage/OneDrive-Blend360/Consumer Reports/Documents/2023 Bayesian Modeling - Phase II/CR-Bayesian/code/8 - Hierarchical Bayesian logistic probability model with Age and Miles/"
working_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code/8 - Hierarchical Bayesian logistic probability model with Age and Miles/"
setwd(working_directory)

source(paste(working_directory, "bayes_hier_LP_M8.R", sep = ""))


# TO DO:
# - for loop where input is the number of samples
# - extract Nissan Rogue from data and sample it
# - append the sampled Nissan Rogue back to training data
# - run code as normal
# - extract Nissan Rogue predictions and append to new dataframe
# - return the new dataframe

# set a seed for replicating data
set.seed(1231)

# initialize empty res_df
complete_res_df = data.frame()

# break up the data to sample from Nissan Rogue group
nissan_rogue_df = df[df$MMT == "Nissan Rogue", ]

# set up sampling parameter
sample_list = c(300, 200, 100)
#sample_list = c(600)

make = "Nissan"
iter = 5000
chains = 12

# run the for loop
for(sample_num in sample_list){
  # break up training data to remove the nissan rogue rows
  training_df = df[df$MMT != "Nissan Rogue", ]
  
  # sample from the nissan_rogue_df
  nissan_rogue_sample = nissan_rogue_df[sample(nrow(nissan_rogue_df), sample_num), ]
  
  # combine the training_df with the sample of nissan rogue
  training_df = rbind(training_df, nissan_rogue_sample)
  
  # run the model 
  fit = run_model(training_df, make = "Nissan", iter=iter, chains=chains, problem_area = "q19_2", save_fit=FALSE)
  
  # Save the model
  saveRDS(fit, paste("models/fit_", make, "_M8_", as.character(iter), "_", as.character(chains), "_rogue_ds", "_", sample_num,".rds", sep=""))
  
  # collect the prediction scores for the nissan model
  res_df = pred_prob(training_df, 
                     fit_model_name =  paste("models/fit_", make, "_M8_", as.character(iter), "_", as.character(chains), "_rogue_ds", "_", sample_num,".rds", sep=""), 
                     problem_area = "q19_2", 
                     coef_mode="mode")
  
  # filter predictions to get only nissan rogue
  nissan_res_df = res_df[res_df$MMT == "Nissan Rogue", ]
  
  # append a column to record sample size
  nissan_res_df = cbind(nissan_res_df, sample_num = sample_num)
  
  # append the result to the res_df
  complete_res_df <- rbind(complete_res_df, nissan_res_df)
}

# check the file
view(complete_res_df)

# write file to csv
write.csv(complete_res_df, "nissan_rogue_downsampling.csv")

