working_directory = "/Users/nick.mertens/Library/CloudStorage/OneDrive-Blend360/Consumer Reports/Documents/2023 Bayesian Modeling - Phase II/CR-Bayesian/code/8_v2 - Hierarchical Bayesian logistic probability model with Age and Miles/"
# working_directory = "C:/Users/JaeHunLee/OneDrive - Blend 360/Desktop/CR/Bayesian_git/code/8 - Hierarchical Bayesian logistic probability model with Age and Miles/"
setwd(working_directory)
getwd()

source(paste(working_directory, "bayes_hier_LP_M8_v2.R", sep = ""))


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

make = "Nissan"
iter = 5000
chains = 12
sample_num = 300
runs = 5

# run the for loop
for(run in 1:runs){
  # break up training data to remove the nissan rogue rows
  training_df = df[df$MMT != "Nissan Rogue", ]
  
  # sample from the nissan_rogue_df
  nissan_rogue_sample = nissan_rogue_df[sample(nrow(nissan_rogue_df), sample_num), ]
  
  # combine the training_df with the sample of nissan rogue
  training_df = rbind(training_df, nissan_rogue_sample)
  
  fit = run_model(training_df, make = "Nissan", iter = iter, chains = chains, problem_area = "q19_2", save_fit=FALSE)
  
  # Save the model
  saveRDS(fit, paste("models/fit_", make, "_M8_v2_", as.character(iter), "_", as.character(chains), "_rogue_ds", "_", sample_num,".rds", sep=""))
  
  # collect the prediction scores for the nissan model
  res_df = pred_prob(training_df, 
                     fit_model_name = paste("models/fit_", make, "_M8_v2_", as.character(iter), "_", as.character(chains), "_rogue_ds", "_", sample_num,".rds", sep=""), 
                     problem_area = "q19_2", 
                     coef_mode="mode")
  
  write.csv(res_df, paste("./downsample_results/", make, "_downsampling_", sample_num, "M8_v2_run", run, ".csv", sep = ""))
  
  # filter predictions to get only nissan rogue
  nissan_res_df = res_df[res_df$MMT == "Nissan Rogue", ]
  
  # append a column to record sample size
  nissan_res_df = cbind(nissan_res_df, sample_num = sample_num, run = run)
  
  # append the result to the res_df
  complete_res_df <- rbind(complete_res_df, nissan_res_df)
}

# check the file
view(complete_res_df)

# save final file
write.csv(complete_res_df, paste("./downsample_results/", make, "_downsampling_", sample_num, "_M8_v2_allruns.csv", sep = ""))

# create a list of all csv files in the working directory
files <- list.files(path = "./downsample_results/", pattern = "*.csv")

# get the files for the specific downsampled results
files_300 <- files[2:6]
files_600 <- files[8:12]

# create an empty dataframe to store the combined data
combined_data <- data.frame()

# loop through all csv files and combine them into one dataframe
# for (file in files_300) {
for (file in files_600) {
  print(file)
  data <- read.csv(paste("./downsample_results/", file, sep =""))
  combined_data <- rbind(combined_data, data)
}

# group by MakeName, MMT, and MY and calculate the mean of q19_2_mean and y_pred_mean
grouped_data <- combined_data %>% 
  group_by(MakeName, MMT, MY) %>% 
  summarise(q19_2_mean_mean = mean(q19_2_mean),
            y_pred_mean_mean = mean(y_pred_mean))

# return the final dataframe
grouped_data

