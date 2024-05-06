args <- commandArgs(trailingOnly = TRUE)
sites_path <- args[1]
targets_path <- args[2]
output_path <- args[3]

pkn_path <- '../db/phonemesPKN.csv'
# pkn_path <- '../db/phonemes_PKN_KSN.csv'


carnival_options <- PHONEMeS::default_carnival_options(solver = "cplex")

carnival_options$solverPath <- '../CPLEX/cplex'

targets_df <- read.csv(targets_path, header = FALSE)
targets_vector <- setNames(targets_df[[2]], targets_df[[1]])

sites_df <- read.csv(sites_path)
sites_vector <- setNames(sites_df[[2]], sites_df[[1]])


phonemes_result <- PHONEMeS::run_phonemes(
  inputObj = targets_vector,
  measObj = sites_vector,
  n_steps_pruning = 2,
  netObj = read.csv(file = pkn_path),
  carnival_options = carnival_options)

#We only want to return a network of proteins. Therefore, replace all p-sites by the proteins they sit on
phonemes_result$res$weightedSIF$Node2 <- sub("_.*", "", phonemes_result$res$weightedSIF$Node2)
#Drop duplicates and limit to the columns that we actually use
phonemes_result$res$weightedSIF <- dplyr::distinct(phonemes_result$res$weightedSIF,Node1,Node2)
#Drop Self-Links, PTMNavigator cannot display them
phonemes_result$res$weightedSIF <- dplyr::filter(phonemes_result$res$weightedSIF, Node1 != Node2)

readr::write_csv(phonemes_result$res$weightedSIF, output_path)