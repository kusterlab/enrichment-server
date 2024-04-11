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
  n_steps_pruning = 3,
  netObj = read.csv(file = pkn_path),
  carnival_options = carnival_options)


phonemes_result_protein <- PHONEMeS::get_protein_network(phonemes_result)

readr::write_csv(phonemes_result_protein$weightedSIF, output_path)