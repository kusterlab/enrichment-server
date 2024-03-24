
#TODO: Untested code!

args <- commandArgs(trailingOnly = TRUE)
sites_path <- args[1]
targets_path <- args[2]
output_path <- args[3]

pkn_path <- '../db/phonemes_PKN_KSN.csv'




carnival_options <- PHONEMeS::default_carnival_options(solver = "cplex")
carnival_options$solverPath <- "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/cplex"

phonemes_result <- PHONEMeS::run_phonemes(inputObj = read.csv(file = targets_path),
                                  measObj = read.csv(file = sites_path),
                                  n_steps_pruning = 3,
                                  netObj = read.csv(file = pkn_path),
                                  carnival_options = carnival_options)


#TODO: What happens if I just skip this? Does this just attach the stuff that I don't want?
#phonemes_result_pps <- PHONEMeS::reattach_psites(phonemes_result)

#TODO: Does it have to be the weighted SIF? Can I use an unweighted one instead, I don't need the weights?
readr::write_csv(phonemes_result_pps$res$weightedSIF, output_path)

#TODO: What does this do? How is it different from leaving out the reattach_psites? Could I just use it?
phonemes_result_protein <- get_protein_network(phonemes_result)

