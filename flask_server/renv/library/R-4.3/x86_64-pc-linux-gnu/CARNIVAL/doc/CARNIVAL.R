## ----setup, echo=FALSE, results='hide', warning=FALSE, error=FALSE, message=FALSE, cache=FALSE----
library(knitr)
opts_chunk$set(
  cache = FALSE,
  echo = TRUE,
  warning = FALSE,
  error = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(CARNIVAL)
load(file = system.file("toy_perturbations_ex1.RData",
                        package = "CARNIVAL"))
load(file = system.file("toy_measurements_ex1.RData",
                        package = "CARNIVAL"))
load(file = system.file("toy_network_ex1.RData",
                        package = "CARNIVAL"))

carnivalOptions <- defaultLpSolveCarnivalOptions()

# Output dir
dir.create("ToyExample1", showWarnings = FALSE)
carnivalOptions$outputFolder <- "ToyExample1"

# lpSolve
resultsLpSolve <- runVanillaCarnival( perturbations = toy_perturbations_ex1, 
                                      measurements = toy_measurements_ex1,
                                      priorKnowledgeNetwork = toy_network_ex1, 
                                      carnivalOptions = carnivalOptions)


## -----------------------------------------------------------------------------
print(resultsLpSolve)

## -----------------------------------------------------------------------------
saveRDS(resultsLpSolve, "ToyExample1/network_solution.Rds")

## -----------------------------------------------------------------------------
load(file = system.file("toy_measurements_ex2.RData",
                        package="CARNIVAL"))
load(file = system.file("toy_network_ex2.RData",
                        package="CARNIVAL"))

carnivalOptions <- defaultLpSolveCarnivalOptions()

# Output dir
dir.create("ToyExample2", showWarnings = FALSE)
carnivalOptions$outputFolder <- "ToyExample2"

# lpSolve
toy_network_ex2 <- as.data.frame(toy_network_ex2)
colnames(toy_network_ex2) <- c("source", "interaction", "target")
toy_measurements_ex2 <- c("M1" = 1, "M2" = 1, "M3" = 1)

inverseCarnivalResults <- runInverseCarnival( measurements = c(toy_measurements_ex2), 
                                              priorKnowledgeNetwork = toy_network_ex2, 
                                              carnivalOptions = carnivalOptions)

## -----------------------------------------------------------------------------
print(inverseCarnivalResults)

## -----------------------------------------------------------------------------

load(file = system.file("toy_perturbations_ex1.RData",
                        package = "CARNIVAL"))
load(file = system.file("toy_measurements_ex1.RData",
                        package = "CARNIVAL"))
load(file = system.file("toy_network_ex1.RData",
                        package = "CARNIVAL"))

carnivalOptions <- defaultLpSolveCarnivalOptions()

# Output dir
dir.create("ToyExample3", showWarnings = FALSE)
carnivalOptions$outputFolder <- "ToyExample3"

# lpSolve
generateLpFileCarnival(perturbations = toy_perturbations_ex1,
                       measurements = toy_measurements_ex1, 
                       priorKnowledgeNetwork = toy_network_ex1, 
                       carnivalOptions = carnivalOptions)


## -----------------------------------------------------------------------------

lpFilename <- "toy_files_vignettes/lpFile_t18_01_53d28_04_2021n4.lp"
parsedDataFile <- "toy_files_vignettes/parsedData_t18_01_53d28_04_2021n4.RData"
resultsFromLp <- runFromLpCarnival(lpFile = lpFilename,
                                   parsedDataFile = parsedDataFile,
                                   carnivalOptions = defaultLpSolveCarnivalOptions())

## -----------------------------------------------------------------------------
resultsFromLp

## -----------------------------------------------------------------------------
sessionInfo()

