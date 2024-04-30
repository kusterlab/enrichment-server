library(Matrix)
source("../RokaiApp/compute_pvalues.R")
source("../RokaiApp/rokai_core.R")
source("../RokaiApp/rokai_circuit.R")
source("../RokaiApp/rokai_weights.R")

### Parse arguments
args <- commandArgs(trailingOnly = TRUE)

input_csv <- args[1]
output_csv <- args[2]

### Load the network
network_file <- '../RokaiApp/data/rokai_network_data_uniprotkb_human.rds'
NetworkData <- readRDS(network_file)
NetworkData$Kinase$Type <- "Kinase"
nKinase <- nrow(NetworkData$Kinase)
nSite <- nrow(NetworkData$Site)

NetworkData$net$Wkin2site.depod <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                                        dims = c(nKinase, nSite))
Phosphatase <- data.frame(
  KinaseID = NetworkData$Phosphatase$ID,
  KinaseName = paste("Phospha-", NetworkData$Phosphatase$Gene, sep = ""),
  Gene = NetworkData$Phosphatase$Gene,
  Type = "Phosphatase"
)

nKinase <- nrow(NetworkData$Kinase)
nPhosphatase <- nrow(Phosphatase)
Wphospha2site <- NetworkData$net$Wphospha2site
NetworkData$Kinase <- rbind(NetworkData$Kinase, Phosphatase)
NetworkData$Wkin2site <- rbind(NetworkData$Wkin2site, Wphospha2site)
NetworkData$net$Wkin2site <- rbind(NetworkData$net$Wkin2site, Wphospha2site)
NetworkData$net$Wkin2site.psp <- rbind(NetworkData$net$Wkin2site.psp, Wphospha2site)
NetworkData$net$Wkin2site.psp.base <- rbind(NetworkData$net$Wkin2site.psp.base, Wphospha2site)
NetworkData$net$Wkin2site.signor <- rbind(NetworkData$net$Wkin2site.signor, Wphospha2site)
NetworkData$net$Wkin2kin <- NetworkData$net$Wkin2kin.phospha
Wphospha2kinx <- Matrix::sparseMatrix(i = 1:nPhosphatase, j = nKinase + (1:nPhosphatase), dims = c(nPhosphatase, nKinase + nPhosphatase))
NetworkData$net$Wkin2site.depod <- (Matrix::t(Wphospha2kinx) %*% NetworkData$net$Wphospha2site)

### Parse the input csv file
phospho_data_all <- read.csv(input_csv)
experiment_names <- colnames(phospho_data_all)[2:length(colnames(phospho_data_all))]

phospho_data_all$ID <- gsub('_\\D', '_', phospho_data_all$Site)

rokai_result_all <- list()
for (experiment in experiment_names) {
  #Fix: Skip experiments with too few valid values
  tryCatch({
    #Preprocess
    valids <- !is.na(phospho_data_all[, experiment])
    phospho_data <- phospho_data_all[valids, c('Site', 'ID', experiment)]
    indices <- match(phospho_data$ID, NetworkData$Site$Identifier)
    valids <- !is.na(indices);
    X <- rep(NA, nrow(NetworkData$Site))
    X[indices[valids]] <- phospho_data[valids, experiment]
    validSites <- !is.na(X)
    Xv <- X[validSites]
    #Normalize
    Xv <- (Xv - mean(Xv)) / sd(Xv)
    Sx <- rep(sd(Xv), length(Xv))
    ds <- (list("Xv" = Xv, "Sx" = Sx, "validSites" = validSites))
    ### Run RoKAI
    Wk2s <- NetworkData$net$Wkin2site.psp
    nSite <- ncol(Wk2s) #I think it was already set to that value but let's be on the safe side
    wk2s <- Wk2s[, validSites];
    nSubs <- (wk2s %*% rep(1, length(Xv)))
    #Use all components
    ropts <- list("ppi" = TRUE, "sd" = TRUE, "coev" = TRUE)
    Wk2k <- NetworkData$net$Wkin2kin * 1e-3
    Ws2s <- Matrix::sparseMatrix(
      i = c(),
      j = c(),
      x = TRUE,
      dims = c(nSite, nSite)
    )
    Ws2s <- Ws2s | NetworkData$net$Wsite2site.sd
    Ws2s <- Ws2s | NetworkData$net$Wsite2site.coev
    Ws2s <- Ws2s[validSites, validSites]
    rc <- rokai_core(Xv, Sx, wk2s, Wk2k, Ws2s)
    rokai_result_experiment <- data.frame(phospho_data[valids, 'Site'][order(indices[valids])], rc$Xs)
    names(rokai_result_experiment) <- c('Site', experiment)
    rokai_result_all[[length(rokai_result_all) + 1]] <- rokai_result_experiment
  }, error = function(e) e)
}

rokai_result_singledf <- Reduce(function(x, y) merge(x, y, by = 'Site', all = TRUE), rokai_result_all)

write.csv(rokai_result_singledf, output_csv, quote = F, row.names = F)