library(Matrix)
source("../RokaiApp/compute_pvalues.R")
source("../RokaiApp/rokai_core.R")
source("../RokaiApp/rokai_circuit.R")
source("../RokaiApp/rokai_weights.R")
source("../RokaiApp/rokai_kinase_weights.R")
source("../RokaiApp/rokai_inference.R")
source("../RokaiApp/compute_pvalues.R")

### Parse arguments
args <- commandArgs(trailingOnly = TRUE)

input_csv <- args[1]
output_csv <- args[2]
only_refine_phospho_profiles_str <- args[3]
datanorm <- args[4]
include_signor_str <- args[5]
include_ppi_str <- args[6]
include_sd_str <- args[7]
include_coev_str <- args[8]
network_file <- args[9]

only_refine_phospho_profiles <- tolower(only_refine_phospho_profiles_str) %in% c('true', 't')
include_signor <- tolower(include_signor_str) %in% c('true', 't')
include_ppi <- tolower(include_ppi_str) %in% c('true', 't')
include_sd <- tolower(include_sd_str) %in% c('true', 't')
include_coev <- tolower(include_coev_str) %in% c('true', 't')

### Load the network
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
phospho_data_all <- read.csv(input_csv, sep='\t')
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
    switch (datanorm,
        "Centered" = Xv <- (Xv - mean(Xv)),
        "Normalized" = Xv <- (Xv - mean(Xv)) / sd(Xv),
        "Raw" = Xv <- Xv)
    Sx <- rep(sd(Xv), length(Xv))
    ds <- (list("Xv" = Xv, "Sx" = Sx, "validSites" = validSites))
    ### Run RoKAI
    Wk2s <- NetworkData$net$Wkin2site.psp

    if(include_signor){
      Wk2s <- Wk2s | NetworkData$net$Wkin2site.signor
    }

    nSite <- ncol(Wk2s) #I think it was already set to that value but let's be on the safe side
    wk2s <- Wk2s[, validSites];
    nSubs <- (wk2s %*% rep(1, length(Xv)))

    #Add 'ppi' network
    if(include_ppi){
    Wk2k <- NetworkData$net$Wkin2kin * 1e-3
    }else{
      Wk2k <- NULL
    }

    Ws2s <- Matrix::sparseMatrix(
      i = c(),
      j = c(),
      x = TRUE,
      dims = c(nSite, nSite)
    )
    #Add 'sd' network
    if(include_sd){
    Ws2s <- Ws2s | NetworkData$net$Wsite2site.sd
    }

    #Add 'coev' network
    if(include_coev){
    Ws2s <- Ws2s | NetworkData$net$Wsite2site.coev
    }
    Ws2s <- Ws2s[validSites, validSites]
    rc <- rokai_core(Xv, Sx, wk2s, Wk2k, Ws2s)


    if (only_refine_phospho_profiles) {
      rokai_result_experiment <- data.frame(phospho_data[valids, 'Site'][order(indices[valids])], rc$Xs)
      names(rokai_result_experiment) <- c('Site', experiment)
    }else {
      #Run the kinase activity inference part of RoKAI
      Fk <- rokai_kinase_weights(Xv, wk2s, rc$F)
      ri <- rokai_inference(Xv, Sx, Fk)
      A <- ri$A
      S <- ri$S
      Z <- ri$Z
      res <- compute_pvalues(as.matrix(Z))
      K <- NetworkData$Kinase
      #Define Column names experiment-specific
      activity_col <- paste0('Activity (', experiment, ')')
      std_err_col <- paste0('StdErr (', experiment, ')')
      zscore_col <- paste0('ZScore (', experiment, ')')
      fdr_col <- paste0('FDR (', experiment, ')')
      K[activity_col] <- as.matrix(A)
      K[std_err_col] <- as.matrix(S)
      K[zscore_col] <- as.matrix(Z)
      K[fdr_col] <- res$QValues
      isPhosphatase <- K$Type == "Phosphatase"
      K[isPhosphatase, activity_col] <- -1 * K[isPhosphatase, activity_col]
      K[isPhosphatase, zscore_col] <- -1 * K[isPhosphatase, zscore_col]

      rokai_result_experiment <- K[complete.cases(K), c('Gene', activity_col, std_err_col, zscore_col, fdr_col)]
    }
    rokai_result_all[[length(rokai_result_all) + 1]] <- rokai_result_experiment
  }, error = function(e) e)
}

if (only_refine_phospho_profiles) {
  rokai_result_singledf <- Reduce(function(x, y) merge(x, y, by = 'Site', all = TRUE), rokai_result_all)
} else {
  rokai_result_singledf <- Reduce(function(x, y) merge(x, y, by = 'Gene', all = TRUE), rokai_result_all)
}
write.table(rokai_result_singledf, output_csv, quote = F, row.names = F, sep='\t')