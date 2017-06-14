summarizeAMRlevels <- function(X, amrLevel) {
  amrResultsSummarised <- X %>% group_by_('Sample', amrLevel) %>% summarise(Hits=sum(Hits_Seen))
  return(amrResultsSummarised)
}

widenAMR <- function(summarizedAMR) {
  amrLevelWide <- summarizedAMR %>% spread_('Sample', 'Hits', fill = 0)
  return(amrLevelWide)
}

matrixAMR <- function(amrLevelWide) {
  amrLevelMat <- amrLevelWide[,2:ncol(widenAMR)]
  return(amrLevelMat)
}

amr_experiment <- function(amrLevelMat)
  

normalizeAMR <- function(matrixAMR){
  amr_class_experiment <- newMRexperiment(amr_results_class_wide)
  amr_class_norm <- data.frame(MRcounts(amr_class_experiment, norm=T))
  amr_class_norm <- t(amr_class_norm)
  row.names(amr_class_norm) <- str_replace(row.names(amr_class_norm), "X", "")
}

alphaRarefactionROP <- function(X, step=0.05, method='invsimpson') {
    S <- specnumber(X, MARGIN=2)
    raremax <- min(colSums(X))
    tot <- colSums(X)
    nr <- ncol(X)
     <- lapply(seq_len(nr), function(i){
      n <- seq(1, tot[i], by = tot[i]*step)
      if (n[length(n)] != tot[i])
        n <- c(n, tot[i])
      drop(rarefy(X[i, ], n))
      })
    Srare <- rarefy(X, raremax, MARGIN=2)
    Xrare <- t(rrarefy(t(X), raremax))
    alphadiv <- diversity(Xrare, index=method, MARGIN=2)
    return(list(raw_species_abundance=S,
                #rarefied_species_abundance=Srare,
                #rarefied_data=Xrare,
                alphadiv=alphadiv,
                =out))
} 
