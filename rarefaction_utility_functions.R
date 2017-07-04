# Utility functions for rarefaction analysis

#' Summarize AMR results by levels
#' 
#' @param X A tidy dataset.
#' @param amrLevel A level of the MEGARes database. For example: class, group 
#' or mechanism
#' 
#' @return Summarized results of \code{X} by \code{amrLevel} 
#' 
#'
#' @examples 

summarizeAMRlevels <- function(X, amrLevel) {
  amrResultsSummarised <- X %>% group_by_('Sample', amrLevel) %>% 
    summarise(Hits=sum(Hits_Seen))
  return(amrResultsSummarised)
}

#' Title
#'
#' @param X 
#'
#' @return
#' @export
#'
#' @examples Summarized results of \code{X} by \code{}
 
summarizeAMRbySample <- function(X) {
  amrResultsSummarised <- X %>% 
    group_by_('Sample','LevelName') %>% 
    summarise(Hits=sum(Hits_Seen))
  return(amrResultsSummarised)
}
  
#' Title
#'
#' @param summarizedAMR 
#'
#' @return
#' @export
#'
#' @examples

widenAMR <- function(summarizedAMR) {
  amrLevelWide <- summarizedAMR %>% 
    spread('Sample', 'Hits', fill = 0)
  return(amrLevelWide)
}

#' Title
#'
#' @param amrLevelWide 
#'
#' @return
#' @export
#'
#' @examples
#' 
matrixAMR <- function(amrLevelWide) {
  amrLevelMat <- amrLevelWide[,2:ncol(amrLevelWide)]
  row.names(amrLevelMat) <- row.names(amrLevelWide)
  return(amrLevelMat)
}


#' Title
#'
#' @param matrixAMR 
#'
#' @return
#' @export
#'
#' @examples
normalizeAMR <- function(matrixAMR){
  experiment <- newMRexperiment(matrixAMR)
  cumNorm(experiment)
  extractExperiment <- data.frame(MRcounts(experiment, norm=T))
  extractExperiment <- t(extractExperiment)
  extractExperiment <- round(extractExperiment)
#row.names(extractExperiment) <- str_replace(row.names(matrixAMR), "X", "")
  return(extractExperiment)
}

#' Title
#'
#' @param X 
#' @param step 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
alphaRarefactionROP <- function(X, step=0.05, method='invsimpson') {
    S <- specnumber(X, MARGIN=2)
    raremax <- min(colSums(X))
    tot <- colSums(X)
    nr <- ncol(X)
    rarefy_out <- lapply(seq_len(nr), function(i){
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
                rarefy_out=rarefy_out))
} 

#' Title
#'
#' @param rarefiedData 
#'
#' @return
#' @export
#'
#' @examples
samplesByLevel <- function(rarefiedData) {
  
  sampleNames <- attr(alphaRarefactionROP(normalizeAMR)$raw_species_abundance, "names")
  
}

ggplot_rarefy_df <- lapply(rarecurve_ROP$out, unlist)

ggplot_rarefy_df <- lapply(ggplot_rarefy_df, data.frame)

rarefactionDF <- do.call("rbind",ggplot_rarefy_df[1:32])
names(rarefactionDF) <- "Classes"

rarefactionDF$Samples <- rep(class_samples, each=21)
rarefactionDF$Sampling_size <- str_replace(row.names(rarefactionDF), "N","")
rarefactionDF$Sampling_size <- as.numeric(as.character(rarefactionDF$Sampling_size))

rarefactionDF$Sample_type <- str_replace(rarefactionDF$Samples, "\\d+_", "")

rarefactionDF$Sample_type <- str_replace(rarefactionDF$Sample_type, "\\d$","")

rarefactionDF$Sample_number <- str_replace(rarefactionDF$Samples, "_.*", "")

halves <- rarefactionDF %>% filter(grepl("half",Sample_type))

rarefactionDFGroup <- rarefactionDF %>% group_by(Sample_number, Sample_type) %>% summarise(meanOTUs = mean(Classes), meanSamples=mean(Sampling_size))

ggplot(rarefactionDFGroup, aes(meanSamples, meanOTUs)) + geom_line() + facet_grid(. ~ Sample_type)

ggplot(rarefactionDF, aes(Sampling_size,Classes)) + geom_line(aes(color=Samples)) + facet_grid(Sample_type~ Sample_number)
