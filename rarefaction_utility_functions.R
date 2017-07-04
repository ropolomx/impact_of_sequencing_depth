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
#' @param rarefiedData 
#'
#' @return
#' @export
#'
#' @examples
samplesByLevel <- function(rarefiedData) {
  
  sampleNames <- attr(alphaRarefactionROP(normalizeAMR)$raw_species_abundance, "names")
  
}

