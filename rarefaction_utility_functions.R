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

#' Summarize by depth
#'
#' @param X 
#'
#' @return
#' @export
#'
#' @examples Summarized results of \code{X} by sample depth
 
summarizeAMRbyDepth <- function(X) {
  amrResultsSummarised <- X %>% 
    group_by_('Sample', 'CategoryName') %>% 
    summarise(SumHits=sum(Hits))
  return(amrResultsSummarised)
}

#' Summarize by depth
#'
#' @param X 
#'
#' @return
#' @export
#'
#' @examples Summarized results of \code{X} by sample depth
 
summarizeKrakenbyDepth <- function(X) {
  krakenResultsSummarised <- X %>% 
    group_by_('Sample', 'TaxID') %>% 
    summarise(SumHits=sum(CladeReads))
  return(krakenResultsSummarised)
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
    spread_('Sample', 'SumHits', fill = 0)
  return(amrLevelWide)
}

#' Title
#'
#' @param summarizedKraken
#'
#' @return
#' @export
#'
#' @examples

widenKraken <- function(summarizedKraken) {
  krakenLevelWide <- summarizedKraken %>% 
    spread_('Sample', 'SumHits', fill = 0)
  return(krakenLevelWide)
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
  row.names(amrLevelMat) <- amrLevelWide$CategoryName
  return(amrLevelMat)
}
#' Title
#'
#' @param krakenLevelWide 
#'
#' @return 
#' @export
#'
#' @examples
#' 
#' 
matrixKraken <- function(krakenLevelWide) {
  krakenLevelMat <- krakenLevelWide[,c(2:33)]
  row.names(krakenLevelMat) <- krakenLevelWide$TaxID
  return(krakenLevelMat)
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

#' Title
#'
#' @param taxSubset 
#'
#' @return
#' @export
#'
#' @examples
krakenRarefactionCurve <- function(taxSubset){
    rarefactionCurve <- ggplot(taxSubset, aes(Number_of_Reads, Counts, color=Depth)) +
    geom_point(alpha=1/50, size=2) + 
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=52, hjust=0.5)) +
    xlab("\nNumber of reads") +
    ylab(paste('Number ', 'of ', taxSubset$krakenLevel, '\n')) +
    ggtitle(paste(taxSubset$krakenLevel, '\n')) +
    scale_color_manual(values=vennPalette) +
    scale_x_continuous(label=scientific) +
    facet_grid(. ~ Depth)
}

#' Title
#'
#' @param amrSubset 
#'
#' @return
#' @export
#'
#' @examples
amrRarefactionCurve <- function(amrSubset){
    rarefactionCurve <- ggplot(amrSubset, aes(Subsample, value, color=Depth)) +
    geom_line(aes(group=SampleID, alpha=0.4, size=4)) + 
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=52, hjust=0.5)) +
      xlab("\nNumber of reads") +
      ylab(paste('Number ', 'of ', amrSubset$AMRLevel, '\n')) +
      ggtitle(paste(amrSubset$AMRLevel, '\n')) +
      scale_color_manual(values=vennPalette) +
      scale_x_continuous(label=scientific) +
      facet_grid(. ~ Depth)
}


#' Title
#'
#' @param taxSubset 
#'
#' @return
#' @export
#'
#' @examples

krakenAlphaDiv <- function(taxSubset){
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, AlphaDiv, color=Depth)) +
    geom_boxplot() +
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Inverse Simpson's Index") +
    ggtitle('Alpha Diversity by Depth for Rarefied data\nInverse Simpson Index') +
    scale_color_manual(values=vennPalette) +
    facet_wrap( ~ Level, nrow=2, scales = "free_y")
}

#' Title
#'
#' @param taxSubset 
#'
#' @return
#' @export
#'
#' @examples
amrAlphaDiv <- function(amrSubset){
  alphaDivBoxPlot<- ggplot(amrSubset, aes(Depth, AlphaDiv, color=Depth)) +
    geom_boxplot() +
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Inverse Simpson's Index") +
    ggtitle('Alpha Diversity by Depth for Rarefied data\nInverse Simpson Index') +
    scale_color_manual(values=vennPalette) +
    facet_wrap( ~ Level, nrow=2, scales = "free_y")
}





#' Title
#'
#' @param taxSubset 
#'
#' @return
#' @export
#'
#' @examples
amrRawSpeciesRich <- function(amrSubset){
  alphaDivBoxPlot<- ggplot(amrSubset, aes(Depth, RawSpeciesAbundance, color=Depth)) +
    geom_boxplot(size=1) + 
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Unique Taxa") +
    ggtitle('Species Richness by Depth for Raw Data') +
    scale_color_manual(values=vennPalette) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}

krakenRawSpeciesRich <- function(taxSubset){
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, RawSpeciesAbundance, color=Depth)) +
    geom_boxplot(size=1) + 
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Unique Taxa\n") +
    ggtitle('Species Richness by Depth for Raw Data\n') +
    scale_color_manual(values=vennPalette) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}

#' alpha_rarefaction: function that returns species count, 
#' rarefied species count, and alpha diversity measures
#' for each sample in the m x n matrix, m = features, n = samples
#' 
#' Written by Steven Lakin Colorado State University
#'
#' @param X 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
alpha_rarefaction <- function(X, minlevel, method='invsimpson') {
  S <- specnumber(X, MARGIN=1)
  raremax <- min(rowSums(X))
  if( raremax < minlevel ) raremax <- minlevel
  Srare <- rarefy(X, raremax, MARGIN=1)
  Xrare <- t(rrarefy(t(X), raremax))
  alphadiv <- diversity(Xrare, index=method, MARGIN=1)
  H <- diversity(Xrare, index="shannon", MARGIN=1)
  J <- H/log(S)
  return(list(raw_species_abundance=S,
              rarefied_species_abundance=Srare,
              rarefied_data=Xrare,
              alphadiv=alphadiv,
              shannon=H,
              evenness=J))
}


krakenEvenness <- function(taxSubset){
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, Evenness, color=Depth)) +
    geom_boxplot(size=1) +
    geom_point() +
    geom_jitter() +
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Pielou's evenness") +
    ggtitle('Evenness by Depth for Rarefied data\nInverse Simpson Index') +
    scale_color_manual(values=vennPalette) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}

krakenShannon <- function(taxSubset){
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, Shannon, color=Depth)) +
    geom_boxplot(size=1) +
    geom_point() +
    geom_jitter() +
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Shannon's index of diversity") +
    ggtitle('Alpha Diversity by Depth for Rarefied data\nShannon Index of Diversity') +
    scale_color_manual(values=vennPalette) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}
