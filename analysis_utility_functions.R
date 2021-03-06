#' Summarize AMR results by levels
#' 
#' @param X A tidy AMR Coverage Sampler dataset.
#' @param amrLevel An AMR level of the MEGARes database. For example: class, group 
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

#' Summarize AMR results by depth
#'
#' @param X A tidy dataset with AMR results from Coverage Sampler
#'
#' @return Summarized results of \code{X} by AMR category
#' @export
#'
#' @examples 
 
summarizeAMRbyCategory <- function(X) {
  amrResultsSummarised <- X %>% 
    group_by_('Sample', 'CategoryName') %>% 
    summarise(SumHits=sum(Hits))
  return(amrResultsSummarised)
}

#' Summarize Kraken by depth
#'
#' @param X A tidy Kraken report dataset 
#'
#' @return Summarized dataframe with the sum of hits by Sample and TaxID
#' @export
#'
#' @examples
 
summarizeKrakenbyTaxID <- function(X) {
  krakenResultsSummarised <- X %>% 
    group_by_('Sample', 'TaxID') %>% 
    summarise(SumHits=sum(CladeReads))
  return(krakenResultsSummarised)
}
  
#' Title
#'
#' @param summarizedAMR 
#'
#' @return Summarized AMR results in wide format
#' @export
#'

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
#' @param amrLevelWide 
#'
#' @return
#' @export
#'
#' @examples
#' 
matrixAMRanalytical <- function(amrLevelWide) {
  amrLevelMat <- as.matrix(amrLevelWide[,2:ncol(amrLevelWide)])
  row.names(amrLevelMat) <- amrLevelWide$Gene
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
  krakenLevelMat <- as.matrix(krakenLevelWide[,c(2:ncol(krakenLevelWide))])
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
    rarefactionCurve <- ggplot(taxSubset, aes(Subsample, value, color=Depth)) +
      geom_point(aes(group=Sample), alpha=0.8, size=4) + 
      theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=42),
          legend.position="right",
          legend.title=element_text(size=34),
          legend.text=element_text(size=34, vjust=0.5),
          legend.key.size = unit(2, "lines"),
          legend.spacing = unit(0.2,"lines"),
          plot.title=element_text(size=52, hjust=0.5)) +
    xlab("\nNumber of reads") +
    ylab(paste('Number ', 'of ', taxSubset$krakenLevel, '\n')) +
    scale_color_manual(values=rev(cbPalette)) +
    #scale_y_log10() +
    scale_x_continuous(labels=scientific)
    #facet_grid(. ~ Depth)
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
    geom_point(aes(group=SampleID), alpha=0.2, size=4) + 
    theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=42),
          legend.position="right",
          legend.title=element_text(size=36),
          legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=52, hjust=0.5)) +
      xlab("\nNumber of reads") +
      ylab(paste('Number ', 'of ', amrSubset$AMRLevel, '\n')) +
      ggtitle(paste(amrSubset$AMRLevel, '\n')) +
      scale_color_manual(values=rev(cbPalette)) +
      scale_x_continuous(labels=scientific) 
      #scale_x_log10() +
      #facet_grid(. ~ Depth)
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
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, InvSimpson, color=Depth)) +
    geom_boxplot(size=1.25) +
    theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=38),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=42),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Inverse Simpson's Index\n") +
    #ggtitle('Alpha Diversity by Depth for Rarefied data\nInverse Simpson Index') +
    scale_color_manual(values=rev(cbPalette)) +
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
  alphaDivBoxPlot<- ggplot(amrSubset, aes(Depth, InvSimpson, color=Depth)) +
    geom_boxplot(size=1.25) +
    theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=42),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Inverse Simpson's Index\n") + 
    #ggtitle('Alpha Diversity by Depth for Rarefied data\nInverse Simpson Index') +
    scale_color_manual(values=rev(cbPalette)) +
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
    geom_boxplot(size=1.25)+ 
    theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=42),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5)) +
    xlab("Depth") +
    ylab("Number of Unique\nAMR Categories\n") +
    #ggtitle('AMR Category Richness by Depth for Raw Data') +
    scale_color_manual(values=rev(cbPalette)) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}

#' Title
#'
#' @param taxSubset 
#'
#' @return
#' @export
#'
#' @examples
krakenRawSpeciesRich <- function(taxSubset){
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, RawSpeciesAbundance, color=Depth)) +
    geom_boxplot(size=1.25) + 
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5),
          panel.background = element_rect(fill = "grey90", colour = "grey80")) +
    xlab("Depth") +
    ylab("Number of Unique Taxa\n") +
    #ggtitle('Species Richness by Depth for Raw Data\n') +
    scale_color_manual(values=rev(cbPalette)) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}

#' alpha_rarefaction: function that returns species count, 
#' rarefied species count, and alpha diversity measures
#' for each sample in the m x n matrix, m = features, n = samples
#' 
#' Written by Steven Lakin (Colorado State University)
#'
#' @param X Microbiome or resistome results in matrix format
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


#' Boxplots of Pielou's evenness with Kraken results
#'
#' @param taxSubset 
#'
#' @return
#' @export
#'
#' @examples
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

#' krakenShannon: Boxplots of Shannon's Index of diversity calculated for
#' Kraken results
#'
#' @param taxSubset A dataset of Kraken hits for a given taxon
#'
#' @return
#' @export
#'
#' @examples
krakenShannon <- function(taxSubset){
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, Shannon, color=Depth)) +
    geom_boxplot(size=1.25) + 
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=38),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5),
          panel.background = element_rect(fill = "grey90", colour = "grey80")) +
    xlab("Depth") +
    ylab("Shannon's Index\n") +
    #ggtitle('Alpha Diversity by Depth for Normalized data\nShannon Index of Diversity') +
    scale_color_manual(values=rev(cbPalette)) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}

#' amrShannon: Boxplots of Shannon's Index of diversity calculated for
#' AMR results
#'
#' @param taxSubset A dataset of AMR hits for a given taxon
#'
#' @return
#' @export
#'
#' @examples
amrShannon <- function(taxSubset){
  alphaDivBoxPlot<- ggplot(taxSubset, aes(Depth, Shannon, color=Depth)) +
    geom_boxplot(size=1.25) + 
    theme(strip.text.x=element_text(size=30),
          axis.text.y=element_text(size=38),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=44),
          axis.title.y=element_text(size=44),
          legend.position="none",
          #legend.title=element_text(size=36),
          #legend.text=element_text(size=36, vjust=0.5),
          plot.title=element_text(size=50, hjust=0.5),
          panel.background = element_rect(fill = "grey90", colour = "grey80")) +
    xlab("Depth") +
    ylab("Shannon's Index\n") +
    #ggtitle('Alpha Diversity by Depth for Normalized data\nShannon Index of Diversity') +
    scale_color_manual(values=rev(cbPalette)) +
    facet_wrap(~ Level, nrow=2, scales = "free_y")
}

#' amrScaledNonSmooth: AMR rarefaction curves from Rarefaction Analyzer data
#'
#' @param amrLevel 
#'
#' @return
#' @export
#'
#' @examples
amrScaledNonSmooth <- function(amrLevel){
    rarefactionCurve <- ggplot(amrLevel, aes(ScaledReads, Counts, color=Depth)) +
      geom_line(aes(group=SampleID), alpha=0.6,size=3) + 
      theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=32),
          axis.text.x=element_text(size=32, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=36, vjust = -0.5),
          legend.position="right",
          legend.title=element_text(size=34),
          legend.text=element_text(size=34, vjust=0.5),
          legend.key.size = unit(2, "lines"),
          legend.spacing = unit(0.2,"lines"),
          plot.title=element_text(size=52, hjust=0.5)) +
    xlab("\nNumber of reads") +
    ylab(paste('Number ', 'of ', amrLevel$Level, '\n')) +
    scale_color_manual(values=rev(cbPalette)) +
      scale_x_continuous(labels=scientific) 
}

#' Title
#'
#' @param taxSubset 
#'
#' @return
#' @export
#'
#' @examples
krakenRarCurveNML <- function(taxSubset){
    rarefactionCurve <- ggplot(taxSubset, aes(Reads, taxonCount, color=Depth)) +
      geom_line(aes(group=Sample), alpha=0.5, size=4) + 
      #geom_point(aes(group=Sample), size=4) +
      theme(strip.text.x=element_text(size=35),
          axis.text.y=element_text(size=40),
          axis.text.x=element_text(size=35, angle=90, vjust=0.3),
          axis.title.x=element_text(size=42),
          axis.title.y=element_text(size=42),
          ylim(0,40),
          legend.position="right",
          legend.title=element_text(size=34),
          legend.text=element_text(size=34, vjust=0.5),
          legend.key.size = unit(2, "lines"),
          legend.spacing = unit(0.2,"lines"),
          plot.title=element_text(size=52, hjust=0.5)) +
    xlab("\nNumber of reads") +
    ylab(paste('Number ', 'of ', taxSubset$Level, '\n')) +
    scale_color_manual(values=rev(cbPalette)) +
    #scale_y_log10() +
    scale_x_continuous(labels=scientific)
    #facet_grid(. ~ Depth)
}
