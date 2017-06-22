# Functions for rarefaction analysis

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

summarizeAMRbySample <- function(X) {
  amrResultsSummarised <- X %>% group_by_('Sample','LevelName') %>% 
    summarise(Hits=sum(Hits_Seen))
  return(amrResultsSummarised)
}
  
widenAMR <- function(summarizedAMR) {
  amrLevelWide <- summarizedAMR %>% 
    spread('Sample', 'Hits', fill = 0)
  return(amrLevelWide)
}

matrixAMR <- function(amrLevelWide) {
  amrLevelMat <- amrLevelWide[,2:ncol(amrLevelWide)]
  row.names(amrLevelMat) <- row.names(amrLevelWide)
  return(amrLevelMat)
}

normalizeAMR <- function(matrixAMR){
  experiment <- newMRexperiment(matrixAMR)
  cumNorm(experiment)
  extractExperiment <- data.frame(MRcounts(experiment, norm=T))
  extractExperiment <- t(extractExperiment)
  extractExperiment <- round(extractExperiment)
#row.names(extractExperiment) <- str_replace(row.names(matrixAMR), "X", "")
  return(extractExperiment)
}

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

samplesByLevel <- function(rarefiedData) {
  
  sampleNames <- attr(alphaRarefactionROP(normalizeAMR)$raw_species_abundance, "names")
  
}

ggplot_rarefy_df <- lapply(ggplot_rarefy$out, unlist)

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


# from Stack Overflow

rarec <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Genes", 
                   label = TRUE, cols = c(rep('red', nrow(x) / 2), rep('blue', nrow(x) / 2)), ...) {
  tot <- rowSums(x)
  S <- specnumber(x)
  nr <- nrow(x)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample)
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    abline(h = rare, lwd = 0.5)
  }
  for (ln in seq_len(length(out))) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = cols[ln], ...)
  }
  if (label) {
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}
