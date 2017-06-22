library(dplyr)
library(vegan)

amrResultsClass <- amrResults %>% 
  group_by(Sample, Class) %>% 
  summarise(Hits=sum(Hits_Seen))

amrResultsClass$Depth <- str_replace(amrResultsClass$Sample, "\\d+_","")
amrResultsClass$Depth <- str_replace(amrResultsClass$Depth, "\\d$","")
amrResultsClass$Depth <- str_replace(amrResultsClass$Depth, "full", "D1")
amrResultsClass$Depth <- str_replace(amrResultsClass$Depth, "half", "D0.5")
amrResultsClass$Depth <- str_replace(amrResultsClass$Depth, "quar*", "D0.25")

krakenPhylumMean <- kraken %>% 
  filter(Level.of.Classification == "P") %>% 
  group_by(Type, Name) %>% 
  summarise(Hits=mean(Hits.at.taxon))

krakenPhylumMeanTopFiltered <- as.data.frame(krakenPhylumMean) %>% 
  filter(Phylum %in% topPhyla)

names(krakenPhylumMean) <- str_replace(names(krakenPhylumMean), "Type", "Depth")

names(krakenPhylumMean) <- str_replace(names(krakenPhylumMean), "Name", "Phylum")

newCrayolaPalette6 <- c(
  "#CB7119", 
  "#404E5A",
  "#5F4F3A",
  "#D6AEDD",
  "#FEBAAD",
  "#6F7285",
  "#803790",
  "#0095B6",
  "#FFCD48",
  "#F653A6"
  )

newCrayolaPalette7 <- c(
  "#FE6F5E",
  "#346114",
  "#CB7119",
  "#D6AEDD",
  "#0095B6",
  "#F653A6",
  "#803790",
  "#FFCD48",
  "#631F41",
  "#404E5A"
  ) 

krakenPhylumMeanTop <- as.data.frame(krakenPhylumMean) %>% 
  group_by(Name) %>% 
  summarise(SumHits=sum(Hits)) %>% 
  arrange(desc(SumHits)) %>% 
  slice(1:5)

amrClassMeanTop <- as.data.frame(amrResultsClass) %>% 
  group_by(Class) %>% 
  summarise(SumHits=sum(Hits)) %>% 
  arrange(desc(SumHits)) %>% 
  slice(1:5)

topAMRClass <- as.vector(amrClassMeanTop$Class)

topPhyla <- as.vector(krakenPhylumMeanTop$Name)

krakenPhylumMeanTopFiltered <- as.data.frame(krakenPhylumMean) %>% 
  filter(Name %in% topPhyla)

amrClassTopFiltered <- as.data.frame(amrResultsClass) %>% 
  filter(Class %in% topAMRClass)

krakenPhylumMeanTopFiltered <- as.data.frame()

# Plots for poster

ggplot(amrClassTopFiltered, aes(x=Depth, y=Hits, fill=Class)) + 
  geom_bar(stat="identity", position="fill") + 
  scale_fill_manual(values=newCrayolaPalette9) +
  theme(
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    axis.text.x = element_text(size = 28),
    axis.text.y = element_text(size = 28),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 28),
    legend.key = element_rect(size = 2),
    legend.key.size = unit(2, "lines"),
    legend.spacing = unit(0.2,"lines")
  ) +
  labs(x='\nDepth', 
       y='Hits\n',
       fill = 'Class\n')

ggplot(krakenPhylumMeanTopFiltered, aes(x=Depth, y=Hits, fill=Phylum)) +
  geom_bar(stat="identity", position="fill") + 
  scale_fill_manual(values=newCrayolaPalette8) +
  theme(
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    axis.text.x = element_text(size = 28),
    axis.text.y = element_text(size = 28),
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 28),
    legend.key = element_rect(size = 2),
    legend.key.size = unit(2, "lines"),
    legend.spacing = unit(2,"lines")
  ) +
  labs(x='\nDepth', 
       y='Hits\n',
       fill = 'Phylum\n')
  #xlab('\nDepth') +
  #ylab('Hits\n')

ggplot(amrResultsSumName, aes(Sample_type, Hits)) + geom_bar(stat="identity")

# Rarefaction for list of datasets
amrResultsTidy <- amrResultsTidy <- amrResults %>% gather(Level, LevelName, c(1,6:8))

amrResultsTidy$Depth <- str_replace(amrResultsTidy$Sample, "\\d+_","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "\\d$","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "full", "D1")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "half", "D0.5")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "quar*", "D0.25")
amrResultsTidy$LevelName <- as.factor(amrResultsTidy$LevelName)

amrResultsList <- split(amrResultsTidy,amrResultsTidy$Level)

amrResultsSummary <- lapply(amrResultsList, function(x){
  summarizeAMRbySample(x)
})

amrResultsWide <- lapply(amrResultsSummary, function(x){
  widenAMR(x)
})

amrResultsMat <- lapply(amrResultsWide, function(x){
  matrixAMR(x)
})

# Awesome!!!

amrResultsMat2 <- lapply(amrResultsMat, function(x){
  
  t(x)
  
})

# One specific rarefaction curve

rarec(amrResultsMat2[['Name']], step=100, sample=min(rowSums(amrResultsMat2[['Name']])))

# Rarefaction curves for all members of the list

amrRarefy <- lapply(amrResultsMat2, function(x){
raremax <- min(rowSums(x))
rarecurve(x, step=5, sample=raremax)
})

# Unlisting rarefied data and isolating DFs
amrRarefyDF2 <- lapply(amrRarefy, unlist)

amrRarefyDF2 <- lapply(amrRarefyDF2, function(x){
  data.frame(otus=x,subsample=attr(x, "names"))
})

# Isolate Class DF
amrRarefyClass <- amrRarefyDF2[['Class']]

# Avoid the work below by making sure we have a named list

amrRarefyClassDF$Test <- ifelse(amrRarefyClassDF$subsample == "N1", 
                                amrRarefyClassDF$Test == sapply(sampleNames, function(x){x}), NA)


amr_results_class_wide <- amr_results_class %>%
    spread(Class, Hits, fill = 0)

amr_results_class_wide <- amr_results_class %>%
  spread(Sample, Hits, fill = 0)

table(annotations$class)

amr_results_class_wide <- amr_results_class %>% 
  spread(Class, Hits, fill = 0)

amr_results_class_mat <- amr_results_class_wide[,2:ncol(amr_results_class_wide)]
amr_class_norm <- data.frame(MRcounts(amr_class_experiment, norm=T))
amr_class_norm <- t(amr_class_norm)
row.names(amr_class_norm) <- str_replace(row.names(amr_class_norm), "X", "")


rarecurveROP <- function (x, step, sample, xlab = "Sample Size", ylab = "Species", 
          label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")?invisible
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)

  return(list(Nmax, Smax))
    
}


rarecurve_ROP <- function (x, step = 5, sample=raremax, label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  
  #plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
  #    type = "n", ...)
  #if (!missing(sample)) {
  #  abline(v = sample)
  #  rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
  #                                         y = z, xout = sample, rule = 1)$y)
  #  abline(h = rare, lwd = 0.5)
  #}
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    #lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  
  if (label) {
    samples <- as.character(row.names(x))
  }
  
  return(list(N=N, Out=out[[ln]], Samples=samples))
  #invisible(out)
}

rarefy_ggplot <- rarecurve_ROP(amr_results_class_mat)
