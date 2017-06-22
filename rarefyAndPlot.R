library(dplyr)
library(tidyr)
library(metagenomeSeq)
library(vegan)
library(stringr)

source('rarefaction_functions.R')

amrResults <- read.csv('AMR/amr_new_dataframe_ROP.csv')

amrResultsTidy$Sample_type <- str_replace(amrResultsTidy$Sample, 
                                          "\\d+_","")
amrResultsTidy$Sample_type <- str_replace(amrResultsTidy$Sample_type, 
                                          "\\d$","")
amrResultsTidy$Sample_type <- str_replace(amrResultsTidy$Sample_type, 
                                          "full", "D1")
amrResultsTidy$Sample_type <- str_replace(amrResultsTidy$Sample_type, 
                                          "half", "D0.5")
amrResultsTidy$Sample_type <- str_replace(amrResultsTidy$Sample_type,
                                          "quar", "D0.25")
amrResultsTidy <- amrResultsTidy %>% filter(Coverage_Ratio >= 0.80)

amrLevels <- c("Class", "Mechanism", "Group")

amrSummaries <- lapply(amrLevels, function(i){summarizeAMRlevels(amrResults, i)})
amrSummariesWide <- lapply(amrSummaries, function(i){widenAMR(i)})

amrSummariesMat <- lapply(amrSummariesWide, function(i){matrixAMR(i)})

amrExperiment <- lapply(amrSummariesMat, function(i){normalizeAMR(i)})

amrRarefaction <- lapply(amrExperiment, function(X){alphaRarefactionROP(X)})


#alpha_rarefaction <- function(X, step=0.05, method='invsimpson') {
#    S <- specnumber(X, MARGIN=2)
#    raremax <- min(colSums(X))
#    tot <- colSums(X)
#    nr <- ncol(X)
#    out <- lapply(seq_len(nr), function(i){
#      n <- seq(1, tot[i], by = tot[i]*step)
#      if (n[length(n)] != tot[i])
#        n <- c(n, tot[i])
#      drop(rarefy(X[i, ], n))
#      })
#    Srare <- rarefy(X, raremax, MARGIN=2)
#    Xrare <- t(rrarefy(t(X), raremax))
#    alphadiv <- diversity(Xrare, index=method, MARGIN=2)
#    return(list(raw_species_abundance=S,
#                #rarefied_species_abundance=Srare,
#                #rarefied_data=Xrare,
#                alphadiv=alphadiv,
#                out=out))
#}

sampleNames <- lapply(amrExperiment, function(i){attr(i$raw_species_abundance, "names")})
unlistExperiments <- lapply(amrRarefaction$rarefy_out, unlist)
experimentDF <- lapply(unlistExperiments, data.frame)
rarefactionDF <- do.call("rbind",ggplot_rarefy_df[1:32])
names(rarefactionDF) <- "Classes"
rarefactionDF$Samples <- rep(Sample, each=21)
rarefactionDF$Sampling_size <- str_replace(row.names(rarefactionDF), "N","")
rarefactionDF$Sampling_size <- as.numeric(as.character(rarefactionDF$Sampling_size))

rarefactionDF$Sample_type <- str_replace(rarefactionDF$Samples, 
                                         "\\d+_", "")

rarefactionDF$Sample_type <- str_replace(rarefactionDF$Sample_type, 
                                         "\\d$","")

rarefactionDF$Sample_number <- str_replace(rarefactionDF$Samples, 
                                           "_.*", "")

halves <- rarefactionDF %>% filter(grepl("half",Sample_type))

rarefactionDFGroup <- rarefactionDF %>% 
  group_by(Sample_number, Sample_type) %>% 
  summarise(meanOTUs = mean(Classes), meanSamples=mean(Sampling_size))


# ggplot code for rarefaction plots

ggplot(rarefactionDFGroup, aes(meanSamples, meanOTUs)) + 
  geom_line() + 
  facet_grid(. ~ Sample_type)

ggplot(rarefactionDF, aes(Sampling_size,Classes)) + 
  geom_line(aes(color=Samples)) + 
  facet_grid(Sample_type~ Sample_number)

