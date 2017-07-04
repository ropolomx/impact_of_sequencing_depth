# Rarefaction analysis of 2-4-8 study data
# Processing of AMR and Kraken data

# Load necessary packages

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(metagenomeSeq)
library(vegan)

# Source script of utility functions

source('rarefaction_utility_functions.R')

# Read AMR and Kraken data

# TODO: Need to add a column with sample name to the Coverage Sampler output
#amrResults <- read_tsv('FC_N013.tabular_parsed.tab')

amrResults <- read_csv('noelle/AMR/amr_new_dataframe_ROP.csv')

krakenResults <- read_csv('noelle/Kraken/kraken_new_dataframe.csv')

# TODO: If output is generated from Rarefaction Analyzer,
# see how it can be applied to the analysis of Kraken data, 
# although it seems that program can only be applied to process SAM files

# Convert AMR results to "tidy" format
# Kraken table is already in tidy format.

amrResultsTidy <- amrResults %>% 
  gather(Level, LevelName, c(1,6:8))

# Change the names of the AMR sequencing depths columns

amrResultsTidy$Depth <- str_replace(amrResultsTidy$Sample, "\\d+_","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "\\d$","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "full", "D1")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "half", "D0.5")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "quar*", "D0.25")

# Change the names of the Kraken sequencing depths columns

krakenResults$Depth <- str_replace(krakenResults$Sample, "_filtered2_report$","")
krakenResults$Depth <- str_replace(krakenResults$Depth, "^\\d+_","")
krakenResults$Depth <- str_replace(krakenResults$Depth, "\\d$","")
krakenResults$Depth <- str_replace(krakenResults$Depth, "full", "D1")
krakenResults$Depth <- str_replace(krakenResults$Depth, "half", "D0.5")
krakenResults$Depth <- str_replace(krakenResults$Depth, "quar*", "D0.25")

# Convert the LevelName to Factor

amrResultsTidy$LevelName <- as.factor(amrResultsTidy$LevelName)

# Keep results with over 80 % coverage ratio

amrResultsTidy <- amrResultsTidy %>% 
  filter(Coverage_Ratio >= 0.80)

# Vectors containing amr Levels and taxon levels to analyze
#amrLevels <- c("Class", "Mechanism", "Group", "Gene Id")

amrResultsTidy$Level <- factor(amrResultsTidy$Level, 
                               levels = c('Class', 'Mechanism', 'Group', 'Name'))

amrLevels <- levels(amrResultsTidy$Level)

taxonLevels <- c("D", "K", "C", "O", "F", "G", "S", "-")


# Split tidy data frame of AMR results according to the Level vector

amrResultsList <- amrResultsTidy %>% 
  split(.$Level) %>% #splitting dataframe using base R but with purrr's piping
  set_names(nm=amrLevels)

# Summarize results: add counts by Sample and AMR Level

amrResultsSummary <- lapply(amrResultsList, function(x){
  summarizeAMRbySample(x)
})

# Convert AMR results to wide format

amrResultsWide <- lapply(amrResultsSummary, function(x){
  widenAMR(x)
})

# Convert AMR results from wide format to matrix

amrResultsMat <- lapply(amrResultsWide, function(x){
  matrixAMR(x)
})


# Transpose matrices
# Awesome!!!

amrResultsMat <- lapply(amrResultsMat, function(x){
  
  t(x)
  
})

# Example: one specific rarefaction curve

#rarec(amrResultsMat[['Name']], 
#      step=100, 
#      sample=min(rowSums(amrResultsMat[['Name']])))


# Rarefaction curves for all members of the list
# This step requires optimization (parallel mapping or compiling, maybe?)

# purrr version

#amrRarefy <- map(amrResultsMat, function(x){
#  raremax <- min(rowSums(x))
#  rarecurve_ROP(x, step=5, sample=raremax)
#})

# mclapply version

amrRarefy <- mclapply(amrResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve_ROP(x, step=5, sample=raremax)
},mc.cores=3)

# Use microbenchmark to compare between the two approaches

mbm <- microbenchmark(
  mapping = map(amrResultsMat, function(x){
    raremax <- min(rowSums(x))
    rarecurve_ROP(x, step=5, sample=raremax)}),
  multicore = mclapply(amrResultsMat, function(x){
    raremax <- min(rowSums(x)) 
    rarecurve_ROP(x, step=5, sample=raremax)}, mc.cores=12),
  times=2
)

# Results
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# mapping 432.4910 432.4910 439.0898 439.0898 445.6887 445.6887     2
# multicore 140.0049 140.0049 140.8782 140.8782 141.7515 141.7515     2

# mclapply wins!

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
