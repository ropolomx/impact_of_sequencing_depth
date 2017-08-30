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

# Get the paths of the files that we want to process

# Try to fix the code below if reading and concatenating all
# rarefaction files with R

# However, it can be done much more easily with Python-Pandas

amrRarefiedConcat <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/rarefiedConcat.csv')

amrFiltered <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/cov_sampler_parsed/amrFiltered_75_genefrac.csv')

# Remove seqtk data from dataframe

amrFiltered <- amrFiltered %>% filter(Sample_type != "D0.25_seqtk")

# Steps below are for further data cleaning. This might not be needed if the 
# dataset is already filtered.

# Need to fix the halves

fixSampleName <- amrRarefiedDF[amrRarefiedDF$Sample %in% c("H_006", "H_007", "H_008"),]

fixSampleName$Sample <- str_extract(fixSampleName$SampleName, "H_00[6-8]_\\w{4}")

amrRarefiedDF <- amrRarefiedDF %>% filter(!Sample %in% c("H_006", "H_007", "H_008"))

amrRarefied <- rbind(amrRarefiedDF, fixSampleName)

# Update path accordingly

amrRarefiedFiles <- Sys.glob(file.path("~",
                                      "amr",
                                      "2-4-8_results",
                                      "2_4_8_study_RZ",
                                      "Results_Aug2017_75_gene_frac",
                                      "*",
                                      "*rarefied*.tab*"))


# Split the path names and extract the sample name only

amrResultsNames <- str_split(amrRarefiedFiles, pattern = "\\/")

# Extraction of sample name is being done by extracting the 8th element out of each list element.

# TODO: Explore if this can be done with a print statement

amrResultsNames <- amrResultsNames %>% 
  map(function(x){
    filename <- x[9]
  })

# Split filename strings by dot character. Exclude extension.

amrResultsNames <- amrResultsNames %>%
  map(function(x){
   rarSampleName <- str_split(x, pattern="\\.")
   rarSampleName
  })

# Remove one level of hierarchy to list

amrResultsNames <- flatten(amrResultsNames)

# Extract only first element (sample name and AMR level)

amrResultsNames <- amrResultsNames %>%
  map(function(x){
    x <- x[1]
    x
  })

amrResultsNames <- unlist(amrResultsNames)

# Let's now read all the Coverage Sampler tabular files
# We are using readr (read_tsv)
# We are also using the list of sample names extracted in the previous function
# to set the names of the list elements
# This will make life so much easier!

amrResults <- amrRarefiedFiles %>%
  map(read_tsv, col_names = c("Sample", "Counts")) %>%
  set_names(nm=amrResultsNames)

# Need to add a column with sample name to the Coverage Sampler output
# Also need to add a column with sample depth

#TODO: Attempt to use map2 or pmap

amrArgs <- list(amrResults, amrResultsNames)

amrResults <- amrArgs %>%
  map2(function(x) {
   x$Sample <- rep(amrResultsNames, nrow(x))
   x
  })


# Join all the datasets into one dataframe that will be analyzed

amrResults <- do.call(amrResults, "rbind")

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

# Start here if AMR results have been filtered already

# Tidy the dataset

amrResultsTidy <- amrFiltered %>% 
  gather(Category, CategoryName, c(1:4))

amrResultsTidy$Category <- factor(amrResultsTidy$Category, 
                               levels = c('Class', 'Mechanism', 'Group', 'Gene'))

amrCategories <- levels(amrResultsTidy$Category)

#taxonLevels <- c("D", "K", "C", "O", "F", "G", "S", "-")

# Split tidy data frame of AMR results according to the Level vector

amrResultsList <- amrResultsTidy %>% 
  split(.$Category) %>% #splitting dataframe using base R but with purrr's piping
  set_names(nm=amrCategories)

# Summarize results: calculate mean by Depth

amrResultsSummary <- lapply(amrResultsList, function(x){
  summarizeAMRbyDepth(x)
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
  rarecurve(x, step=5, sample=raremax)
},mc.cores=3)

# Use microbenchmark to compare between the two approaches

#mbm <- microbenchmark(
#  mapping = map(amrResultsMat, function(x){
#    raremax <- min(rowSums(x))
#    rarecurve(x, step=5, sample=raremax)}),
#  multicore = mclapply(amrResultsMat, function(x){
#    raremax <- min(rowSums(x)) 
#    rarecurve(x, step=5, sample=raremax)}, mc.cores=12),
#  times=2
#)

# Results
# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# mapping 432.4910 432.4910 439.0898 439.0898 445.6887 445.6887     2
# multicore 140.0049 140.0049 140.8782 140.8782 141.7515 141.7515     2

# mclapply wins!

# Rename list of rarefied data using the sample names
# Need to think of a better way to extract the sample names

samples <- rownames(amrResultsMat[["Class"]])

# Review syntax here
# Might need to add function to utility function list

amrRarefy <- mclapply(amrRarefy, function(x) set_names(x,samples), mc.cores=3)

# Generate list of dataframes
# Hacky, but it works. Need to make it cleaner and faster.

amrRarefy <- amrRarefy %>% 
  map(as_vector)

amrRarefyDF <- map(amrRarefy, function(x) as_tibble(x, attr(x, "names")))

# Split rownames in order to generate columns with useful information

amrRarefyDF <- mclapply(amrRarefyDF, function(x) {
  x$Sample <- row.names(x)
  x$Subsample <- str_extract(x$Sample, "N\\d+")
  x$Subsample <- as.numeric(str_replace(x$Subsample, "N",""))
  x$Depth <- str_replace(x$Sample, "\\.N.*$", "")
  x$Depth <- str_replace(x$Depth, "^\\d+_", "")
  x$SampleNumber <- str_replace(x$Sample, "_.*", "")
  x$Sample <- str_replace(x$Sample, "\\.N.*$", "")
  x
}, 
mc.cores=3
)

# Generate one single dataframe and create column for AMR Level

amrRarefyDF <- do.call("rbind", amrRarefyDF)
amrRarefyDF$AMRLevel <- row.names(amrRarefyDF)
amrRarefyDF$AMRLevel <- str_extract(amrRarefyDF$AMRLevel, "^\\w+")

# Plot rarefaction curves using Rarefaction Analyzer data

# Slice dataset by AMR level

# Remove seqtk data (not needed for now)

amrRarefiedConcat <- amrRarefiedConcat %>%
  filter(Depth != "seqtk")

amrRarefiedConcat$Depth <- str_replace(amrRarefiedConcat$Depth, 'half[1-2]', 'half')

amrRarefiedMean <- amrRarefiedConcat %>%
  group_by(Depth,PercentSampling, Level) %>%
  summarise(MeanCounts = mean(Counts))

amrRarefiedClass <- amrRarefiedMean[amrRarefiedMean$Level == "class",]

amrRarefiedMechanism <- amrRarefiedConcat[amrRarefiedConcat$Level == "mechanism",]

amrRarefiedGroup <- amrRarefiedConcat[amrRarefiedConcat$Level == "group",] 
  
amrRarefiedGene <- amrRarefiedConcat[amrRarefiedConcat$Level == "gene",]


amrClassRarCurve <- amrRarefiedClass %>% 
  ggplot(aes(PercentSampling, Counts, color=Depth)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(5, 100, 5)) +
  facet_wrap( ~ SampleID, nrow = 2, ncol = 4, scales = "free")

amrClassMeanCurve <- amrRarefiedClass %>% 
  ggplot(aes(PercentSampling, MeanCounts, color=Depth)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(5, 100, 5))
  #facet_wrap( ~ SampleID, nrow = 2, ncol = 4, scales = "free")

amrMechanismRarCurve <- amrRarefiedMechanism %>% 
  ggplot(aes(PercentSampling, Counts, color=Depth)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(5, 100, 5)) +
  facet_wrap( ~ SampleID, nrow = 2, ncol = 4, scales = "free")

amrGroupRarCurve <- amrRarefiedGroup %>% 
  ggplot(aes(PercentSampling, Counts, color=Depth)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(5, 100, 5)) +
  facet_wrap( ~ SampleID, nrow = 2, ncol = 4, scales = "free")

amrGeneRarCurve <- amrRarefiedGene %>% 
  ggplot(aes(PercentSampling, Counts, color=Depth)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(5, 100, 5)) +
  facet_wrap( ~ SampleID, nrow = 2, ncol = 4, scales = "free")
