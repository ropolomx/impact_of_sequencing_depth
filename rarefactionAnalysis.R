# Rarefaction analysis of 2-4-8 study data
# Analysis of alpha-diversity and species richness

# Load necessary packages

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)
library(metagenomeSeq)
library(vegan)
library(scales)

# Source script with utility functions

# Will create an R package in the future

source('rarefaction_utility_functions.R')

# Read AMR and Kraken data

# Another alternative
# Results generated with Coverage Sampler and filtered with Python-Pandas
# Filtering involved keeping results with gene fraction >= 75% and 
# removing all those results with genes that require SNP confirmation.

amrResultsFiltered <- read_csv('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/cov_sampler_parsed/amrFiltered_75_genefrac.csv')

amrReadstoHitRatio <- read_tsv('~/amr/2-4-8_results/2_4_8_study_RZ/hitToReadRatios.tsv')

# Read Kraken concatenated and filtered file (no Eukaryotes, and no PhiX)

krakenResultsFiltered <- read.table('~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/allKraken_FHQ/kraken_filtered/krakenConcat.tsv', sep="\t", header=TRUE)


# Remove D0.25_seqtk data from filtered data frame

amrResultsFiltered <- amrResultsFiltered %>% 
  filter(Sample_type != "D0.25_seqtk")

# Replace "Name" with "Gene" in headers of the filtered data frame

names(amrResultsFiltered) <- str_replace(names(amrResultsFiltered), "Name", "Gene")

# Convert AMR filtered results to "tidy" format in order to have a column with 
# all AMR levels

amrResultsTidy <- amrResultsFiltered %>% 
  gather(Level, LevelName, c(1,6:8))

# Change the names of the AMR sequencing depths columns

amrResultsTidy$Depth <- str_replace(amrResultsTidy$Sample, "\\d+_","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "\\d$","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "full", "D1")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "half", "D0.5")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "quar*", "D0.25")

#### Starting with Filtered data ####

# Start here if AMR results have been filtered already

# Tidy the dataset

amrResultsTidy <- amrResultsFiltered %>% 
  gather(Category, CategoryName, c(1:4))

amrResultsTidy$Category <- factor(amrResultsTidy$Category, 
                                  levels = c('Class', 'Mechanism', 'Group', 'Gene'))

amrCategories <- levels(amrResultsTidy$Category)

krakenTaxa <- levels(krakenResultsFiltered$TaxRank)

# Split tidy data frame of AMR results according to the Level vector

amrResultsList <- amrResultsTidy %>% 
  split(.$Category) %>% #splitting dataframe using base R but with purrr's piping
  set_names(nm=amrCategories)

krakenResultsList <- krakenResultsFiltered %>% 
  split(.$TaxRank) %>% #splitting dataframe using base R but with purrr's piping
  set_names(nm=krakenTaxa)
# Summarize results: calculate mean by Depth

amrResultsSummary <- lapply(amrResultsList, function(x){
  summarizeAMRbyDepth(x)
})

krakenResultsSummary <- lapply(krakenResultsList, function(x){
  summarizeKrakenbyDepth(x)
})

# Convert AMR results to wide format

amrResultsWide <- lapply(amrResultsSummary, function(x){
  widenAMR(x)
})


krakenResultsWide <- lapply(krakenResultsSummary, function(x){
  widenKraken(x)
})

# Convert AMR results from wide format to matrix

amrResultsMat <- lapply(amrResultsWide, function(x){
  matrixAMR(x)
})


krakenResultsMat <- lapply(krakenResultsWide, function(x){
  matrixKraken(x)
})

# Transpose matrices
# Awesome!!!

amrResultsMat <- lapply(amrResultsMat, function(x){
  
  t(x)
  
})

krakenResultsMat <- lapply(krakenResultsMat, function(x){
  
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

amrRarCurve <- mclapply(amrResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=50, sample=raremax)
},mc.cores=10)

krakenRarCurve <- mclapply(krakenResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=1000, sample=raremax)
},mc.cores=10)

# Use mclapply function
krakenAlphaRarefaction <- mclapply(krakenResultsMat, function(x){
  alpha_rarefaction(x, minlevel=0)
},mc.cores = 10)

amrAlphaRarefaction <- mclapply(amrResultsMat, function(x){
  alpha_rarefaction(x, minlevel=0)
}, mc.cores=10)

krakenAlphaRarefaction <- lapply(krakenAlphaRarefaction2, function(x) data.table(
  ID=names(x$raw_species_abundance), 
  RawSpeciesAbundance=as.numeric(x$raw_species_abundance), 
  RarSpeciesAbundance=as.numeric(x$rarefied_species_abundance), 
  AlphaDiv=as.numeric(x$alphadiv), 
  Shannon=as.numeric(x$shannon), 
  Evenness=as.numeric(x$evenness)
))

amrAlphaRarefaction <- lapply(amrAlphaRarefaction, function(x) data.table(
  ID=names(x$raw_species_abundance), 
  RawSpeciesAbundance=as.numeric(x$raw_species_abundance), 
  RarSpeciesAbundance=as.numeric(x$rarefied_species_abundance), 
  AlphaDiv=as.numeric(x$alphadiv), 
  Shannon=as.numeric(x$shannon), 
  Evenness=as.numeric(x$evenness)
))

amrCategories <- levels(amrResultsTidy$Category)

amrAlphaRarefaction[[1]]$Level <- rep(amrCategories[[1]], length(amrAlphaRarefaction[[1]]$ID))
amrAlphaRarefaction[[2]]$Level <- rep(amrCategories[[2]], length(amrAlphaRarefaction[[2]]$ID))
amrAlphaRarefaction[[3]]$Level <- rep(amrCategories[[3]], length(amrAlphaRarefaction[[3]]$ID))
amrAlphaRarefaction[[4]]$Level <- rep(amrCategories[[4]], length(amrAlphaRarefaction[[4]]$ID))

krakenAlphaRarefaction[[1]]$Level <- rep(krakenTaxa[[1]], length(krakenAlphaRarefaction[[1]]$ID))
krakenAlphaRarefaction[[2]]$Level <- rep(krakenTaxa[[2]], length(krakenAlphaRarefaction[[2]]$ID))
krakenAlphaRarefaction[[3]]$Level <- rep(krakenTaxa[[3]], length(krakenAlphaRarefaction[[3]]$ID))
krakenAlphaRarefaction[[4]]$Level <- rep(krakenTaxa[[4]], length(krakenAlphaRarefaction[[4]]$ID))
krakenAlphaRarefaction[[5]]$Level <- rep(krakenTaxa[[5]], length(krakenAlphaRarefaction[[5]]$ID))
krakenAlphaRarefaction[[6]]$Level <- rep(krakenTaxa[[6]], length(krakenAlphaRarefaction[[6]]$ID))
krakenAlphaRarefaction[[7]]$Level <- rep(krakenTaxa[[7]], length(krakenAlphaRarefaction[[7]]$ID))
krakenAlphaRarefaction[[8]]$Level <- rep(krakenTaxa[[8]], length(krakenAlphaRarefaction[[8]]$ID))

krakenAlphaRarefaction2[[1]]$Level <- rep(krakenTaxa[[1]], length(krakenAlphaRarefaction2[[1]]$ID))
krakenAlphaRarefaction2[[2]]$Level <- rep(krakenTaxa[[2]], length(krakenAlphaRarefaction2[[2]]$ID))
krakenAlphaRarefaction2[[3]]$Level <- rep(krakenTaxa[[3]], length(krakenAlphaRarefaction2[[3]]$ID))
krakenAlphaRarefaction2[[4]]$Level <- rep(krakenTaxa[[4]], length(krakenAlphaRarefaction2[[4]]$ID))
krakenAlphaRarefaction2[[5]]$Level <- rep(krakenTaxa[[5]], length(krakenAlphaRarefaction2[[5]]$ID))
krakenAlphaRarefaction2[[6]]$Level <- rep(krakenTaxa[[6]], length(krakenAlphaRarefaction2[[6]]$ID))
krakenAlphaRarefaction2[[7]]$Level <- rep(krakenTaxa[[7]], length(krakenAlphaRarefaction2[[7]]$ID))
krakenAlphaRarefaction2[[8]]$Level <- rep(krakenTaxa[[8]], length(krakenAlphaRarefaction2[[8]]$ID))

amrAlphaRarefactionDF <- lapply(amrAlphaRarefaction, function(x){
  x <- as.data.frame(x)
  x
})

krakenAlphaRarefactionDF <- lapply(krakenAlphaRarefaction, function(x){
  x <- as.data.frame(x)
  x
})

krakenAlphaRarefaction2DF <- lapply(krakenAlphaRarefaction2, function(x){
  x <- as.data.frame(x)
  x
})

amrAlphaRarefactionDF <- do.call("rbind", amrAlphaRarefactionDF)
  
krakenAlphaRarefactionDF <- do.call("rbind", krakenAlphaRarefactionDF)

krakenAlphaRarefaction2DF <- do.call("rbind", krakenAlphaRarefaction2DF)

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

amrSamples <- rownames(amrResultsMat[["Class"]])

krakenSamples <- rownames(krakenResultsMat[["D"]])

# Review syntax here
# Might need to add function to utility function list

amrRarCurve <- mclapply(amrRarCurve, function(x){
  set_names(x,amrSamples)}, mc.cores=3)

krakenRarCurve <- mclapply(krakenRarCurve, function(x){
  set_names(x,krakenSamples)}, mc.cores=10)
# Generate list of dataframes
# Hacky, but it works. Need to make it cleaner and faster.

amrRarCurve <- amrRarCurve %>% 
  map(as_vector)

krakenRarCurve <- krakenRarCurve %>% 
  map(as_vector)

krakenAlphaDiversity <- krakenAlphaDiversity %>%
  map(as_vector)

amrRarCurveDF <- map(amrRarCurve, function(x) as_tibble(x, attr(x, "names")))

krakenRarefyDF <- map(krakenRarCurve, function(x) as_tibble(x, attr(x, "names")))

krakenAlphaDivDF <- map(krakenAlphaDiversity, function(x) as_tibble(x, attr(x, "names")))

# Split rownames in order to generate columns with useful information

amrRarCurveDF <- mclapply(amrRarCurveDF, function(x) {
  x$Sample <- row.names(x)
  x$Subsample <- str_extract(x$Sample, "N\\d+")
  x$Subsample <- as.numeric(str_replace(x$Subsample, "N",""))
  x$SampleID <- str_replace(x$Sample, "\\.N.*$", "")
  x$Depth <- str_replace(x$Sample, "_.*", "")
  x
}, 
mc.cores=3
)

krakenRarefyDF <- mclapply(krakenRarefyDF, function(x) {
  x$Sample <- row.names(x)
  x$Subsample <- str_extract(x$Sample, "\\.N\\d+")
  x$Subsample <- as.numeric(str_replace(x$Subsample, "\\.N",""))
  x$Depth <- str_replace(x$Sample, "_.*", "")
  x$SampleID <- str_extract(x$Sample, "_.*\\.")
  x$SampleID <- str_replace(x$SampleID, "_", "")
  x$SampleID <- str_replace(x$SampleID, "\\.", "")
  x$Sample <- str_replace(x$Sample, "\\.N.*$", "")
  x
}, 
mc.cores=10
)

amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$ID, "_.*", "")

krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$ID, "_.*", "")


krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$ID, "_.*", "")
# Generate one single dataframe and create column for AMR Level

amrRarCurveDF <- do.call("rbind", amrRarCurveDF)
amrRarCurveDF$AMRLevel <- row.names(amrRarCurveDF)
amrRarCurveDF$AMRLevel <- str_extract(amrRarCurveDF$AMRLevel, "^\\w+")
amrRarCurveDF$Depth <- str_replace(amrRarCurveDF$Depth, "F", "D1")
amrRarCurveDF$Depth <- str_replace(amrRarCurveDF$Depth, "H", "D0.5")
amrRarCurveDF$Depth <- str_replace(amrRarCurveDF$Depth, "QD", "D0.25")

krakenRarefyDF <- do.call("rbind", krakenRarefyDF)
krakenRarefyDF$krakenLevel <- row.names(krakenRarefyDF)
krakenRarefyDF$krakenLevel <- str_extract(krakenRarefyDF$krakenLevel, "^.\\.")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "\\.", "")
krakenRarefyDF$Depth <- str_replace(krakenRarefyDF$Depth, "F", "D1")
krakenRarefyDF$Depth <- str_replace(krakenRarefyDF$Depth, "H", "D0.5")
krakenRarefyDF$Depth <- str_replace(krakenRarefyDF$Depth, "QD", "D0.25")


krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "F", "D1")
krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "H", "D0.5")
krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "QD", "D0.25")

amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "F", "D1")
amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "H", "D0.5")
amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "QD", "D0.25")

krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "F", "D1")
krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "H", "D0.5")
krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "QD", "D0.25")


# Plot rarefaction curves using Rarefaction Analyzer data

# Slice dataset by AMR level

# Remove seqtk datacolor=Sample_type (not needed for now)

amrRarefiedConcat <- amrRarefiedConcat %>%
  filter(Depth != "seqtk")

amrRarefiedConcat$Depth <- str_replace(amrRarefiedConcat$Depth, 'half[1-2]', 'half')

amrRarefiedMean <- amrRarefiedConcat %>%
  group_by(Depth,PercentSampling, Level) %>%
  summarise(MeanCounts = mean(Counts))

amrRarefiedClass <- amrRarefiedConcat[amrRarefiedConcat$Level == "class",]

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

# Plot aggregated rarefaction curves

krakenRarPhylum <- krakenRarefiedMean[krakenRarefiedMean$krakenLevel == "P",]
krakenRarClass <- krakenRarefiedMean[krakenRarefiedMean$krakenLevel == "C",]
krakenRarOrder <- krakenRarefiedMean[krakenRarefiedMean$krakenLevel == "O",]
krakenRarFamily <- krakenRarefiedMean[krakenRarefiedMean$krakenLevel == "F",]
krakenRarGenus <- krakenRarefiedMean[krakenRarefiedMean$krakenLevel == "G",]
krakenRarSpecies <- krakenRarefiedMean[krakenRarefiedMean$krakenLevel == "S",]

krakenAllPhylum <- krakenRarefyDF[krakenRarefyDF$krakenLevel == "P",]
krakenAllClass <- krakenRarefyDF[krakenRarefyDF$krakenLevel == "C",]
krakenAllOrder <- krakenRarefyDF[krakenRarefyDF$krakenLevel == "O",]
krakenAllFamily <- krakenRarefyDF[krakenRarefyDF$krakenLevel == "F",]
krakenAllGenus <- krakenRarefyDF[krakenRarefyDF$krakenLevel == "G",]
krakenAllSpecies <- krakenRarefyDF[krakenRarefyDF$krakenLevel == "S",]

# Same colour palette as the matplotlib_venn default color palette

vennPalette <- c("#F8766D", "#00BA38", "#619CFF")
vennPalette <- rev(vennPalette)

# Write a function!!!

# Make all of this functional programming!!!

krakenAllPhylumRarCurve <- krakenAllPhylum %>%
  ggplot(aes(Subsample, value, color=Depth)) +
  geom_point(alpha=0.6) + 
  xlab = "Number of Phyla" +
  scale_color_manual(vennPalette) +
  facet_grid(. ~ Depth) +
  tras_coord(x="log10")
  
  
krakenAllClassRarCurve <- krakenAllClass %>%
  ggplot(aes(Number_of_Reads, Counts, color=Depth)) +
  geom_point(alpha=0.6) +
  scale_fill_manual(c(vennPalette))

krakenAllOrderRarCurve <- krakenAllOrder %>%
  ggplot(aes(Subsample, value, color=Depth)) +
  geom_point(alpha=0.7) +

krakenAllFamilyRarCurve <- krakenAllFamily %>%
  ggplot(aes(Subsample, value, color=Depth)) +
  geom_point(alpha=0.7) +

krakenAllGenusRarCurve <- krakenAllGenus %>%
  ggplot(aes(Subsample, value, color=Depth)) +
  geom_point(alpha=0.7) +

krakenAllSpeciesRarCurve <- krakenAllSpecies %>%
  ggplot(aes(Subsample, value, color=Depth)) +
  geom_point(alpha=0.7) +

# Attempt to more functional programming
 
krakenRarefyDF <- krakenRarefyDF %>% filter(!krakenLevel %in% c("-", "D"))

krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "P", "Phyla")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "C", "Classes")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "O", "Orders")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "F", "Families")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "G", "Genera")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "S", "Species")

krakenAllRarCurveList <- krakenRarefyDF %>% 
  split(.$krakenLevel)

amrAllList <- amrRarCurveDF %>% 
  split(.$AMRLevel)

# Alpha Diversity

krakenAlphaRarefactionDF <- krakenAlphaRarefactionDF %>% filter(!Level %in% c("-", "D"))
krakenAlphaRarefaction2DF <- krakenAlphaRarefaction2DF %>% filter(!Level %in% c("-", "D"))

krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "P", "Phyla")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "C", "Classes")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "O", "Orders")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "F", "Families")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "G", "Genera")
krakenAlphaRarefactionDF$Level <- str_replace(krakenAlphaRarefactionDF$Level, "S", "Species")

krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "P", "Phyla")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "C", "Classes")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "O", "Orders")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "F", "Families")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "G", "Genera")
krakenAlphaRarefaction2DF$Level <- str_replace(krakenAlphaRarefaction2DF$Level, "S", "Species")

krakenAlphaRarefactionDF$Level <- factor(krakenAlphaRarefactionDF$Level, 
                                  levels = c('Phyla', 'Classes', 'Orders', 'Families', 'Genera', 'Species'))


amrRarCurveDF$AMRLevel <- factor(amrRarCurveDF$AMRLevel, 
                                  levels = c('Class', 'Mechanism', 'Group', 'Gene'))

amrAlphaRarefactionDF$Level <- factor(amrAlphaRarefactionDF$Level, 
                                 levels = c('Class', 'Mechanism', 'Group', 'Gene'))


krakenAllAlphaDivList <- krakenAlphaDivDF %>% 
  split(.$krakenLevel)

amrAllRarCurves <- amrAllList %>%
  map(function(x){
    amrRarefactionCurve(x)
  })

krakenAllRarCurvesLarger <- krakenAllRarCurveList %>%
  map(function(x){
    krakenRarefactionCurve(x)
  })


amrAllAlphaBoxPlots <- amrAlphaRarefactionDF %>%
  amrAlphaDiv()

ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity/amrAlphaDiversityCB.png', 
       width = 10.50,
       height = 8.50,
       units = "in")

amrAllSpRawBoxPlots <- amrAlphaRarefactionDF %>%
    amrRawSpeciesRich()

ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity/amrSpeciesRichnessCB.png',
       plot = amrAllSpRawBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenAllAlphaBoxPlots <- krakenAlphaRarefactionDF %>%
    krakenAlphaDiv()

ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/')

krakenAllAlphaBoxPlots2 <- krakenAlphaRarefaction2DF %>%
    krakenAlphaDiv()

krakenAllSpRawBoxPlots <- krakenAlphaRarefaction2DF %>%
  krakenRawSpeciesRich()



# AMR rarefaction curves

# TODO: change code to ggsave

png(filename = "amrClassRarefaction2.png", width=1962, height = 1297)
print(amrAllRarCurves[[1]])
dev.off()

png(filename = "amrGeneRarefaction2.png", width=1962, height = 1297)
print(amrAllRarCurves[[2]])
dev.off()


png(filename = "amrGroupRarefaction2.png", width=1962, height = 1297)
print(amrAllRarCurves[[3]])
dev.off()

png(filename = "amrMechanismRarefaction2.png", width=1962, height = 1297)
print(amrAllRarCurves[[4]])
dev.off()

# Correlation plot: amr reads to hits

amrReadstoHitRatio$Sample_type <- str_replace(amrReadstoHitRatio$Sample_type, "F", "D1")

amrReadstoHitRatio$Sample_type <- str_replace(amrReadstoHitRatio$Sample_type, "H", "D0.5")

amrReadstoHitRatio$Sample_type <- str_replace(amrReadstoHitRatio$Sample_type, "Q", "D0.25")

amrReadsvsHits <- cor(x = amrReadstoHitRatio$Number_of_reads, amrReadstoHitRatio$AMR_hits, method = "spearman")

amrReadsvsHitsCorTest <- cor.test(x = amrReadstoHitRatio$Number_of_reads, amrReadstoHitRatio$AMR_hits, method = "spearman")

amrReadsvsHitsCor <- ggplot(amrReadstoHitRatio, aes(Number_of_reads, AMR_hits)) + 
  geom_point(aes(fill=Sample_type), 
             alpha=0.6, 
             size=10, 
             pch=21, 
             color="grey")+ 
  geom_smooth(aes(group=1, weight=0.2), 
              method="lm", 
              se=FALSE, 
              colour="grey", 
              alpha=0.5) 

amrReadsvsHitsCor +
  ylab("Number of AMR Hits\n") +
  xlab("\nNumber of reads") +
  theme(axis.text.y=element_text(size=35),
        axis.title.y=element_text(size=44),
        axis.text.x=element_text(size=35),
        axis.title.x=element_text(size=44),
        legend.title=element_text(size=36),
        legend.text=element_text(size=36, vjust=0.5),
        legend.key = element_rect(size = 2),
        legend.key.size = unit(2, "lines"),
        legend.spacing = unit(0.2,"lines"),
        panel.background = element_rect(fill = "grey90", colour = "grey80")) +
  scale_fill_manual(values=rev(cbPalette),
                     name="Sequencing Depth\n")
ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/amrReadsvsHitsCorCBScheme.png', 
       width=14, 
       height=8.50,
       units="in")

amrReadsvsHitsCor <- ggplot(amrReadstoHitRatio, aes(Number_of_reads, AMR_hits)) + 
  geom_point(alpha=0.4, size=4) + 
  geom_smooth(aes(group=1, weight=0.2), method="lm", se=FALSE, colour="grey", alpha=0.5)

amrReadsvsHitsCor + 
  ylab("Number of AMR Hits\n") + 
  xlab("\nNumber of reads") + 
  theme(axis.text.y=element_text(size=35),
        axis.title.y=element_text(size=44),
          axis.text.x=element_text(size=35),
          axis.title.x=element_text(size=44),
          legend.title=element_text(size=36),
          legend.text=element_text(size=36, vjust=0.5),
        legend.key = element_rect(size = 2),
        legend.key.size = unit(2, "lines"),
        legend.spacing = unit(0.2,"lines")) + 
  scale_color_manual(values=vennPalette, 
                    name="Sample type\n")

# Kruskal-Wallis tests

# Resistome

amrAlphaRarefactionLevels <- amrAlphaRarefactionDF %>%
  split(.$Level)

amrRichnessKruskal <- amrAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(RawSpeciesAbundance ~ Depth, data = x)
    })

amrRichnessPosthoc <- amrAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(RawSpeciesAbundance ~ Depth, data=x, dist="Chisq")
  })

amrAlphaKruskal <- amrAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(AlphaDiv ~ Depth, data = x)
  })

amrAlphaPosthoc <- amrAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(AlphaDiv ~ Depth, data = x, dist="Chisq")
  })

# Microbiome

krakenAlphaRarefaction2DF$Depth <- as.factor(krakenAlphaRarefaction2DF$Depth)

krakenAlphaRarefactionLevels <- krakenAlphaRarefaction2DF %>%
  split(.$Level)

krakenRichnessKruskal <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(RawSpeccolor=Sample_typeiesAbundance ~ Depth, data = x)
    })

krakenRichnessPosthoc <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(RawSpeciesAbundance ~ Depth, data=x, dist="Chisq")
  })

krakenAlphaKruskal <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    kruskal.test(AlphaDiv ~ Depth, data = x)
  })

krakenAlphaPosthoc <- krakenAlphaRarefactionLevels %>%
  map(function(x){
    posthoc.kruskal.nemenyi.test(AlphaDiv ~ Depth, data = x, dist="Chisq")
  })

# Generating kraken correlation plot

krakenPhylumResults <- krakenResultsList[["P"]]

krakenPhylumResults$SampleID <- str_replace(krakenPhylumResults$Sample, "_00[6-8]")

krakenPhylumSums <- krakenPhylumResults %>%
  group_by(TaxID, SampleID, Sample_Type) %>%
  summarise(MeanHits = mean(CladeReads)) %>%
  group_by(SampleID) %>%
  summarise(SumHits = sum(MeanHits))


cbPalette <- c("#FFCD48", colorschemes$Categorical.12[8], colorschemes$Categorical.12[12])

krakenReadsvsHitsCorTest <- cor.test(x = krakenReadstoHitRatio$Reads, krakenReadstoHitRatio$KrakenHits, method = "spearman")

krakenReadsvsHitsCor <- ggplot(krakenReadstoHitRatio, aes(Reads, KrakenHits)) + 
  geom_point(aes(fill=Sample_type), alpha=0.6, size=10, pch=21, color="grey")+ 
  geom_smooth(aes(group=1, weight=0.2), method="lm", se=FALSE, colour="grey", alpha=0.5) 

krakenReadsvsHitsCor +
  ylab("Number of Kraken Hits\n") +
  xlab("\nNumber of reads") +
  theme(axis.text.y=element_text(size=35),
        axis.title.y=element_text(size=44),
        axis.text.x=element_text(size=35),
        axis.title.x=element_text(size=44),
        legend.title=element_text(size=36),
        legend.text=element_text(size=36, vjust=0.5),
        legend.key = element_rect(size = 2),
        legend.key.size = unit(2, "lines"),
        legend.spacing = unit(0.2,"lines"),
        panel.background = element_rect(fill = "grey90", colour = "grey80")) +
  scale_fill_manual(values=rev(cbPalette),
                     name="Sequencing Depth\n")
ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/krakenReadsvsHitsCorCBScheme.png', width=14, height=8.50)
