# Rarefaction analysis of 2-4-8 study data
# Analysis of alpha-diversity and species richness

# Load packages -----------------------------------------------------------

library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(foreach)
library(ggplot2)
library(metagenomeSeq)
library(vegan)
library(scales)

# Source utility functions ------------------------------------------------

source('rarefaction_utility_functions.R')


# Colour palettes ---------------------------------------------------------

# Same colour palette as the matplotlib_venn default color palette
# Consider adding them to the utility functions script

vennPalette <- c("#F8766D", "#00BA38", "#619CFF")
vennPalette <- rev(vennPalette)

# Colourblind friendly palette

cbPalette <- c("#FFCD48", # Mango from Crayola palette
               dichromat::colorschemes$Categorical.12[8], # Blue from dichromat package
                 dichromat::colorschemes$Categorical.12[12]) # Red from dichromat package

# Load and filter AMR and Kraken data -------------------------------------

# Results generated with Coverage Sampler and filtered with Python-Pandas
# Filtering involved keeping results with gene fraction >= 75% and 
# removing all those results with genes that require SNP confirmation.

amrResultsFiltered <- read_csv(Sys.glob(file.path(
  '~',
  'aafc',
  'amr',
  'amrplusplus_rarefaction_analysis',
  '2_4_8_study_RZ',
  'Results_Aug2017',
  'Parsed_Aug2017',
  'amrFiltered_75_genefrac.csv'))
  )

amrReadstoHitRatio <- read_tsv('~/amr/2-4-8_results/2_4_8_study_RZ/hitToReadRatios.tsv')

# Read Kraken concatenated and filtered file (no Eukaryotes, and no PhiX)

krakenResultsFiltered <- read.delim(Sys.glob(file.path(
  '~',
  'aafc',
  'amr',
  'amrplusplus_rarefaction_analysis',
  '2_4_8_study_RZ',
  'Results_Aug2017',
  'krakenConcat.tsv')),
  stringsAsFactors = FALSE
  )

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
# Make this more functional

amrResultsTidy$Depth <- str_replace(amrResultsTidy$Sample, "\\d+_","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "\\d$","")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "full", "D1")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "half", "D0.5")
amrResultsTidy$Depth <- str_replace(amrResultsTidy$Depth, "quar*", "D0.25")

# Load filtered AMR data --------------------------------------------------

# Start here if AMR results have been filtered already

# Tidy the dataset

amrResultsTidy <- amrResultsFiltered %>% 
  gather(Category, CategoryName, c(1:4))

amrResultsTidy$Category <- factor(amrResultsTidy$Category, 
                                  levels = c('Class', 
                                             'Mechanism', 
                                             'Group', 
                                             'Gene'))

amrCategories <- levels(amrResultsTidy$Category)

krakenTaxa <- levels(krakenResultsFiltered$TaxRank)

# Split results by categories ---------------------------------------------

# Will adopt the split, apply, combine strategy
# Using the names of the AMR categories and taxonomic ranks as names of the 
# elements of the lists

amrResultsList <- amrResultsTidy %>% 
  split(.$Category) %>% #splitting dataframe using base R but with purrr's magrittr piping
  set_names(nm=amrCategories)


krakenResultsList <- krakenResultsFiltered %>% 
  split(.$TaxRank) %>% 
  set_names(nm=krakenTaxa)


# Summarize results(sum) --------------------------------------------------

amrResultsSummary <- lapply(amrResultsList, function(x){
  summarizeAMRbyCategory(x)
})

krakenResultsSummary <- lapply(krakenResultsList, function(x){
  summarizeKrakenbyTaxID(x)
})

# Convert results to wide format ------------------------------------------

# Potentially can use the same function and a list/vector of AMR and Kraken 
# results to avoid duplication

amrResultsWide <- lapply(amrResultsSummary, function(x){
  widenAMR(x)
})

krakenResultsWide <- lapply(krakenResultsSummary, function(x){
  widenKraken(x)
})


# Convert wide to matrix --------------------------------------------------

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


# CSS Normalization -------------------------------------------------------

# Transform dataframes into wide format
# Fill NAs with a value of zero

amrResultsAnalytical <- amrResultsFiltered %>%
  select(Gene, Hits, Sample) %>%
  spread(Sample, Hits, fill = 0, convert = TRUE)

krakenResultsAnalytical <- krakenResultsFiltered %>%
  select(TaxLineage, CladeReads, Sample) %>%
  spread(Sample, CladeReads, fill = 0, convert = TRUE)

amrAnalyticalMatrix <- matrixAMRanalytical(amrResultsAnalytical)

krakenAnalyticalMatrix <- matrixKraken(krakenResultsAnalytical)

amrExp <- newMRexperiment(amrAnalyticalMatrix[rowSums(amrAnalyticalMatrix) > 0, ])

krakenExp <- newMRexperiment(krakenAnalyticalMatrix[rowSums(krakenAnalyticalMatrix) > 0,])

cumNorm(amrExp)

cumNorm(krakenExp)

amrRaw <- data.frame(MRcounts(amrExp, norm=F))

amrNorm <- data.frame(MRcounts(amrExp, norm=T))

krakenRaw <- data.frame(MRcounts(krakenExp, norm=F))

krakenNorm <- data.frame(MRcounts(krakenExp, norm=T))

amrGenes <- row.names(amrAnalyticalMatrix) 

amrAnnotations <- read_tsv('amr_genes.tabular_parsed.tab')

amrNorm <- cbind(amrNorm, amrAnnotations)

krakenNorm$lineage <- row.names(krakenAnalyticalMatrix)

krakenNorm <- krakenNorm %>% separate(lineage, c('Domain',
                                   'Phylum',
                                    'Class',
                                    'Order',
                                    'Family',
                                    'Genus',
                                    'Species'),
                        sep = "|", fill = "right")




amr_class <- amr_norm[, lapply(.SD, sum), by='class', .SDcols=!c('header', 'mechanism', 'group')] 

amr_class_analytic <- newMRexperiment(counts=amr_class[, .SD, .SDcols=!'class']) rownames(amr_class_analytic) <- amr_class$class

# Construction of rarefaction curves --------------------------------------

# This step requires optimization (parallel mapping or compiling, maybe?)

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

# purrr version

#amrRarefy <- map(amrResultsMat, function(x){
#  raremax <- min(rowSums(x))
#  rarecurve_ROP(x, step=5, sample=raremax)
#})

# mclapply version (and the winner!)

amrRarCurve <- mclapply(amrResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=50, sample=raremax)
},mc.cores=10)

krakenRarCurve <- mclapply(krakenResultsMat, function(x){
  raremax <- min(rowSums(x))
  rarecurve(x, step=1000, sample=raremax)
},mc.cores=10)

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

amrRarCurveDF <- map(amrRarCurve, function(x) as_tibble(x, attr(x, "names")))

krakenRarefyDF <- map(krakenRarCurve, function(x) as_tibble(x, attr(x, "names")))

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

amrRarCurveDF$AMRLevel <- factor(amrRarCurveDF$AMRLevel, 
                                  levels = c('Class', 'Mechanism', 'Group', 'Gene'))

amrAlphaRarefactionDF$Level <- factor(amrAlphaRarefactionDF$Level, 
                                 levels = c('Class', 'Mechanism', 'Group', 'Gene'))


# Attempt to more functional programming
 
krakenRarefyDF <- krakenRarefyDF %>% filter(!krakenLevel %in% c("-", "D"))

krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "P", "Phyla")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "C", "Classes")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "O", "Orders")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "F", "Families")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "G", "Genera")
krakenRarefyDF$krakenLevel <- str_replace(krakenRarefyDF$krakenLevel, "S", "Species")

amrAllList <- amrRarCurveDF %>% 
  split(.$AMRLevel)

amrAllRarCurves <- amrAllList %>%
  map(function(x){
    amrRarefactionCurve(x)
  })

krakenAllRarCurveList <- krakenRarefyDF %>% 
  split(.$krakenLevel)

krakenAllRarCurvesLarger <- krakenAllRarCurveList %>%
  map(function(x){
    krakenRarefactionCurve(x)
  })

foreach(i=krakenAllRarCurvesLarger) %do%
  ggsave(filename=paste('rarefaction',unique(i$data$krakenLevel),'CB','noFacet','.png', sep='', collapse=''),
         path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/rarefaction',
         plot = i,
         height=8.50,
         width=10.50,
         units="in",
         device="png")

amrAllAlphaBoxPlots <- amrAlphaRarefactionDF %>%
  amrAlphaDiv()

ggsave('~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity/amrAlphaDiversityCB.png',
       plot = amrAllAlphaBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")


# Alpha Diversity (rarefied) and Species Richness calculations -------------

krakenAlphaRarefaction <- mclapply(krakenResultsMat, function(x){
  alpha_rarefaction(x, minlevel=0)
},mc.cores = 10)

amrAlphaRarefaction <- mclapply(amrResultsMat, function(x){
  alpha_rarefaction(x, minlevel=0)
}, mc.cores=10)

# Consider change to data frame/tibble 

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

foreach(i=amrAlphaRarefaction, j=amrCategories) %do%
  rep(j,length(i$ID)) -> i$Level

#amrAlphaRarefaction[[1]]$Level <- rep(amrCategories[[1]], length(amrAlphaRarefaction[[1]]$ID))
#amrAlphaRarefaction[[2]]$Level <- rep(amrCategories[[2]], length(amrAlphaRarefaction[[2]]$ID))
#amrAlphaRarefaction[[3]]$Level <- rep(amrCategories[[3]], length(amrAlphaRarefaction[[3]]$ID))
#amrAlphaRarefaction[[4]]$Level <- rep(amrCategories[[4]], length(amrAlphaRarefaction[[4]]$ID))

foreach(i=krakenAlphaRarefaction, j=krakenTaxa) %do%
  rep(j,length(i$ID)) -> i$Level


foreach(i=krakenAlphaRarefaction2, j=krakenTaxa) %do%
  rep(j,length(i$ID)) -> i$Level

# Seems like unnecessary repetition to convert these to data frames when they
# could have been created as dataframes

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

amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$ID, "_.*", "")

krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$ID, "_.*", "")


krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$ID, "_.*", "")
# Generate one single dataframe and create column for AMR Level

krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "F", "D1")
krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "H", "D0.5")
krakenAlphaRarefactionDF$Depth <- str_replace(krakenAlphaRarefactionDF$Depth, "QD", "D0.25")

amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "F", "D1")
amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "H", "D0.5")
amrAlphaRarefactionDF$Depth <- str_replace(amrAlphaRarefactionDF$Depth, "QD", "D0.25")

krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "F", "D1")
krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "H", "D0.5")
krakenAlphaRarefaction2DF$Depth <- str_replace(krakenAlphaRarefaction2DF$Depth, "QD", "D0.25")


# Alpha Diversity (rarefied) and Species Richness boxplots -------------------

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
                                  levels = c('Phyla', 
                                             'Classes', 
                                             'Orders', 
                                             'Families', 
                                             'Genera', 
                                             'Species'))

amrAllSpRawBoxPlots <- amrAlphaRarefactionDF %>%
    amrRawSpeciesRich()

ggsave(filename = 'amrSpeciesRichnessCB.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/amrResults_Aug2017_75_gene_frac/alphaDiversity',
       plot = amrAllSpRawBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenAllAlphaBoxPlots <- krakenAlphaRarefaction2DF %>%
    krakenAlphaDiv()

ggsave(filename = 'krakenAlphaDivCB.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenAllAlphaBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenAllShannonBoxPlots <- krakenAlphaRarefaction2DF %>%
  krakenShannon()

ggsave(filename = 'krakenAlphaDivCB.png',
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenAllAlphaBoxPlots,
       width = 10.50,
       height = 8.50,
       units = "in")

krakenAllSpRawBoxPlots <- krakenAlphaRarefaction2DF %>%
  krakenRawSpeciesRich()

ggsave(filename = 'krakenSpRichnessCB.png' ,
       path = '~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017/alphaDiversity',
       plot = krakenAllSpRawBoxPlots,
       width = 10.50,
       height = 8.50,
       units="in")

# Correlation plots -------------------------------------------------------

# Should make into a function for the package

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
             color="grey") 
  # geom_smooth(aes(group=1, weight=0.2), 
  #             method="lm", 
  #             se=FALSE, 
  #             colour="grey", 
  #             alpha=0.5) 

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

# Generating kraken correlation plot

krakenReadsvsHitsCorTest <- cor.test(x = krakenReadstoHitRatio$Reads, krakenReadstoHitRatio$KrakenHits, method = "spearman")

krakenReadsvsHitsCor <- ggplot(krakenReadstoHitRatio, aes(Reads, KrakenHits)) + 
  geom_point(aes(fill=Sample_type), alpha=0.6, size=10, pch=21, color="grey") 
  # geom_smooth(aes(group=1, weight=0.2), method="lm", se=FALSE, colour="grey", alpha=0.5) 

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

ggsave(filename='krakenReadsvsHitsCorCBScheme.png', 
       path='~/amr/2-4-8_results/2_4_8_study_RZ/krakenResults_Aug2017',
       width=14, 
       height=8.50,
       units="in")

# Kruskal-Wallis tests ----------------------------------------------------

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

# Computing aggregated sums of hits ---------------------------------------

krakenPhylumResults <- krakenResultsList[["P"]]

krakenPhylumResults$SampleID <- str_replace(krakenPhylumResults$Sample, "_00[6-8]")

krakenPhylumSums <- krakenPhylumResults %>%
  group_by(TaxID, SampleID, Sample_Type) %>%
  summarise(MeanHits = mean(CladeReads)) %>%
  group_by(SampleID) %>%
  summarise(SumHits = sum(MeanHits))
