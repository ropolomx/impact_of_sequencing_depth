library(dplyr)
library(tidyr)
library(metagenomeSeq)
library(vegan)
library(stringr)

source('rarefaction_functions.R')

amrResults <- read.csv('AMR/amr_new_dataframe_ROP.csv')

amrLevels <- c("Class", "Mechanism", "Group")

  
#amr_results_class <- amr_results %>% group_by(Sample, Class) %>% summarise(Hits=sum(Hits_Seen))

#amr_results_class_wide <- amr_results_class %>% spread(Class, Hits, fill = 0)
#amr_results_class_wide <- amr_results_class %>% spread(Sample, Hits, fill = 0)
#table(annotations$class)


#amr_results_class_wide <- amr_results_class %>% spread(Class, Hits, fill = 0)
#
#amr_results_class_mat <- amr_results_class_wide[,2:ncol(amr_results_class_wide)]

amr_class_norm <- data.frame(MRcounts(amr_class_experiment, norm=T))
amr_class_norm <- t(amr_class_norm)
row.names(amr_class_norm) <- str_replace(row.names(amr_class_norm), "X", "")

alpha_rarefaction <- function(X, step=0.05, method='invsimpson') {
    S <- specnumber(X, MARGIN=2)
    raremax <- min(colSums(X))
    tot <- colSums(X)
    nr <- ncol(X)
    out <- lapply(seq_len(nr), function(i){
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
                out=out))
}

class_samples <- attr(ggplot_rarefy$raw_species_abundance, "names")
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
