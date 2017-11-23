# contingency plan


amrSummariesMechanism <- summarizeAMRlevels(amrResults, amrLevel = "Mechanism")

amrSummariesGroup <- summarizeAMRlevels(amrResults, amrLevel = "Group")

amrMatGroup <- matrixAMR(amrWideGroup)


# Rarefaction strategy ----------------------------------------------------

rarefyGroupDF <- lapply(rarefyGroup$rarefy_out, unlist)
rarefyGroupDF <- lapply(rarefyGroupDF, data.frame)
rarefyGroupDF <- do.call("rbind", rarefyGroupDF[1:32])


names(rarefyGroupDF) <- "Groups"
rarefyGroupDF$Samples <- rep(sampleNames, each=21)
rarefyGroupDF$Sampling_size <- str_replace(row.names(rarefyGroupDF), "N","")
rarefyGroupDF$Sampling_size <- as.numeric(as.character
                                          (rarefyGroupDF$Sampling_size))

rarefyGroupDF$Sample_type <- str_replace(rarefyGroupDF$Samples, "\\d+_", "")

rarefyGroupDF$Sample_type <- str_replace(rarefyGroupDF$Sample_type, "\\d$","")

rarefyGroupDF$Sample_number <- str_replace(rarefyGroupDF$Samples, "_.*", "")


ggplot(rarefyGroupDF, aes(Sampling_size,Groups)) + geom_line(aes(color=Samples)) + facet_grid(Sample_type~ Sample_number)



meg_bar <- ggplot(bar_subset, aes_string(x=group_var, y='Normalized_Count', fill='Name')) +
        geom_bar(stat='identity') + 
        scale_fill_brewer(palette="Spectral") +
        theme(strip.text.x=element_text(size=18),
              axis.text.y=element_text(size=20),
              axis.text.x=element_text(size=22, vjust=1, hjust=1, angle=33),
              axis.title.x=element_text(size=24),
              axis.title.y=element_text(size=24),
              legend.title=element_text(size=20),
              legend.text=element_text(size=18),
              plot.title=element_text(size=26, hjust=0.25)) +
        xlab(group_var) +
        ylab('Mean of Normalized Count\n') +
        ggtitle(paste('Mean ', data_type, ' ', level_var, ' Normalized Count by ', group_var, '\n',
                      sep='', collapse=''))
    png(filename=paste(outdir, '/', data_type, '_', level_var, '_BarPlot_by_', group_var, '.png',
                       sep='', collapse=''), width=1024, height=768)
    print(meg_bar)
    dev.off()
}