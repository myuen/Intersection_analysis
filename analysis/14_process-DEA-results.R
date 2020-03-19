sidf <- read.delim("results/stats-results.long.22Feb.tsv", 
                   colClasses = c(rep("character", 3), rep("numeric", 6)))
str(sidf)
# 'data.frame':	270977 obs. of  9 variables:


# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfc) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05


# Get all DE IDs
de <- subset(sidf, abs(sidf$logFC) >= lfcCutoff & sidf$adj.P.Val <= pCutoff)
dim(de)
# [1] 17507     9

# write out DE results
write.table(de, "results/all-DE.stats.txt", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# write out id of all differentially expressed contigs
write(x = unique(sort(de$cds)), "results/all-DE.id.txt")
