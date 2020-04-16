library(dplyr)

# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfc) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05


sidf <- read.delim("results/stats-results.long.22Feb.tsv", 
                   colClasses = c(rep("character", 3), rep("numeric", 6)))
str(sidf)
# 'data.frame':	270977 obs. of  9 variables:


# Get all DE IDs
de <- sidf %>% 
  filter(abs(logFC) >= lfcCutoff & adj.P.Val <= pCutoff) %>%
  select(-contig, -AveExpr, -t, -B)

# We are focusing on only constDiff and weevilInd_Q903 
# for this study
de <- de %>% 
  filter(focus_term == "constDiff" | focus_term == "weevilInd_Q903")

str(de)
# 'data.frame':	10187 obs. of  5 variables:

# write out DE results
write.table(de, "results/sigDE.stats.txt", quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# write out id of all differentially expressed contigs
# write(x = unique(sort(de$cds)), "results/sigDE.ctgId.txt")
