library(dplyr)
library(tibble)

# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfc) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05


# Read differential expression statistics
de <- read.table("results/stats-results.wide.22Feb.tsv",
                 header = TRUE, stringsAsFactors = FALSE ,
                 colClasses = c(rep("character", 2), rep("numeric", 14)))

de <- de %>% rownames_to_column("cds")


# Focus will be on comparing constDiff and weevilInd_Q903
de <- de %>% select(cds, logFC.constDiff, adj.P.Val.constDiff,
                    logFC.weevilInd_Q903, adj.P.Val.weevilInd_Q903)

str(de)
# 'data.frame':	38711 obs. of  5 variables:

cd_vs_wiq <- de %>% 
  filter((abs(logFC.constDiff) >= lfcCutoff & adj.P.Val.constDiff <= pCutoff) | 
           (abs(logFC.weevilInd_Q903) >= lfcCutoff & adj.P.Val.weevilInd_Q903 <= pCutoff))

str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  5 variables:


# Select contigs with higher expression in resistance genotype (H898)
cdUp <- cd_vs_wiq %>% 
  filter(logFC.constDiff >= lfcCutoff & adj.P.Val.constDiff <= pCutoff) %>% 
  select(cds)

cdUp <- as.character(cdUp$cds)

length(cdUp)
# [1] 2254


# Select contigs with higher expression in susceptible genotype (Q903)
cdDown <- cd_vs_wiq %>% 
  filter(logFC.constDiff <= (-1 *lfcCutoff) & adj.P.Val.constDiff <= pCutoff) %>% 
  select(cds)

cdDown <- as.character(cdDown$cds)

length(cdDown)
# [1] 2214


# Select contigs with up-regulated expression in susceptible genotype (Q903) 
# during insect feeding
wiqUp <- cd_vs_wiq %>% 
  filter(logFC.weevilInd_Q903 >= lfcCutoff & adj.P.Val.weevilInd_Q903 <= pCutoff) %>%
  select(cds)

wiqUp <- as.character(wiqUp$cds)

length(wiqUp)
# [1] 2009


# Select contigs with down-regulated expression in susceptible genotype (Q903) 
# during insect feeding
wiqDown <- cd_vs_wiq %>% 
  filter(logFC.weevilInd_Q903 <= (-1 * lfcCutoff) & adj.P.Val.weevilInd_Q903 <= pCutoff) %>%
  select(cds)

wiqDown <- as.character(wiqDown$cds)

length(wiqDown)
# [1] 3710


# The size of the total set
length(unique(c(cdUp, cdDown, wiqUp, wiqDown)))
# [1] 9298


# Split the data into different sets

# Set 1 = Higher expression in resistance genotype (H898)
# but not DE in insect feeding
set1 <- setdiff(cdUp, union(wiqUp, wiqDown))
cd_vs_wiq[cd_vs_wiq$cds %in% set1, "set"] <- 1

# Set 2 = Higher expression in susceptible genotype (Q903)
# but not DE in insect feeding
set2 <- setdiff(cdDown, union(wiqUp, wiqDown))
cd_vs_wiq[cd_vs_wiq$cds %in% set2, "set"] <- 2

# Set 3 = Up-regulated in susceptible genotype (Q903)
# during insect feeding but not DE constitutively
set3 <- setdiff(wiqUp, union(cdUp, cdDown))
cd_vs_wiq[cd_vs_wiq$cds %in% set3, "set"] <- 3

# Set 4 = Down-regulated in susceptible genotype (Q903)
# during insect feeding but not DE constitutively
set4 <- setdiff(wiqDown, union(cdUp, cdDown))
cd_vs_wiq[cd_vs_wiq$cds %in% set4, "set"] <- 4

# Set 5 = Up-regulated in susceptible genotype (Q903) during 
# insect feeding AND higher expression constitutively in 
# resistance genotype (H898)
set5 <- intersect(cdUp, wiqUp)
cd_vs_wiq[cd_vs_wiq$cds %in% set5, "set"] <- 5

# Set 6 = Up-regulated in susceptible genotype (Q903) during 
# insect feeding AND higher expression constitutively in 
# susceptible genotype (Q903)
set6 <- intersect(cdDown, wiqUp)
cd_vs_wiq[cd_vs_wiq$cds %in% set6, "set"] <- 6

# Set 7 = Down-regulated in susceptible genotype (Q903) during 
# insect feeding AND higher expression constitutively in 
# resistance genotype (H898)
set7 <- intersect(cdUp, wiqDown)
cd_vs_wiq[cd_vs_wiq$cds %in% set7, "set"] <- 7

# Set 8 = Down-regulated in susceptible genotype (Q903) during 
# insect feeding AND higher expression constitutively in
# susceptible genotype (Q903)
set8 <- intersect(cdDown, wiqDown)
cd_vs_wiq[cd_vs_wiq$cds %in% set8, "set"] <- 8


table(cd_vs_wiq$set)
#    1    2    3    4    5    6    7    8 
# 2065 1514 1855 2975   85   69  104  631 


write.table(cd_vs_wiq, "results/cd_vs_wiq.txt", quote = FALSE, 
            sep = "\t", col.names = TRUE, row.names = FALSE)
