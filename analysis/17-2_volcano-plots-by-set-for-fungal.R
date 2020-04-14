library(dplyr)
library(stringr)
library(tidyr)

source("analysis/99_make-volcano-plot.R")

# abs(logFC) log fold change cut-off.  Anything greater 
# than (-1 x lfc) and less than lfc will be deemed 
# biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed 
# statistically insignificant.
pCutoff <- 0.05


# Read top BLAST hit and taxonomy lineage
topBlastHit <- read.delim("results/all-DE.blastpUniProt.topHit.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)
fungalCtgIds <- topBlastHit %>% filter(taxonomy == "Fungi") %>%
  select(qseqid)
fungalCtgIds <- as.character(fungalCtgIds$qseqid)
length(fungalCtgIds)
# [1] 1363


# Get all the statistic results from DE analysis
allStats <- read.delim("results/stats-results.long.22Feb.tsv", header = TRUE, 
                       row.names = NULL, stringsAsFactors = FALSE)
allStats <- allStats %>% select(-contig, -AveExpr, -t, -P.Value, -B) %>%
  filter(focus_term == "constDiff" | focus_term == "weevilInd_Q903") %>%
  filter(cds %in% fungalCtgIds)



# Focus will be on cd_vs_wiq
cd_vs_wiq <- read.delim("results/cd_vs_wiq.txt", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE) %>% select(cds, set)
str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  2 variables:


set12 <- cd_vs_wiq %>% filter(set == 1 | set == 2)
set12 <- subset(allStats, allStats$cds %in% set12$cds)

set12_vp <- makeVolcanoPlot(set12, pCutoff, lfcCutoff, "Set 1/2 Fungal Draft")
ggsave("results/figures/set12_volcano_plot_fungal.draft.svg", 
       plot = set12_vp, height = 4, width = 6)

# set 3/4
set34 <- cd_vs_wiq %>% filter(set == 3 | set == 4)
set34 <- subset(allStats, allStats$cds %in% set34$cds)

set34_vp <- makeVolcanoPlot(set34, pCutoff, lfcCutoff, "Set 3/4 Fungal Draft")
ggsave("results/figures/set34_volcano_plot_fungal.draft.svg", 
       plot = set34_vp, height = 4, width = 6)

#set5
set5 <- cd_vs_wiq %>% filter(set == 5)
set5 <- subset(allStats, allStats$cds %in% set5$cds)

set5_vp <- makeVolcanoPlot(set5, pCutoff, lfcCutoff, "Set 5 Fungal Draft")
ggsave("results/figures/set5_volcano_plot_fungal.draft.svg", 
       plot = set5_vp, height = 4, width = 6)
