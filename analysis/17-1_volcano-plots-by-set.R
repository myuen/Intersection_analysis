library(dplyr)
library(stringr)
library(tidyr)

source("analysis/helper04_make-volcano-plot.R")


# abs(logFC) log fold change cut-off.  Anything greater 
# than (-1 x lfc) and less than lfc will be deemed 
# biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed 
# statistically insignificant.
pCutoff <- 0.05


# Get all the statistic results from DE analysis
allStats <- read.delim("results/stats-results.long.22Feb.tsv", header = TRUE, 
                       row.names = NULL, stringsAsFactors = FALSE)

allStats <- allStats %>% select(-contig, -AveExpr, -t, -P.Value, -B) %>%
  filter(focus_term == "constDiff" | focus_term == "weevilInd_Q903")


# Focus will be on cd_vs_wiq
cd_vs_wiq <- read.delim("results/cd_vs_wiq.txt", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE) %>% select(cds, set)
str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  2 variables:


map(c(1:8), function(x){

  set <- cd_vs_wiq %>% filter(set == x)
  set <- subset(allStats, allStats$cds %in% set$cds)
  
  set_vp <- makeVolcanoPlot(set, pCutoff, lfcCutoff, 
                            paste("Set", x, "Draft, 1 Oct"))
  
  ggsave(paste0("results/figures/set", x, "_volcano_plot.1Oct.svg"), 
         plot = set_vp, height = 12, width = 12)
})
