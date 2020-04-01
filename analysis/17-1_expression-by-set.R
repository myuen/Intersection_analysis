library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# abs(logFC) log fold change cut-off.  Anything greater 
# than (-1 x lfc) and less than lfc will be deemed 
# biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed 
# statistically insignificant.
pCutoff <- 0.05


allStats <- read.delim("results/stats-results.long.22Feb.tsv", header = TRUE, 
                       row.names = NULL, stringsAsFactors = FALSE)
allStats <- allStats %>% select(-contig, -AveExpr, -t, -P.Value, -B) %>%
  filter(focus_term == "constDiff" | focus_term == "weevilInd_Q903")


# Focus will be on cd_vs_wiq
cd_vs_wiq <- read.delim("results/cd_vs_wiq.txt", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE) %>% select(cds, set)
str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  2 variables:

table(cd_vs_wiq$set)
#    1    2    3    4    5    6    7    8 
# 2065 1514 1855 2975   85   69  104  631 



# Function to make Volcano plot
volcanoPlot <- function(df, pTitle) {
  
  minX <- round(min(df$logFC))
  maxX <- round(max(df$logFC))
  
  g <- ggplot(
    df, 
    aes(x = logFC, 
        y = -log10(adj.P.Val),
        color = focus_term,
        shape = focus_term)) + 
    geom_point(size = 0.9) + 
    
    scale_colour_manual(values = c("black", "grey")) +
    
    scale_x_continuous(breaks = c(minX, -2, 0, 2, maxX)) +
    
    scale_y_continuous(breaks = c(0, 1.3), labels = c(0, 0.05)) +
    
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    
    geom_abline(aes(intercept = -log10(pCutoff), slope = 0), colour = "blue") + 
    geom_vline(xintercept = lfcCutoff, colour = "blue") + 
    geom_vline(xintercept = -(lfcCutoff), colour = "blue") +
    
    labs(title = pTitle,
         x = expression(paste('Log'[2], " Fold Change")), 
         y = expression(paste('-Log'[10], " Adj. p-Value"))) +
    
    guides(colour = guide_legend("", keywidth = .5),
           shape = guide_legend("", keywidth = .5))
  
  return (g)
}

set12 <- cd_vs_wiq %>% filter(set == 1 | set == 2)
set12 <- subset(allStats, allStats$cds %in% set12$cds)

set12_vp <- volcanoPlot(set12, "Set 1/2 Draft")
ggsave("results/figures/set12_volcano_plot.draft.svg", 
       plot = set12_vp, height = 4, width = 6)


# set 3/4
set34 <- cd_vs_wiq %>% filter(set == 3 | set == 4)
set34 <- subset(allStats, allStats$cds %in% set34$cds)

set34_vp <- volcanoPlot(set34, "Set 3/4 Draft")
ggsave("results/figures/set34_volcano_plot.draft.svg", 
       plot = set34_vp, height = 4, width = 6)


#set5
set5 <- cd_vs_wiq %>% filter(set == 5)
set5 <- subset(allStats, allStats$cds %in% set5$cds)

set5_vp <- volcanoPlot(set5, "Set 5 Draft")
ggsave("results/figures/set5_volcano_plot.draft.svg", 
       plot = set5_vp, height = 4, width = 6)

#set6
set6 <- cd_vs_wiq %>% filter(set == 6)
set6 <- subset(allStats, allStats$cds %in% set6$cds)

set6_vp <- volcanoPlot(set6, "Set 6 Draft")
ggsave("results/figures/set6_volcano_plot.draft.svg", 
       plot = set6_vp, height = 4, width = 6)

#set7
set7 <- cd_vs_wiq %>% filter(set == 7)
set7 <- subset(allStats, allStats$cds %in% set7$cds)

set7_vp <- volcanoPlot(set7, "Set 7 Draft")
ggsave("results/figures/set7_volcano_plot.draft.svg", 
       plot = set7_vp, height = 4, width = 6)

# set8
set8 <- cd_vs_wiq %>% filter(set == 8)
set8 <- subset(allStats, allStats$cds %in% set8$cds)

set8_vp <- volcanoPlot(set8, "Set 8 Draft")
ggsave("results/figures/set8_volcano_plot.draft.svg", 
       plot = set8_vp, height = 4, width = 6)
