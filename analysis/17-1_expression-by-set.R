library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# Read normalized CPM
cpm <- read.table("results/normalized_cpm.22Feb.txt")
str(cpm)
# 'data.frame':	38711 obs. of  24 variables:


# Focus will be on cd_vs_wiq
cd_vs_wiq <- read.delim("results/cd_vs_wiq.txt", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE) %>% select(cds, set)
str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  2 variables:

table(cd_vs_wiq$set)
#    1    2    3    4    5    6    7    8 
# 2065 1514 1855 2975   85   69  104  631 

cd_vs_wiq.cpm <- merge(cd_vs_wiq, cpm, by.x = "cds", by.y = 0)


# Need to reshape data into a format that can be read by ggplot.
# One variable per column, One reading per row
by_set <- cd_vs_wiq.cpm %>% gather("library", "cpm", -set, -"cds") %>% 
  group_by(set, library)

# Create genoptyep column by reading library name
# H898 = R; Q903 = S
by_set <- by_set %>% mutate(
  genotype = case_when(
    str_detect(library, "H898") ~ "R", 
    str_detect(library, "Q903") ~ "S"
  )
)
by_set$genotype <- factor(by_set$genotype)
by_set$genotype <- relevel(by_set$genotype, ref = "S")


# Create treatment column by reading library name
by_set$treatment <- 
  factor(str_sub(by_set$library, 5, 5), levels = c("C", "W", "G"))
by_set$treatment <- relevel(by_set$treatment, ref = "C")


makeBoxplot <- function(data, filename) {
  ggplot(data, aes(x = treatment, y = cpm, color = genotype)) + 
    geom_boxplot(outlier.shape = 1) + theme_classic() + 
    xlab("Treatment") + ylab("CPM") + labs(color='Genoytpe') + 
    scale_colour_manual(values = c("grey", "black")) + 
    ggsave(paste0("results/figures/", filename, ".svg"),
           width = 20, height = 20, units = c("cm"), dpi = 300)
}

makeBoxplot(by_set %>% filter(set == 1), "set1_exp_boxplot.24Mar")

makeBoxplot(by_set %>% filter(set == 2), "set2_exp_boxplot.24Mar")

makeBoxplot(by_set %>% filter(set == 3), "set3_exp_boxplot.24Mar")

makeBoxplot(by_set %>% filter(set == 4), "set4_exp_boxplot.24Mar")

makeBoxplot(by_set %>% filter(set == 5), "set5_exp_boxplot.24Mar")

makeBoxplot(by_set %>% filter(set == 6), "set6_exp_boxplot.24Mar")

makeBoxplot(by_set %>% filter(set == 7), "set7_exp_boxplot.24Mar")

makeBoxplot(by_set %>% filter(set == 8), "set8_exp_boxplot.24Mar")
