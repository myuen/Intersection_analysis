library(dplyr)
library(gplots)
library(tibble)
library(tidyr)
library(RColorBrewer)


# Read set analysis
cd_vs_wiq <- read.delim("results/cd_vs_wiq.txt", sep = "\t", 
                        header = TRUE, stringsAsFactors = FALSE)
cd_vs_wiq <- cd_vs_wiq %>% select(cds, set)
str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  2 variables:


# Read GO Terms
goTerms <- read.delim("results/blastHitGoTerms.txt", sep = "\t",
                      header = TRUE, stringsAsFactors = FALSE)
str(goTerms)
# 'data.frame':	3662 obs. of  5 variables:


cd_vs_wiq.tax.go <- inner_join(cd_vs_wiq, goTerms)
str(cd_vs_wiq.tax.go)
# 'data.frame':	2824 obs. of  6 variables:

write.table(cd_vs_wiq.tax.go, "results/cd_vs_wiq.go.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)


# Summarize GO term count
goCount <- cd_vs_wiq.tax.go  %>% group_by(set, TERM, Taxonomy) %>% 
  summarise(count = n()) %>% 
  filter(Taxonomy == "Viridiplantae" | Taxonomy == "Fungi")
str(goCount)
# Classes ‘grouped_df’, ‘tbl_df’, ‘tbl’ and 'data.frame':	935 obs. of  4 variables:


goCountSum <- goCount %>% group_by(set, Taxonomy) %>% 
  summarise(sum = sum(count))


goCount <- inner_join(goCount, goCountSum)
# Joining, by = c("set", "Taxonomy")
str(goCount)
# Classes ‘grouped_df’, ‘tbl_df’, ‘tbl’ and 'data.frame':	935 obs. of  5 variables:


# Normalize go term counts
goCount <- goCount %>% mutate(perc = count / sum * 100) %>%
  select(-count, -sum)


# Set heatmap color
coul <- colorRampPalette(brewer.pal(8, "RdYlBu"))(25)

makeHM <- function(data, filename) {
  svg(file = paste0("results/figures/", filename, ".svg"),
      width = 12, height = 12,
      pointsize = 16, family = "helvetica")
  
  heatmap.2(
    as.matrix(data), 
    Colv = NA, 
    scale = "row", 
    dendrogram = c("row"), 
    trace = c("none"),
    col = coul,
    # cexRow = 0.75, cexCol = 1, # Row and Column label size
    cexRow = 0.1, cexCol = 1, # Row and Column label size
    srtCol = 45, # Col label angle
    margins = c(3, 20),
    keysize = 1, # Legend size
    density.info = c("none"),
    lhei = c(0.5, 3),
    lwid = c(1, 4)
  )
  
  dev.off()
}


# Select GO terms from set 1 and 2 from Viridiplantae
# Reformat data to be feed into heatmap.2
# One term per row, two columns from normalize count from two sets
goTerms_set_1_2_plants <- goCount %>%
  filter(Taxonomy == "Viridiplantae" & (set == 1 | set == 2)) %>%
  select(-Taxonomy) %>%
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>% 
  column_to_rownames("TERM")

# write.table(goTerms_set_1_2_plants, "results/goTerms_set_1_2_plants.txt", 
#             quote = FALSE, sep = "\t")

makeHM(goTerms_set_1_2_plants, "go_set12_plant_heatmap.30Mar")


# Repeat for GO terms from set 3/4
goTerms_set_3_4_plants <- goCount %>%
  filter(Taxonomy == "Viridiplantae" & (set == 3 | set == 4)) %>%
  select(-Taxonomy) %>%
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>%
  column_to_rownames("TERM")

# write.table(goTerms_set_3_4_plants, "results/goTerms_set_3_4_plants.txt", 
#             quote = FALSE, sep = "\t")

makeHM(goTerms_set_3_4_plants, "go_set34_plant_heatmap.30Mar")


# Repeat for GO terms from set 5-8
goTerms_set_5_8_plants <- goCount %>%
  filter(Taxonomy == "Viridiplantae" & (set >= 5 & set <= 8)) %>%
  select(-Taxonomy) %>% 
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>%
  column_to_rownames("TERM")

# write.table(goTerms_set_5_8_plants, "results/goTerms_set_5_8_plants.txt", 
#             quote = FALSE, sep = "\t")

makeHM(goTerms_set_5_8_plants, "go_set58_plant_heatmap.30Mar")


# Make a heatmap for GOslim terms for all sets for fungi
goTerms_all_sets_fungi <- goCount %>%
  filter(Taxonomy == "Fungi") %>%
  select(-Taxonomy) %>%
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>%
  column_to_rownames("TERM")

# write.table(goTerms_all_sets_fungi, "results/goTerms_all_sets_fungi.txt", 
#             quote = FALSE, sep = "\t")

makeHM(goTerms_all_sets_fungi, "go_all_sets_fungi_heatmap.30Mar")
