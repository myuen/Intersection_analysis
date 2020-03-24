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

table(cd_vs_wiq$set)
#    1    2    3    4    5    6    7    8 
# 2065 1514 1855 2975   85   69  104  631 


# Read top BLAST hit and taxonomy lineage
topBlastHit <- read.delim("results/all-DE.blastpUniProt.topHit.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)
topBlastHit <- topBlastHit %>% select(qseqid, taxonomy)
str(topBlastHit)
# 'data.frame':	9714 obs. of  2 variables:

cd_vs_wiq.tax <- merge(cd_vs_wiq, topBlastHit, by.x = "cds", by.y = "qseqid")
str(cd_vs_wiq.tax)
# 'data.frame':	7928 obs. of  3 variables:
# Some sequence don't have annotation, hence size not the same as cd_vs_wiq


# Read top GOslim annotations
topTairGo <- read.delim("results/all-DE.blastpTAIR.topHit.txt", sep = "\t", 
                        header = TRUE, stringsAsFactors = FALSE)
topTairGo <- topTairGo %>% select(qseqid, GOslim_term)
str(topTairGo)
# 'data.frame':	18655 obs. of  2 variables:

cd_vs_wiq.go <- merge(cd_vs_wiq.tax, topTairGo, by.x = "cds", by.y = "qseqid")
str(cd_vs_wiq.go)
# 'data.frame':	15309 obs. of  4 variables:
# Some sequence have more than 1 go term, hence size not the 
# same as cd_vs_wiq.tax


# Summarize taxonomy count
goCount <- cd_vs_wiq.go  %>% group_by(set, GOslim_term, taxonomy) %>% 
  summarise(count = n())

# Calculate the sum of go term counts by set
goCountSum <- goCount %>% group_by(set) %>% summarise(sum = sum(count))

goCount <- merge(goCount, goCountSum, by.x = "set", by.y = "set")

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
    as.matrix(data), Colv = NA, 
    scale = "row", dendrogram = c("row"), trace = c("none"),
    col = coul,
    # density.info = c("none"),
    keysize = 1, # Legend size
    cexRow = 0.75, cexCol = 0.85, # Row and Column label size
    srtCol = 45, # Col label angle
    margins = c(5, 20)
  )
  
  dev.off()
}


# Make heatmap for GOslim terms in set 1/2 for plant and fungi separately

# Select GOslim terms from set 1 and 2 from Viridiplantae
goTerms_set_1_2_plants <- goCount %>% filter(set == 1 | set == 2) %>%
  filter(taxonomy == "Viridiplantae") %>% select(-taxonomy)

# Reformat data to be feed into heatmap.2
# One term per row, two columns from normalize count from two sets
goTerms_set_1_2_plants <- goTerms_set_1_2_plants %>%
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>%
  column_to_rownames("GOslim_term")

makeHM(goTerms_set_1_2_plants, "go_set12_plant_heatmap.24Mar")

# Set 1/2 Fungi
goTerms_set_1_2_fungi <- goCount %>% filter(set == 1 | set == 2) %>% 
  filter(taxonomy == "Fungi") %>% select(-taxonomy) 

goTerms_set_1_2_fungi <- goTerms_set_1_2_fungi %>% 
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>% 
  column_to_rownames("GOslim_term")

makeHM(goTerms_set_1_2_fungi, "go_set12_fungi_heatmap.24Mar")


# Make heatmap for GOslim terms in set 3/4 for plant and fungi separately
goTerms_set_3_4_plants <- goCount %>% filter(set == 3 | set == 4) %>% 
  filter(taxonomy == "Viridiplantae") %>% select(-taxonomy) 

goTerms_set_3_4_plants <- goTerms_set_3_4_plants %>% 
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>% 
  column_to_rownames("GOslim_term")

makeHM(goTerms_set_3_4_plants, "go_set34_plant_heatmap.24Mar")

# Set 3/4 Fungi
goTerms_set_3_4_fungi <- goCount %>% filter(set == 3 | set == 4) %>% 
  filter(taxonomy == "Fungi") %>% select(-taxonomy) 

goTerms_set_3_4_fungi <- goTerms_set_3_4_fungi %>% 
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>% 
  column_to_rownames("GOslim_term")

makeHM(goTerms_set_3_4_fungi, "go_set34_fungi_heatmap.24Mar")


# Make heatmap for GOslim terms in set 5 to 8 for plant 
# There are no Fungi gene found in set 6-8, so only 1 heatmap for plant

goTerms_set_5_8_plants <- goCount %>% filter(set >= 5 & set <= 8) %>%
  filter(taxonomy == "Viridiplantae") %>% select(-taxonomy) 

goTerms_set_5_8_plants <- goTerms_set_5_8_plants %>%
  pivot_wider(names_from = set, values_from = perc, names_prefix = "set",
              values_fill = list(perc = 0)) %>% 
  column_to_rownames("GOslim_term")

makeHM(goTerms_set_5_8_plants, "go_set58_plant_heatmap.24Mar")
