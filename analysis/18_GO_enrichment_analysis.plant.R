library(dplyr)
library(purrr)
library(stringr)
library(topGO)

# Gene Ontology (GO) term enrichment analysis.

options(digits = 4)

# Read top BLAST hit with GOIDs
cd_vs_wiq.annot <-read.delim("results/cd_vs_wiq.annot.txt", 
                             sep = "\t", header = TRUE, 
                             stringsAsFactors = FALSE)

# Keep only Plant sequences
cd_vs_wiq.annot.plant <-
  cd_vs_wiq.annot %>%
  filter(Taxonomy == "Viridiplantae")

str(cd_vs_wiq.annot.plant)
# 'data.frame':	7425 obs. of  12 variables:


# Reformat into named character vector
cds_2_goid.ncv <-
  map(cd_vs_wiq.annot.plant$cds, function(c){
    cd_vs_wiq.annot.plant %>%
      filter(cds == c) %>%
      dplyr::select(GO) %>%
      str_split(";") %>%
      unlist() %>%
      str_trim()
  })

names(cds_2_goid.ncv) <- cd_vs_wiq.annot.plant$cds

# Remove CDS without GO annotations
cds_2_goid.ncv <- cds_2_goid.ncv[!is.na(cds_2_goid.ncv)]

length(cds_2_goid.ncv)
# [1] 5533

# Define selection function
upRegGenes <- function(x) {
  return(x > 0)
}

downRegGenes <- function(x) {
  return(x < 0)
}


# We are interest in what terms are enriched in Set 1 against Set 2
# Subset data and sorted by desc logFC.constDiff 
# (i.e. higher logFC.constDiff = higher expression in H898)
set12 <- 
  cd_vs_wiq.annot.plant %>% 
  filter(set == 1 | set == 2) %>% 
  arrange(desc(logFC.constDiff)) %>%
  dplyr::select(cds, logFC.constDiff, set)

table(set12$set)
#    1    2 
# 1633 1266 

# Reformat into a named numeric vector (nfv)
set12.nnv <- set12$logFC.constDiff
names(set12.nnv) <- set12$cds

length(set12.nnv)
# [1] 2899

table(upRegGenes(set12.nnv))
# FALSE  TRUE 
#  1266  1633 

set1.GOdata <- 
  new("topGOdata",
      description = "Set 1 (BP)",
      ontology = "BP",
      allGenes = set12.nnv,
      geneSel = upRegGenes,
      annot = annFUN.gene2GO,
      gene2GO = cds_2_goid.ncv,
      nodeSize = 10)

set1.weight01.fisher <- 
  runTest(set1.GOdata,
          algorithm = "weight01",
          statistic = "fisher")
# Description: Set 1 (BP) 
# Ontology: BP 
# 'weight01' algorithm with the 'fisher' test
# 204 GO terms scored: 2 terms with p < 0.01
# Annotation data:
#     Annotated genes: 735 
#     Significant genes: 387 
#     Min. no. of genes annotated to a GO: 10 
#     Nontrivial nodes: 204 

set1.ge <- 
  GenTable(set1.GOdata, 
           p.value = set1.weight01.fisher,
           topNodes = 100, 
           numChar = 1000)

set1.ge$p.value <- as.numeric(set1.ge$p.value)
set1.ge <- set1.ge[set1.ge$p.value <= 0.01,]
str(set1.ge)
# 'data.frame':	2 obs. of  6 variables:

set1.ge$set <- 1

# We are now comparing what terms are enriched in Set 2 against Set 1
# Subset data and sorted by logFC.constDiff 
# (i.e. higher logFC.constDiff = higher expression in H898)
table(downRegGenes(set12.nnv))
# FALSE  TRUE 
#  1633  1266 

set2.GOdata <- 
  new("topGOdata",
      description = "Set 2 (BP)",
      ontology = "BP",
      allGenes = set12.nnv,
      geneSel = downRegGenes,
      annot = annFUN.gene2GO,
      gene2GO = cds_2_goid.ncv,
      nodeSize = 10)

set2.weight01.fisher <- 
  runTest(set2.GOdata,
          algorithm = "weight01",
          statistic = "fisher")
# Description: Set 2 (BP) 
# Ontology: BP 
# 'weight01' algorithm with the 'fisher' test
# 204 GO terms scored: 2 terms with p < 0.01
# Annotation data:
#     Annotated genes: 735 
#     Significant genes: 348 
#     Min. no. of genes annotated to a GO: 10 
#     Nontrivial nodes: 203 

set2.ge <- 
  GenTable(set2.GOdata, 
           p.value = set2.weight01.fisher,
           topNodes = 100, 
           numChar = 1000)

set2.ge$p.value <- as.numeric(set2.ge$p.value)
set2.ge <- set2.ge[set2.ge$p.value <= 0.01,]
str(set2.ge)
# 'data.frame':	2 obs. of  6 variables:

set2.ge$set <- 2


# Repeat analysis for set3 against set4
set34 <-
  cd_vs_wiq.annot.plant %>%
  filter(set == 3 | set == 4) %>%
  arrange(desc(logFC.weevilInd_Q903)) %>%
  dplyr::select(cds, logFC.weevilInd_Q903, set)

table(set34$set)
#    3    4 
# 1329 2513 

# Reformat into a named numeric vector (nfv)
set34.nnv <- set34$logFC.weevilInd_Q903
names(set34.nnv) <- set34$cds

length(set34.nnv)
# [1] 3842

table(upRegGenes(set34.nnv))
# FALSE  TRUE 
#  2513  1329 

set3.GOdata <- 
  new("topGOdata",
      description = "Set 3 (BP)",
      ontology = "BP",
      allGenes = set34.nnv,
      geneSel = upRegGenes,
      annot = annFUN.gene2GO,
      gene2GO = cds_2_goid.ncv,
      nodeSize = 10)

set3.weight01.fisher <- 
  runTest(set3.GOdata,
          algorithm = "weight01",
          statistic = "fisher")
# Description: Set 3 (BP) 
# Ontology: BP 
# 'weight01' algorithm with the 'fisher' test
# 297 GO terms scored: 16 terms with p < 0.01
# Annotation data:
#     Annotated genes: 1095 
#     Significant genes: 357 
#     Min. no. of genes annotated to a GO: 10 
#     Nontrivial nodes: 245 

set3.ge <- 
  GenTable(set3.GOdata, 
           p.value = set3.weight01.fisher,
           topNodes = 100, 
           numChar = 1000)

set3.ge$p.value <- as.numeric(set3.ge$p.value)
set3.ge <- set3.ge[set3.ge$p.value <= 0.01,]
str(set3.ge)
# 'data.frame':	16 obs. of  6 variables:

set3.ge$set <- 3


# Set 4
table(downRegGenes(set34.nnv))
# FALSE  TRUE 
#  1329  2513 

set4.GOdata <- 
  new("topGOdata",
      description = "Set 4 (BP)",
      ontology = "BP",
      allGenes = set34.nnv,
      geneSel = downRegGenes,
      annot = annFUN.gene2GO,
      gene2GO = cds_2_goid.ncv,
      nodeSize = 10)

set4.weight01.fisher <- 
  runTest(set4.GOdata,
          algorithm = "weight01",
          statistic = "fisher")
# Description: Set 4 (BP) 
# Ontology: BP 
# 'weight01' algorithm with the 'fisher' test
# 297 GO terms scored: 21 terms with p < 0.01
# Annotation data:
#     Annotated genes: 1095 
#     Significant genes: 738 
#     Min. no. of genes annotated to a GO: 10 
#     Nontrivial nodes: 297 

set4.ge <- 
  GenTable(set4.GOdata, 
           p.value = set4.weight01.fisher,
           topNodes = 100, 
           numChar = 1000)

set4.ge$p.value <- as.numeric(set4.ge$p.value)
set4.ge <- set4.ge[set4.ge$p.value <= 0.01,]
str(set4.ge)
# 'data.frame':	21 obs. of  6 variables:

set4.ge$set <- 4

enrichedGO <- rbind(set1.ge, set2.ge, set3.ge, set4.ge)

enrichedGO <- enrichedGO %>% dplyr::select(-Annotated, -Significant, -Expected)

write.table(enrichedGO, "results/enriched_go_terms.txt",
            sep = "\t", quote = FALSE, 
            col.names = TRUE, row.names = FALSE)
