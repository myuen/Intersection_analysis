library(dplyr)
library(purrr)
library(stringr)
library(tidyr)


# This will take some time to load
source("analysis/helper01_getLineage.R")
source("analysis/helper02_get-GOID.R")

topBlastHit <- read.delim("results/all-DE.blastpUniProt.topHit.txt",
                          sep = "\t", stringsAsFactors = FALSE)

str(topBlastHit)
# 'data.frame':	9714 obs. of  5 variables:


# 1) Add result from set analysis
cd_vs_wiq <- read.table("results/cd_vs_wiq.txt", 
                        header = TRUE, 
                        stringsAsFactors = FALSE)

str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  6 variables:


# Only keeping BLAST results for contigs diff. exp. in constDiff and
# weevilInd_Q903 for downstream analysis
cd_vs_wiq.annot <- left_join(cd_vs_wiq, topBlastHit)
# Joining, by = "cds"


# 2) Extract taxonomy ID from all top blast hit and add lineage
uniq.TaxIDs <- unique(sort(cd_vs_wiq.annot$tax_id))
length(uniq.TaxIDs)
# [1] 381

uniq.TaxIDs <- getLineage(uniq.TaxIDs) %>% 
  select(old_tax_id, Taxonomy)

colnames(uniq.TaxIDs)[1] <- "tax_id"

str(uniq.TaxIDs)
# 'data.frame':	381 obs. of  2 variables:

cd_vs_wiq.annot <- left_join(cd_vs_wiq.annot, uniq.TaxIDs) 
# Joining, by = "tax_id"


# 3) Add GOID
# Extract unique UniProt ID
uniq.uId <- unique(sort(cd_vs_wiq.annot$UniProtKB.AC))
length(uniq.uId)
# [1] 4369

uId_2_goId <- getGOID(uniq.uId)

str(uId_2_goId)
# 'data.frame':	3213 obs. of  2 variables:

cd_vs_wiq.annot <- left_join(cd_vs_wiq.annot, uId_2_goId)
# Joining, by = "UniProtKB.AC"

str(cd_vs_wiq.annot)
# 'data.frame':	9298 obs. of  12 variables:

# Table S2
write.table(cd_vs_wiq.annot, "results/cd_vs_wiq.annot.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
