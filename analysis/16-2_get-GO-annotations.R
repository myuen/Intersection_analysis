library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

source("analysis/helper02_get-GOID.R")
source("analysis/helper03_get-GOTerm.R")


# Read top BLAST hit
topBlastHit <- read.delim("results/all-DE.blastpUniProt.topHit.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)

topBlastHit <- topBlastHit %>% dplyr::select(cds, UniProtKB.AC, Taxonomy)

str(topBlastHit)
# 'data.frame':	9298 obs. of  3 variables:


# Extract unique UniProt ID
uniq.uId <- unique(sort(topBlastHit$UniProtKB.AC))
length(uniq.uId)
# [1] 4369


# We are only keeping those line that the BLAST hit against
uId_2_goId <- getGOID(uniq.uId)
str(uId_2_goId)
# 'data.frame':	8206 obs. of  2 variables:

goId_2_goTerm <- getGOTerm(unique(uId_2_goId$GOID))
str(goId_2_goTerm)
# 'data.frame':	535 obs. of  2 variables:

uId_2_goTerm <- inner_join(uId_2_goId, goId_2_goTerm)
# Joining, by = "GOID"

str(uId_2_goTerm)
# 'data.frame':	1817 obs. of  3 variables:


blastHitGoTerms <- inner_join(topBlastHit, uId_2_goTerm)
# Joining, by = "UniProtKB.AC"

str(blastHitGoTerms)
# 'data.frame':	2824 obs. of  5 variables:

write.table(blastHitGoTerms, "results/blastHitGoTerms.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
