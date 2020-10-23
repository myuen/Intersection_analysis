library(dplyr)
library(purrr)
library(stringr)

# Get UniProt to GO term mapping file from EMBL
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

# Extract 1,2,7 column from idmapping_selected.tab.gz with GO terms on command line
# 1. UniProtKB-AC
# 2. UniProtKB-ID
# 7. GO

uId_2_goId <- 
  read.csv('data/uniqueUniProtId2goId.txt', sep = "\t", 
           header = FALSE, stringsAsFactors = FALSE, 
           colClasses = rep("character", 3))

colnames(uId_2_goId) <- c("UniProtKB.AC", "UniProtKB.ID", "GO")


# Function takes a single or a vector UniProtKB.AC
getGOID <- function(x) {
  
  # Remove unused column
  uId_2_goId <- uId_2_goId %>% 
    dplyr::select(-UniProtKB.ID)
  
  uId_2_goId.subset <- uId_2_goId[uId_2_goId$UniProtKB.AC %in% x,]
  
  return (uId_2_goId.subset)
}
