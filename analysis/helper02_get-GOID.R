library(dplyr)
library(purrr)
library(stringr)


# Get UniProt to GO term mapping file from EMBL
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

# Extract 1,2,7 column from idmapping_selected.tab.gz with GO terms on command line
# 1. UniProtKB-AC
# 2. UniProtKB-ID
# 7. GO

# Extract subject IDs (i.e. UniProt Accession) from BLAST results.
# grep lines with these IDs on command line
uId_2_goId <- 
  read.csv('data/uniqueUniProtId2goId.txt', sep = "\t", 
           header = FALSE, stringsAsFactors = FALSE, 
           colClasses = rep("character", 3))

colnames(uId_2_goId) <- c("UniProtKB.AC", "UniProtKB.ID", "GO")

getGOID <- function(x) {
  uId_2_goId <- subset(uId_2_goId, uId_2_goId$UniProtKB.AC %in% x) %>%
  dplyr::select(-UniProtKB.ID)

  # Split the GOID delimited by ';' into separate rows
  uId_2_goId <- map_dfr(uId_2_goId$UniProtKB.AC, function(x){
    goIds <- uId_2_goId %>% filter(UniProtKB.AC == x) %>% str_split(";")
    
    data.frame("UniProtKB.AC" = x, "GOID" = str_trim(unlist(goIds[[2]])))
  }
  )
}
