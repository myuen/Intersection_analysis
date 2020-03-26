library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)


# Read top BLAST hit
topBlastHit <- read.delim("results/all-DE.blastpUniProt.topHit.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)
topBlastHit <- topBlastHit %>% select(qseqid, sseqid, taxonomy)
str(topBlastHit)
# 'data.frame':	9714 obs. of  3 variables:


# Extract unique UniProt ID
unique.uId <- unique(sort(topBlastHit$sseqid))
length(unique.uId)
# [1] 5621


# Get UniProt to GO term mapping file from EMBL
# ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

# Extract 1,2,7 column from idmapping_selected.tab.gz with GO terms on command line
# 1. UniProtKB-AC
# 2. UniProtKB-ID
# 7. GO

# Extract subject IDs (i.e. UniProt Accession) from BLAST results.
# grep lines with these IDs on command line
uId_2_goId <- read.csv('data/uniqueUniProtId2goId.txt', sep = "\t", header = FALSE,
                       stringsAsFactors = FALSE, 
                       colClasses = rep("character", 3))
colnames(uId_2_goId) <- c("UniProtKB.AC", "UniProtKB.ID", "GO")

str(uId_2_goId)
# 'data.frame':	7450 obs. of  3 variables:

# We are only keeping those line that the BLAST hit against
uId_2_goId <- subset(uId_2_goId, uId_2_goId$UniProtKB.AC %in% unique.uId) %>%
  select(-UniProtKB.ID)
str(uId_2_goId)
# 'data.frame':	4144 obs. of  2 variables:


# Split the GOID delimited by ';' into separate rows
uId_2_goId <- map_dfr(uId_2_goId$UniProtKB.AC, function(x){
  goIds <- uId_2_goId %>% filter(UniProtKB.AC == x) %>% str_split(";")
  data.frame("UniProtKB.AC" = x, "GOID" = str_trim(unlist(goIds[[2]])))
}
)

write.table(uId_2_goId, "results/uId_2_goId.txt", 
            col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)
