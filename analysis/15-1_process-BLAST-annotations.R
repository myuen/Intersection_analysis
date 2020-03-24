library(dplyr)
library(purrr)
library(stringr)
library(tidyr)


#########################
### BLAST annotations ###
#########################

# BLAST results output on tab-delimited format (i.e. run with -outfmt '6' 
# option on comamnd line BLAST).  We ran with 5 max target hits returned 
# per query sequences.

blast <-
  read.delim("results/all-DE.blastpUniProt.txt", header = TRUE, row.names = NULL,
             colClasses = c(rep("character", 2), rep("numeric", 14), rep("character", 5)))

blast <- blast %>% select(qseqid, sseqid, evalue, salltitles)

str(blast)
# 'data.frame':	50830 obs. of  4 variables:


# Extract only the top hit for downstream analysis
topBlastHit <- blast %>% group_by(qseqid) %>% nest()

topBlastHit$data <- map(topBlastHit$data, function(x) {
  # Extract UniProt ID from sseqid
  x[1, "sseqid"] <- str_split(x[1 ,"sseqid"], "\\|") %>% unlist() %>% nth(2)
  # Extract tax id from description
  x[1, "tax_id"] <- str_split(x[1 ,"salltitles"], "=") %>% unlist %>% nth(3) %>% 
    str_replace(" GN", "") %>% str_replace(" PE", "")
  x[1,]
})

topBlastHit <- topBlastHit %>% unnest(c(data))

topBlastHit$tax_id <- as.integer(topBlastHit$tax_id)

str(topBlastHit)
# Classes ‘grouped_df’, ‘tbl_df’, ‘tbl’ and 'data.frame':	9714 obs. of  5 variables:


# Let's extract all unique organism name
uniqTaxIDs <- data.frame("tax_id" = as.integer(unique(sort(topBlastHit$tax_id))))
str(uniqTaxIDs)
# 'data.frame':	510 obs. of  1 variable:


# Read taxonomy lineage file downloaded from NCBI
# It's a database dump delimited by a pipe "|" and each 
# value is delimited by tab on both ends.
# ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

# Read merged id file to convert old tax id to new tax id
mergedTaxIDs <- read.csv("data/ncbi-taxonomy-dump/merged.dmp",
                         stringsAsFactors = FALSE, 
                         sep = "\t", header = FALSE)

mergedTaxIDs <- mergedTaxIDs %>% select(1,3)

colnames(mergedTaxIDs) <- c("old_tax_id", "new_tax_id")

table(uniqTaxIDs$tax_id %in% mergedTaxIDs$old_tax_id)
# FALSE  TRUE 
#   508     2 

uniqTaxIDs <- merge(uniqTaxIDs, mergedTaxIDs, 
                    by.x = "tax_id", by.y = "old_tax_id", all.x = TRUE)

uniqTaxIDs[is.na(uniqTaxIDs$new_tax_id), "new_tax_id"] <- 
  uniqTaxIDs[is.na(uniqTaxIDs$new_tax_id), "tax_id"]


# Read ranked lineage from NCBI tax id
rl <- read.csv("data/ncbi-taxonomy-dump/rankedlineage.dmp",
               stringsAsFactors = FALSE, 
               sep = "\t", header = FALSE)

rl <- rl %>% select(1, 3, 5, 7, 9, 11, 13, 15, 17, 19) 

colnames(rl) <- c("tax_id", "tax_name", "species", "genus", 
                  "family", "order", "class", "phylum", 
                  "kingdom", "superkingdom")


# Only keep those that are found in our BLAST results
uniqTaxIDs <- merge(uniqTaxIDs, rl, 
                    by.x = "new_tax_id", by.y = "tax_id", all.x = TRUE)

table(is.na(uniqTaxIDs$tax_name))
# FALSE 
#   510 

# Set taxonomy if it belongs to Bacteria or Viruses.
# Drill deeper if it belongs Eukaryota and report if it 
# belongs to Fungi or Viridiplantae.  Report the class for 
# all other case.

uniqTaxIDs$taxonomy <-
  case_when(
    uniqTaxIDs$superkingdom == "Bacteria" | uniqTaxIDs$superkingdom == "Viruses" ~ 
      uniqTaxIDs$superkingdom, 
    uniqTaxIDs$kingdom == "Fungi" | uniqTaxIDs$kingdom == "Viridiplantae" ~
      uniqTaxIDs$kingdom,
    uniqTaxIDs$phylum != "" ~ uniqTaxIDs$phylum
)


# For those don't have class name, replace with taxonomy name
uniqTaxIDs[is.na(uniqTaxIDs$taxonomy), "taxonomy"] <- 
  uniqTaxIDs[is.na(uniqTaxIDs$taxonomy), "tax_name"]
str(uniqTaxIDs)
# 'data.frame':	510 obs. of  12 variables:
 

topBlastHit <- left_join(topBlastHit, uniqTaxIDs) %>% 
  select(qseqid, sseqid, evalue, salltitles, tax_id, taxonomy)

write.table(topBlastHit, "results/all-DE.blastpUniProt.topHit.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
