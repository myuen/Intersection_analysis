library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

source("analysis/helper01_getLineage.R")

# BLAST results output on tab-delimited format (i.e. run with -outfmt '6' 
# option on comamnd line BLAST).  We ran with 5 max target hits returned 
# per query sequences.

blast <-
  read.delim("results/all-DE.blastpUniProt.txt", 
             header = TRUE, row.names = NULL,
             colClasses = c(rep("character", 2), 
                            rep("numeric", 14), rep("character", 5)))

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


# Let's extract all unique taxonomy ID from all top blast hit
uniq.TaxIDs <- unique(sort(topBlastHit$tax_id))
length(uniq.TaxIDs)
# [1] 510


uniq.TaxIDs <- getLineage(uniq.TaxIDs) %>% select(old_tax_id, taxonomy)
colnames(uniq.TaxIDs)[1] <- "tax_id"
str(uniq.TaxIDs)
# 'data.frame':	510 obs. of  2 variables:

 
topBlastHit <- left_join(topBlastHit, uniq.TaxIDs) 
str(topBlastHit)
# Classes ‘grouped_df’, ‘tbl_df’, ‘tbl’ and 'data.frame':	9714 obs. of  6 variables:


write.table(topBlastHit, "results/all-DE.blastpUniProt.topHit.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
