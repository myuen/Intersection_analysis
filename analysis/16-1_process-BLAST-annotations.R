library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

# Process BLAST results. Keeping only tophit and extract IDs from description.


# BLAST results output on tab-delimited format (i.e. run with -outfmt '6' option
# on comamnd line BLAST).  We ran with 5 max target hits returned per query
# sequences.
blast <-
  read.delim("results/all-DE.blastpUniProt.txt", 
             header = TRUE, row.names = NULL,
             colClasses = c(rep("character", 2), 
                            rep("numeric", 14), rep("character", 5)))

blast <- blast %>% 
  select(qseqid, sseqid, evalue, salltitles)

str(blast)
# 'data.frame':	50830 obs. of  4 variables:


# Extract only the top hit for downstream analysis
topBlastHit <- blast %>% 
  group_by(qseqid) %>% 
  nest()

topBlastHit$data <- map(topBlastHit$data, function(x) {
  # Extract UniProt ID from sseqid
  x[1, "sseqid"] <- 
    str_split(x[1 ,"sseqid"], "\\|") %>% 
    unlist() %>% 
    nth(2)
  
  # Extract tax id from description
  x[1, "tax_id"] <- str_split(x[1 ,"salltitles"], "=") %>% 
    unlist %>% 
    nth(3) %>% 
    str_replace(" GN", "") %>% 
    str_replace(" PE", "")

  x[1,]
})

topBlastHit <- topBlastHit %>% 
  unnest(c(data))

topBlastHit$tax_id <- as.integer(topBlastHit$tax_id)

# Renamed column for better readability
colnames(topBlastHit) <- 
  c("cds", "UniProtKB.AC", "evalue", "Description", "tax_id")

str(topBlastHit)
# tibble [9,714 × 5] (S3: grouped_df/tbl_df/tbl/data.frame)

write.table(topBlastHit, "results/all-DE.blastpUniProt.topHit.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
