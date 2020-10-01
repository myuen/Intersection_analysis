library(dplyr)

# Download taxonomy lineage file from NCBI.
# It's a database dump delimited by a pipe "|" and each 
# value is delimited by tab on both ends.
# ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

# Read merged id file to convert old tax id to new tax id
mergedTaxIDs <- read.csv("data/ncbi-taxonomy-dump/merged.dmp",
                         stringsAsFactors = FALSE, 
                         sep = "\t", header = FALSE)

mergedTaxIDs <- mergedTaxIDs %>% select(1,3)

colnames(mergedTaxIDs) <- c("old_tax_id", "new_tax_id")


# Read ranked lineage from NCBI tax id
rl <- read.csv("data/ncbi-taxonomy-dump/rankedlineage.dmp",
               stringsAsFactors = FALSE, 
               sep = "\t", header = FALSE)

rl <- rl %>% select(1, 3, 5, 7, 9, 11, 13, 15, 17, 19) 

colnames(rl) <- c("tax_id", "tax_name", "species", "genus", 
                  "family", "order", "class", "phylum", 
                  "kingdom", "superkingdom")


getLineage <- function(x) {
  # populate data frame with x in both col
  df <- data.frame("old_tax_id" = x,
                   "new_tax_id" = x)

  # Check for retired or merged Tax IDs 
  retiredTaxIds <- 
    mergedTaxIDs %>% filter(old_tax_id %in% x)
  
  # Remove rows with retired Tax IDs from df
  df <- df[!df$old_tax_id %in% retiredTaxIds$old_tax_id, ]
  
  df <- rbind(df, retiredTaxIds)
  
  # Only keep those that are found in our BLAST results
  df <- merge(df, rl, by.x = "new_tax_id", by.y = "tax_id", all.x = TRUE)

  # Set taxonomy if it belongs to Bacteria or Viruses.
  # Drill deeper if it belongs Eukaryota and report if it 
  # belongs to Fungi or Viridiplantae.  Report the class for 
  # all other case.
  
  df$Taxonomy <-
    case_when(
      df$superkingdom == "Bacteria" | df$superkingdom == "Viruses" ~ 
        df$superkingdom, 
      df$kingdom == "Fungi" | df$kingdom == "Viridiplantae" ~
        df$kingdom,
      df$phylum != "" ~ df$phylum
    )
  
  # For those don't have class name, replace with taxonomy name
  df[is.na(df$Taxonomy), "Taxonomy"] <- 
    df[is.na(df$Taxonomy), "tax_name"]
  
  return(df)
}
