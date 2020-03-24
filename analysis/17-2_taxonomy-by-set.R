library(dplyr)

# Read top BLAST hit and taxonomy lineage
topBlastHit <- read.delim("results/all-DE.blastpUniProt.topHit.txt", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)
topBlastHit <- topBlastHit %>% select(qseqid, taxonomy)
str(topBlastHit)
# 'data.frame':	9714 obs. of  2 variables:


# Read set analysis
cd_vs_wiq <- read.delim("results/cd_vs_wiq.txt", sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE) 
cd_vs_wiq <- cd_vs_wiq %>% select(cds, set)
str(cd_vs_wiq)
# 'data.frame':	9298 obs. of  2 variables:


table(cd_vs_wiq$set)
#    1    2    3    4    5    6    7    8 
# 2065 1514 1855 2975   85   69  104  631 


cd_vs_wiq.taxonomy <- merge(cd_vs_wiq, topBlastHit, by.x = "cds", by.y = "qseqid")
str(cd_vs_wiq.taxonomy)
# 'data.frame':	7928 obs. of  3 variables:
# Some sequence don't have annotation, hence size not the same as cd_vs_wiq


# Summarize taxonomy count
taxCount <- cd_vs_wiq.taxonomy %>% group_by(set, taxonomy) %>% 
  summarise(count = n())

colnames(taxCount) <- c("set", "Taxonomy", "Count")

write.table(taxCount, "results/taxCountBySet.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE)
