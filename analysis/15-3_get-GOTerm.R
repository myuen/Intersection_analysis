library(GO.db)
library(dplyr)

# Read UniProt ID to GO ID mapping
uId_2_goId <- read.delim("results/uId_2_goId.txt", sep = "\t", 
                        header = TRUE, stringsAsFactors = FALSE)
str(uId_2_goId)
# 'data.frame':	10709 obs. of  2 variables:

go.db <- AnnotationDbi::select(GO.db, keys(GO.db, "GOID"), 
                               c("TERM", "ONTOLOGY"))

uId_2_goTerm <- left_join(uId_2_goId, go.db)

# Only keep Biological Process terms
uId_2_goTerm <- uId_2_goTerm %>% filter(ONTOLOGY == "BP") %>%
  select(-ONTOLOGY)
str(uId_2_goTerm)
# 'data.frame':	2446 obs. of  4 variables:

write.table(uId_2_goTerm, "results/uId_2_goTerm.txt", 
            col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)
