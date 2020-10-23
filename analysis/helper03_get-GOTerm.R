library(GO.db)
library(dplyr)

go.db <- AnnotationDbi::select(GO.db, keys(GO.db, "GOID"), 
                               c("TERM", "ONTOLOGY"))

getGOTerm <- function(x){
  uId_2_goTerm <- subset(go.db, go.db$GOID %in% x)
  
  # Only keep Biological Process terms
  uId_2_goTerm <- uId_2_goTerm %>%
    filter(ONTOLOGY == "BP") %>%
    dplyr::select(-ONTOLOGY)
}