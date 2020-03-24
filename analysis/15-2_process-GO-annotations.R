library(dplyr)
library(purrr)
library(stringr)
library(tidyr)


########################
### TAIR annotations ###
########################

# BLAST results output on tab-delimited format (i.e. run with -outfmt '6' 
# option on comamnd line BLAST).  We ran with 5 max target hits returned 
# per query sequences.

tair <-
  read.delim("results/all-DE.blastpTAIR.txt", header = TRUE, row.names = NULL,
             colClasses = c(rep("character", 2), rep("numeric", 14), rep("character", 5)))

tair <- tair %>% select(qseqid, sseqid, evalue, salltitles)

str(tair)
# 'data.frame':	34636 obs. of  4 variables:


# Extract only the top hit for downstream analysis
topTairHit <- tair %>% group_by(qseqid) %>% nest()

topTairHit$data <- map(topTairHit$data, function(x) {
  # strip version from accession number
  x$sseqid <- str_replace(x$sseqid, "\\.\\d", "")
  x[1,]
})

topTairHit <- topTairHit %>% unnest(c(data))

topTairHit <- topTairHit %>% select(qseqid, sseqid, evalue, salltitles)

str(topTairHit)
# Classes ‘grouped_df’, ‘tbl_df’, ‘tbl’ and 'data.frame':	7012 obs. of  4 variables:


# Collect all unique UniProt IDs from top BLAST hit
tIDs <- unique(sort(topTairHit$sseqid))
length(tIDs)
# [1] 3015


# Read UniProt to KEGG id mapping, downloaded from arabidopsis.org
athGoSlim <- read.csv("data/ATH_GO_GOSLIM.txt", header = FALSE,
                      sep = "\t", stringsAsFactors = FALSE, comment.char = "!")

# Columns in file
# 1. locus name: standard AGI convention name
# 2. TAIR accession:the unique identifier for an object in the TAIR database- 
#   the object type is the prefix, followed by a unique accession number(e.g. gene:12345).  
# 3. object name : the name of the object (gene, protein, locus) being annotated.
# 4. relationship type: the relationship between the annotated object and the GO term
# 5. GO term: the actual string of letters corresponding to the GO ID
# 6. GO ID: the unique identifier for a GO term.  
# 7. TAIR Keyword ID: the unique identifier for a keyword in the TAIR database.
# 8.  Aspect: F=molecular function, C=cellular component, P=biological 13process. 
# 9. GOslim term: high level GO term helps in functional categorization.
# 10. Evidence code: three letter code for evidence types (see: http://www.geneontology.org/GO.evidence.html).
# 11. Evidence description: the analysis that was done to support the annotation
# 12. Evidence with: supporting evidence for IGI, IPI, IC, IEA and ISS annotations
# 13. Reference: Either a TAIR accession for a reference (reference table: reference_id) or reference from PubMed (e.g. PMID:1234).  
# 14. Annotator: TAIR, TIGR, GOC (GO Consortium), UniProt, IntAct or a TAIR community member
# 15. Date annotated: date the annotation was made.

athGoSlim <- athGoSlim %>% select(1,2,3,4,5,6,8,9)

colnames(athGoSlim) <- c("geneId", "locus", "accession", "relationship_type", 
                         "GO_term", "GO_ID", "aspect", "GOslim_term")

# Only keep biological process annotations
athGoSlim <- athGoSlim %>% filter(athGoSlim$aspect == "P")
str(athGoSlim)
# 'data.frame':	148546 obs. of  8 variables:


# Only keep those from BLAST results
athGoSlim <- subset(athGoSlim, athGoSlim$geneId %in% tIDs)
str(athGoSlim)
# 'data.frame':	21547 obs. of  8 variables:


# Remove terms that are uninformative
athGoSlim <- athGoSlim %>% filter(GOslim_term != "biosynthetic process")
athGoSlim <- athGoSlim %>% filter(GOslim_term != "other cellular processes")
athGoSlim <- athGoSlim %>% filter(GOslim_term != "other metabolic processes")
athGoSlim <- athGoSlim %>% filter(GOslim_term != "unknown biological processes")


# Only keep unique geneId, GOslim_term pairs
athGoSlim <- athGoSlim %>% select(geneId, GOslim_term) %>% unique()

str(athGoSlim)
# 'data.frame':	7824 obs. of  2 variables:

topTairGo <- merge(topTairHit, athGoSlim, by.x = "sseqid", by.y = "geneId")

topTairGo <- topTairGo %>% select(qseqid, sseqid, evalue, salltitles, GOslim_term)

write.table(topTairGo, "results/all-DE.blastpTAIR.topHit.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
