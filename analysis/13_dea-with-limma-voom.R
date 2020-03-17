library(dplyr)
library(edgeR)
library(purrr)
library(stringr)
library(testthat) # facilitate tests that will catch changes on re-analysis


### Differential Expression Analysis on Sitka Spruce Weevil 
### Experiment with limma + voom

# abs(logFC) log fold change cut-off.  Anything greater than (-1 x lfc) and less 
# than lfc will be deemed biological insignificant
lfcCutoff <- 2

# p-value cut-off.  Anything > pCutoff will be deemed statistically insignificant.
pCutoff <- 0.05

#' Load counts from Sailfish
x <- read.delim("data/consolidated-Sailfish-results.rRNA-Filtered.lowExp-Filtered.txt",
                colClasses = c("character", rep("numeric", 24)), header = TRUE, row.names = 1)
str(x, list.len = 8) # 38711 obs. of  24 variables

#' Load experimental design
expDes <- read.delim("data/White_Pine_Weevil_exp_design.tsv",
                     stringsAsFactors = FALSE)
expDes <-mutate(expDes,
                gType = factor(gType, levels = c("Q903", "H898")),
                txCode = factor(txCode, levels = c('C', 'W', 'G')),
                tx = factor(tx, levels = c("Control", "Wound", "Gallery")), 
                grp = factor(grp, levels = paste(levels(gType),
                                                 rep(levels(tx), each = 2), sep = ".")))
str(expDes) # 'data.frame':  24 obs. of  6 variables


#' Load counts into DGEList object from edgeR package.
y <- DGEList(counts = x, group = expDes$grp)

#' TMM Normalization by Depth
y <- calcNormFactors(y)

#' Write CPM output for heatmap and contig plotting
write.table(cpm(y), file = "results/normalized_cpm.22Feb.txt", quote = FALSE,
            row.names = TRUE, col.names = TRUE)

#' make model matrix
modMat <- model.matrix(~ gType * tx, expDes)

#' It's hard to believe, but the default names for columns associated with 
#' interaction terms will create fatal errors in makeContrasts below; prevent 
#' that, and a warning about the intercept, by modifying these column names 
#' here; see 90_limma-model-term-name-fiasco.r and .md for more
colnames(modMat)
colnames(modMat) <- gsub(":", "_", colnames(modMat))
colnames(modMat) <- gsub("[()]", "", colnames(modMat))
colnames(modMat)

#' voom transformation
#+ voom-plot
v <- voom(y, modMat, plot = TRUE) # take a couple moments

#' Linear modelling and forming contrasts
fit <- lmFit(v, modMat)

cont_matrix <-
  makeContrasts(
    # A - Between genotype comparisons
    # A1 - constitutive difference (H898 control - Q903 control)
    constDiff = gTypeH898,
    # B - Within genotype comparisons
    # B1 - Wound Response (Wound - Control)
    woundResp_Q903 = txWound,
    woundResp_H898 = txWound + gTypeH898_txWound,
    # B1 - Weevil Induced (Gallery - Wound)
    weevilInd_Q903 = txGallery - txWound,
    weevilInd_H898 = txGallery - txWound +
      gTypeH898_txGallery - gTypeH898_txWound,
    # B2 - Weevil Control (Gallery - Control)
    weevilCtrl_Q903 = txGallery,
    weevilCtrl_H898 = txGallery + gTypeH898_txGallery,
    levels = modMat)

fit2 <- contrasts.fit(fit, cont_matrix)

fit3 <- eBayes(fit2)

summary(decideTests(fit3, p.value = pCutoff, lfc = lfcCutoff))
#    constDiff woundResp_Q903 woundResp_H898 weevilInd_Q903
# -1      2214              1              0           3710
# 0      34243          38710          38711          32992
# 1       2254              0              0           2009

#    weevilInd_H898 weevilCtrl_Q903 weevilCtrl_H898
# -1              0            3878               1
# 0           38711           31395           38708
# 1               0            3438               2


focus_terms <- colnames(cont_matrix)

names(focus_terms) <- focus_terms

#' to get individual t statistics, etc., must use topTable() on each coef
#' separately
statInf_focus_terms <-
  map_df(focus_terms, function(x) {
    tmp <- topTable(fit3, coef = x, number = Inf, sort.by = "none")
    tmp$focus_term <- x
    tmp$cds <- row.names(tmp) 
    tmp$contig <- str_replace(tmp$cds, ":\\d+\\-\\S+", "")
    return(tmp)
    })
str(statInf_focus_terms)
# 'data.frame':	270977 obs. of  9 variables:

statInf_focus_terms$focus_term <- 
  factor(statInf_focus_terms$focus_term, 
         levels = c("constDiff", 
                    "woundResp_Q903", "woundResp_H898", "weevilInd_Q903",
                    "weevilInd_H898", "weevilCtrl_Q903", "weevilCtrl_H898"))

# Reorganize columns
statInf_focus_terms <- 
  select(statInf_focus_terms, cds, contig, 
         focus_term, logFC, AveExpr, t, P.Value, adj.P.Val, B)

test_that("stat inf on the focus terms has 270977 rows",
          expect_equal(38711 * 7, nrow(statInf_focus_terms)))

write.table(statInf_focus_terms,
            "results/stats-results.long.22Feb.tsv",
            quote = FALSE, sep = "\t", row.names = FALSE)

statInf_focus_terms.wide <- 
  reshape(statInf_focus_terms, direction = "wide",
          timevar = "focus_term", idvar = c("cds"),
          new.row.names = unique(statInf_focus_terms$cds),
          drop = c("AveExpr", "t", "P.Value", "B", "contig"))

colnames(statInf_focus_terms.wide)[1] <- "contig"

statInf_focus_terms.wide$contig <- 
  str_replace(row.names(statInf_focus_terms.wide), ":\\d+\\-\\S+", "")

write.table(statInf_focus_terms.wide,
            "results/stats-results.wide.22Feb.tsv",
            quote = FALSE, sep = "\t", row.names = TRUE)
