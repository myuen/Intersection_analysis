library(edgeR)
library(testthat) # facilitate tests that will catch changes on re-analysis


### Experiment with limma + voom

# Load counts from Sailfish
rawSailfishCounts <- read.delim("data/consolidated-Sailfish-results.rRNA-Filtered.txt",
                                colClasses = c("character", rep("numeric", 24)), row.names = 1)
str(rawSailfishCounts)
# 108118 obs. of  24 variables

test_that("Sailfish data has 108118 rows upon import",
          expect_equal(108118, nrow(rawSailfishCounts)))

test_that("Sailfish data has data for exactly 24 samples",
          expect_equal(24, ncol(rawSailfishCounts)))


# Load counts into DGEList object from edgeR package.
y <- DGEList(counts = rawSailfishCounts)
lib_size <- data.frame(raw = y$samples$lib.size)

# Filtering low expression genes
# We are setting an arbitary threshold and only keeping contigs with more than
# 1 count-per-million (cpm) in at least 2 samples
y <- y[(rowSums(cpm(y) > 1) > 2), ]
test_that("After low expression filter, we have 38197 rows",
          expect_equal(38711, nrow(y)))

## write to file
write.table(y$counts, "data/consolidated-Sailfish-results.rRNA-Filtered.lowExp-Filtered.txt",
            sep = "\t", quote = FALSE)
