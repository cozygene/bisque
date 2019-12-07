context("Reference-based decomposition")
library(Biobase)

test_that("Catches expression input not in expressionset", {
  bulk.counts <- matrix(0,nrow=2,ncol=2)
  sc.counts <- matrix(0,nrow=2,ncol=2)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.counts, sc.eset))
})

test_that("Catches expressionset missing given labels", {
  bulk.counts <- matrix(0,nrow=2,ncol=2)
  sc.counts <- matrix(0,nrow=2,ncol=2)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.pheno <- data.frame(SubjectName=c('a', 'b'), CellType=c('a', 'b'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "CellType"),
                        row.names=c("SubjectName", "CellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
  sc.pheno <- data.frame(subjectName=c('a', 'b'), cellType=c('a', 'b'))
  sc.meta <- data.frame(labelDescription=c("subjectName", "cellType"),
                        row.names=c("subjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
})

test_that("Catches single-cell data with only one or two subjects", {
  bulk.counts <- matrix(1:4,nrow=2,ncol=2)
  sc.counts <- matrix(0,nrow=2,ncol=2)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  # One Subject
  sc.pheno <- data.frame(SubjectName=c('a', 'a'), cellType=c('a', 'b'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset,
                                                      use.overlap=FALSE))
  # Two subjects
  sc.counts <- matrix(1:8,nrow=2,ncol=4)
  sc.pheno <- data.frame(SubjectName=c('a', 'a', 'b', 'b'),
                         cellType=c('a', 'b', 'a', 'b'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_warning(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset,
                                                        use.overlap=FALSE),
                 regexp="Only two individuals detected in single-cell data")
})

test_that("Catches single-cell data with only one cell type label", {
  bulk.counts <- matrix(0,nrow=2,ncol=2)
  sc.counts <- matrix(0,nrow=2,ncol=2)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.pheno <- data.frame(SubjectName=c('a', 'b'), cellType=c('a', 'a'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
})

test_that("Catches no overlapping samples in overlap mode", {
  bulk.counts <- matrix(0,nrow=2,ncol=2)
  sc.counts <- matrix(0,nrow=2,ncol=2)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.pheno <- data.frame(SubjectName=c('a', 'b'), cellType=c('a', 'b'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
  colnames(bulk.counts) <- c("a","b")
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
})

test_that("Catches no overlapping genes between bulk and single-cell and markers", {
  markers <- c("a", "b")
  bulk.counts <- matrix(0,nrow=2,ncol=2)
  sc.counts <- matrix(0,nrow=2,ncol=2)
  colnames(bulk.counts) <- c("a", "c")
  rownames(bulk.counts) <- c("a", "b")
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.pheno <- data.frame(SubjectName=c('a', 'b'), cellType=c('a', 'b'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
  rownames(bulk.counts) <- c("1", "2")
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, markers))
  
})

test_that("Catches data with no expressed/zero variance genes", {
  bulk.counts <- matrix(1,nrow=2,ncol=2)
  sc.counts <- matrix(1,nrow=2,ncol=2)
  colnames(bulk.counts) <- c("a", "c")
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.pheno <- data.frame(SubjectName=c('a', 'b'), cellType=c('a', 'b'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
})

test_that("Catches data where zero variance in training data", {
  # all genes good
  bulk.counts <- matrix(1:9,nrow=3,ncol=3)
  sc.counts <- matrix(1:9,nrow=3,ncol=3)
  colnames(bulk.counts) <- c("a", "b", "c")
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.pheno <- data.frame(SubjectName=c('b', 'c', 'd'), cellType=c('a', 'b', 'c'))
  sc.meta <- data.frame(labelDescription=c("SubjectName", "cellType"),
                        row.names=c("SubjectName", "cellType"))
  sc.pdata <- new("AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.meta)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  # all genes good, no overlap
  expect_equal(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset, use.overlap=F)$genes.used,
               c("1","3"))
  # Bulk has some zero variance
  bulk.counts[1,c("b","c")] <- c(0,0)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  expect_equal(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset)$genes.used, c("3"))
  # Single cell has some zero variance
  bulk.counts <- matrix(1:9,nrow=3,ncol=3)
  colnames(bulk.counts) <- c("a", "b", "c")
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.counts[1,c(1,2)] <- c(0,0)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_equal(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset)$genes.used, c("3"))
  # All single cell training data has zero variance
  sc.counts[1,c(1,2)] <- c(1,1)
  sc.counts[2,c(1,2)] <- c(1,1)
  sc.counts[3,c(1,2)] <- c(1,1)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_error(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset))
  # Both have zero variance and undefined coefficient
  sc.counts <- matrix(1:9,nrow=3,ncol=3)
  sc.counts[1,c(1,2)] <- c(0,0)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  bulk.counts[1,c("b","c")] <- c(0,0)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  expect_equal(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset)$genes.used, c("2", "3"))
  # Both have zero variance but defined coefficient
  bulk.counts[1,c("b","c")] <- c(2,2)
  bulk.counts[c(2,3), "c"] <- c(6,5)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  sc.counts[1,c(1,2)] <- c(1,1)
  sc.counts[c(2,3), 2] <- c(3,2)
  sc.eset <- Biobase::ExpressionSet(assayData = sc.counts, phenoData = sc.pdata)
  expect_equal(BisqueRNA::ReferenceBasedDecomposition(bulk.eset, sc.eset)$genes.used, c("1","2", "3"))
  })
