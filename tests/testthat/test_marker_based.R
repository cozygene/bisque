context("Marker-based decomposition")
library(Biobase)

test_that("Catches input errors", {
  # Expression not in eset
  bulk.counts <- matrix(0,nrow=2,ncol=2)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  markers <- data.frame(gene=c("1","2"), cluster=c("a", "b"))
  expect_error(BisqueRNA::MarkerBasedDecomposition(bulk.counts, markers))
  # Min gene greater than max gene
  expect_error(BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers,
                                                   min_gene=6, max_gene=5))
  # Not enough markers
  expect_error(BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers))
  # No overlapping markers
  markers <- data.frame(gene=c("3","4"), cluster=c("a", "b")) 
  expect_error(BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, min_gene=1))
  # No marker genes after filtering for zero variance
  markers <- data.frame(gene=c("1","2"), cluster=c("a", "b"))
  expect_error(BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, min_gene=1))
  # One cell type loses too many markers due to zero variance
  bulk.counts <- matrix(1:36,nrow=6,ncol=6)
  bulk.counts[1,] = rep(0, 6)
  bulk.eset <- Biobase::ExpressionSet(assayData = bulk.counts)
  markers <- data.frame(gene=as.character(1:6), cluster=c("a", "a", "a",
                                                          "b", "b", "b"))
  expect_error(BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, min_gene=3))
})

test_that("Produces output for simulated data", {
  bulk.eset <- Biobase::ExpressionSet(assayData = matrix(1:16, nrow=4, ncol=4))
  markers <- data.frame(gene=as.character(1:4), cluster=rep('a', 4), avg_logFC=rep(.3, 4))
  # weighted
  res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, min_gene=4, weighted=T)
  expect_true("bulk.props" %in% unlist(attributes(res)))
  # unweighted
  res <- BisqueRNA::MarkerBasedDecomposition(bulk.eset, markers, min_gene=4, weighted=F)
  expect_true("bulk.props" %in% unlist(attributes(res)))
})