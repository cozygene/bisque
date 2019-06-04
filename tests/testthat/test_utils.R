context("Expression matrix utlities")
library(Biobase)

test_that("CountsToCPM converts properly for valid data", {
  # 3 by 3
  example.counts <- matrix(1:9, nrow = 3, ncol = 3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  example.eset <- CountsToCPM(example.eset)
  expect_equal(unname(colSums(Biobase::exprs(example.eset))), rep(1000000, 3))
  # 3 by 1
  example.counts <- matrix(1:9, nrow = 3, ncol = 1)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  example.eset <- CountsToCPM(example.eset)
  expect_equal(unname(colSums(Biobase::exprs(example.eset))), rep(1000000, 1))
  # 1 by 3
  example.counts <- matrix(1:9, nrow = 1, ncol = 3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  example.eset <- CountsToCPM(example.eset)
  expect_equal(unname(colSums(Biobase::exprs(example.eset))), rep(1000000, 3))
  # 1 by 1
  example.counts <- matrix(1:9, nrow = 1, ncol = 1)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  example.eset <- CountsToCPM(example.eset)
  expect_equal(unname(colSums(Biobase::exprs(example.eset))), rep(1000000, 1))
})

test_that("CountsToCPM catches invalid samples with no expression", {
  # 3 by 3
  example.counts <- matrix(1:9, nrow = 3, ncol = 3)
  example.counts[,1] <- rep(0,3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_error(CountsToCPM(example.eset))
  # 3 by 1
  example.counts <- matrix(rep(0,3), nrow = 3, ncol = 1)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_error(CountsToCPM(example.eset))
  # 1 by 3
  example.counts <- matrix(1:9, nrow = 1, ncol = 3)
  example.counts[1,1] <- 0
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_error(CountsToCPM(example.eset))
  # 1 by 1
  example.counts <- matrix(0, nrow = 1, ncol = 1)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_error(CountsToCPM(example.eset))
})

test_that("FilterZeroVarianceGenes works", {
  # All have nonzero variance
  example.counts <- matrix(1:9, nrow = 3, ncol = 3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(dim(example.eset), dim(FilterZeroVarianceGenes(example.eset)))
  # Should remove 1 row
  example.counts[1,] <- rep(1,3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(2, unname(nrow(FilterZeroVarianceGenes(example.eset))))
  # If only one sample, should remove all genes
  example.counts <- matrix(1:9, nrow = 3, ncol = 1)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(0, unname(nrow(FilterZeroVarianceGenes(example.eset))))
  # If all have 0 variance, get rid of all genes
  example.counts <- matrix(1, nrow = 3, ncol = 3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(0, unname(nrow(FilterZeroVarianceGenes(example.eset))))
})

test_that("FilterUnexpressedGenes works", {
  # All are expressed
  example.counts <- matrix(1, nrow = 3, ncol = 3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(dim(example.eset), dim(FilterUnexpressedGenes(example.eset)))
  # Should remove 1 row
  example.counts[1,] <- rep(0,3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(2, unname(nrow(FilterUnexpressedGenes(example.eset))))
  # Only one sample, should still remove first row
  example.counts <- matrix(c(0,1,2), nrow = 3, ncol = 1)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(2, unname(nrow(FilterUnexpressedGenes(example.eset))))
  # All are unexpressed, should remove all rows
  example.counts <- matrix(0, nrow = 3, ncol = 3)
  example.eset <- Biobase::ExpressionSet(assayData = example.counts)
  expect_equal(0, unname(nrow(FilterUnexpressedGenes(example.eset))))
})

test_that("SeuratToExpressionSet works", {
  # v2 test
  setClass("testthatSeuratv2", representation(cell.names = "character", 
                                              ident = "character", 
                                              raw.data = "matrix"))
  sc.counts <- matrix(0,nrow=3,ncol=3)
  test.cell.names <- c("a-1", "b-2", "c-3")
  test.ident <- c("a", "b", "c")
  names(test.ident) <- test.cell.names
  colnames(sc.counts) <- test.cell.names
  test.seurat.obj <- new("testthatSeuratv2", 
                         cell.names=test.cell.names,
                         ident=test.ident,
                         raw.data=sc.counts
                         )
  expect_warning({test.eset <- SeuratToExpressionSet(test.seurat.obj, delimiter = '-',
                                     position=2, version = "v2")})
  expect_match(class(test.eset), "ExpressionSet")
})