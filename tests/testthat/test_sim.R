context("Basic simulation framwork")
library(Biobase)

test_that("Simulation provides output", {
  set.seed(42)
  sim.data <- BisqueRNA::SimulateData(3,3,2,c('a','b'), c(.5,.5))
  expect_match(class(sim.data$sc.eset), "ExpressionSet")
  expect_match(class(sim.data$bulk.eset), "ExpressionSet")
  expect_match(class(sim.data$props), "matrix")
  expect_match(class(sim.data$markers), "data.frame")
})