context("test-estDropoutNum")


test_that("Test that estDropoutNum's output length is correct", {
  expect_equal(length(estDropoutNum(c(0,0,0,0,1,1,2,3))), 1L)
  expect_equal(length(estDropoutNum(c(0,0,0,0,1,2,3,3))), 1L)
  expect_equal(length(estDropoutNum(c(0,0,0,0,1,2,3,5))), 1L)
})


test_that("Test that estDropoutNum's dropoutNum output is correct", {
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3)), 2)
  expect_equal(estDropoutNum(c(0,0,0,0,1,2,3,3)), 1)
  expect_equal(estDropoutNum(c(0,0,0,0,1,2,3,5)), 0)
})


test_that("Test that estDropoutNum's dropoutNum output is correct with different depth", {
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3), depth = 2), 1)
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3), depth = 20), 2)
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3), depth = 50), 2)
})


test_that("Test that estDropoutNum's geneNum output is correct", {
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3), return = "geneNumPredict"), 6)
  expect_equal(estDropoutNum(c(0,0,0,0,1,2,3,3), return = "geneNumPredict"), 5)
  expect_equal(estDropoutNum(c(0,0,0,0,1,2,3,5), return = "geneNumPredict"), 4)
})


test_that("Test that estDropoutNum's geneNum output is correct with different depth", {
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3), depth = 2, return = "geneNumPredict"), 5)
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3), depth = 20, return = "geneNumPredict"), 6)
  expect_equal(estDropoutNum(c(0,0,0,0,1,1,2,3), depth = 50, return = "geneNumPredict"), 6)
})




