context("test-normalization")


test_that("Test that normalization's output type is matrix", {
  counts_matrix <- as.matrix(rbind(c(0,1,2), c(3,4,5), c(6,7,8)))
  counts_dataframe <- as.data.frame(rbind(c(0,1,2), c(3,4,5), c(6,7,8)))
  expect_is(normalization(counts_matrix), "matrix")
  expect_is(normalization(counts_dataframe), "matrix")
})


test_that("Test that normalization's output dimension is correct", {
  counts <- rbind(c(0,1,2,3), c(4,5,6,7), c(8,9,10,11))
  expect_equal(dim(normalization(counts)), dim(counts))
  counts_t <- t(counts)
  expect_equal(dim(normalization(counts_t)), dim(counts_t))
})


test_that("Test that normalization's output is correct", {
  counts <- rbind(c(0,1,2), c(3,4,5), c(6,7,8))
  expect_equal(normalization(counts), rbind(c(0,1,2), c(4,4,4), c(8,7,7)))
})




