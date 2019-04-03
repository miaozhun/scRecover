context("test-countsSampling")


test_that("Test that countsSampling samples correct fraction of reads", {
  expect_equal(sum(countsSampling(c(10,20,30), fraction = 0)), 0)
  expect_equal(sum(countsSampling(c(10,20,30), fraction = 0.5)), 30)
  expect_equal(sum(countsSampling(c(10,20,30), fraction = 1)), 60)
})


test_that("Test that countsSampling's output length is correct", {
  expect_equal(length(countsSampling(c(10), fraction = 0.5)), 1L)
  expect_equal(length(countsSampling(c(10,20), fraction = 0.5)), 2L)
  expect_equal(length(countsSampling(c(10,20,30), fraction = 0.5)), 3L)
})


test_that("Test that countsSampling's output is correct", {
  expect_equal(countsSampling(c(10,20,30), fraction = 0), c(0,0,0))
  expect_equal(countsSampling(c(10,20,30), fraction = 1), c(10,20,30))
})




