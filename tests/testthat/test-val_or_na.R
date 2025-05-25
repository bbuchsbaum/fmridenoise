context("ndx_val_or_na")

test_that("returns NA for NULL", {
  expect_true(is.na(ndx_val_or_na(NULL)))
})

test_that("returns value when not NULL", {
  expect_equal(ndx_val_or_na(5), 5)
  expect_equal(ndx_val_or_na("a"), "a")
})

