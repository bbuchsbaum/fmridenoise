context("project_cone_heuristic - HRF coefficient adjustments")

test_that("project_cone_heuristic enforces cone constraints and handles edge cases", {
  # Case 1: typical HRF coefficients with non-zero first element and negative second
  coeffs1 <- c(1, -2, 3)
  adjusted1 <- project_cone_heuristic(coeffs1)
  expect_equal(adjusted1[1], 0)
  expect_equal(adjusted1[2], 0)
  expect_equal(adjusted1[3], 2)

  # Case 2: coefficients become flat after adjustment and get small bump
  coeffs2 <- c(0, -0.1, 0)
  adjusted2 <- project_cone_heuristic(coeffs2)
  expect_true(adjusted2[1] > 0)
  expect_equal(adjusted2[2], 0)
  expect_equal(adjusted2[3], 0)

  # Case 3: zero-length input returns numeric(0)
  adjusted3 <- project_cone_heuristic(numeric(0))
  expect_length(adjusted3, 0)
})
