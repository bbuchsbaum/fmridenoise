context(".grassmann_merge_iterative - Orthonormal merging")

# Helper function to check orthonormality of a matrix
.is_orthonormal <- function(M, tol = 1e-6) {
  if (is.null(M) || ncol(M) == 0) return(FALSE)
  I_est <- crossprod(M)
  all(dim(I_est) == c(ncol(M), ncol(M))) &&
    max(abs(I_est - diag(ncol(M)))) < tol
}

# Construct small orthonormal matrices
V1 <- diag(5)[, 1:2]
V2 <- diag(5)[, 3:4]
V3 <- qr.Q(qr(matrix(rnorm(25), 5))) # 5x5 orthonormal, take last column
V3 <- V3[,5, drop=FALSE]

merged <- .grassmann_merge_iterative(list(V1, V2, V3), k_target_global = 3L)


test_that("Merged matrix has correct dimensions", {
  expect_true(!is.null(merged))
  expect_equal(ncol(merged), 3)
  expect_equal(nrow(merged), 5)
})

test_that("Merged matrix remains orthonormal", {
  expect_true(.is_orthonormal(merged, tol = 1e-6))
})

