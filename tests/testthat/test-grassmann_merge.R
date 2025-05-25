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


# Test equivalence of merging strategies

test_that("iterative and concat_svd merge yield equivalent subspace", {
  set.seed(123)
  V1 <- qr.Q(qr(matrix(rnorm(60), 10)))[,1:3]
  V2 <- qr.Q(qr(matrix(rnorm(60), 10)))[,1:3]
  k_target <- 4

  V_concat <- .ndx_basic_merge_voxel_components(list(V1, V2), k_target, strategy = "concat_svd")
  V_iter <- .ndx_basic_merge_voxel_components(list(V1, V2), k_target, strategy = "iterative")

  expect_equal(ncol(V_concat), k_target)
  expect_equal(ncol(V_iter), k_target)
  expect_equal(nrow(V_concat), 10)
  expect_equal(nrow(V_iter), 10)

  P_concat <- V_concat %*% t(V_concat)
  P_iter <- V_iter %*% t(V_iter)
  frob_diff <- norm(P_concat - P_iter, type = "F")
  expect_lt(frob_diff, 1e-6)
})

