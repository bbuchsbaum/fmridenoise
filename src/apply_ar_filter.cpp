#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix apply_ar_filter_matrix_cpp(const NumericMatrix& M,
                                         const NumericVector& ar_coeffs) {
    int n = M.nrow();
    int p = ar_coeffs.size();
    int m = M.ncol();
    NumericMatrix res(n, m);

    if (p < 1) {
        std::copy(M.begin(), M.end(), res.begin());
        return res;
    }

    bool use_filter = false;
    for (int k = 0; k < p; ++k) {
        double phi = ar_coeffs[k];
        if (!NumericVector::is_na(phi) && phi != 0.0) {
            use_filter = true;
            break;
        }
    }

    if (!use_filter) {
        std::copy(M.begin(), M.end(), res.begin());
        return res;
    }

    #pragma omp parallel for
    for (int col = 0; col < m; ++col) {
        for (int t = 0; t < p && t < n; ++t) {
            res(t, col) = NA_REAL;
        }
        for (int t = p; t < n; ++t) {
            double val = M(t, col);
            for (int k = 0; k < p; ++k) {
                val -= ar_coeffs[k] * M(t - k - 1, col);
            }
            res(t, col) = val;
        }
    }
    return res;
}

// [[Rcpp::export]]
NumericMatrix apply_ar_filter_voxelwise_cpp(const NumericMatrix& Y,
                                            const NumericMatrix& coeffs) {
    int n = Y.nrow();
    int m = Y.ncol();
    int p = coeffs.ncol();
    if (coeffs.nrow() != m)
        stop("Coefficient matrix rows must equal number of columns in Y");

    NumericMatrix res(n, m);
    #pragma omp parallel for
    for (int col = 0; col < m; ++col) {
        bool use_filter = false;
        for (int k = 0; k < p; ++k) {
            double phi = coeffs(col, k);
            if (!NumericVector::is_na(phi) && phi != 0.0) {
                use_filter = true;
                break;
            }
        }
        if (!use_filter) {
            for (int t = 0; t < n; ++t) {
                res(t, col) = Y(t, col);
            }
            continue;
        }
        for (int t = 0; t < p && t < n; ++t) {
            res(t, col) = NA_REAL;
        }
        for (int t = p; t < n; ++t) {
            double val = Y(t, col);
            for (int k = 0; k < p; ++k) {
                val -= coeffs(col, k) * Y(t - k - 1, col);
            }
            res(t, col) = val;
        }
    }
    return res;
}

