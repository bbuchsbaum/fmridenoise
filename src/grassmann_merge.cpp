#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List grassmann_merge_iterative_cpp(List V_list, int k_target_global) {
    int n_list = V_list.size();
    if (n_list == 0 || k_target_global <= 0) {
        return List::create(Named("V_global") = NumericMatrix(0,0),
                            Named("svd_times") = NumericVector(0));
    }

    mat V_global = as<mat>(V_list[0]);
    if ((int)V_global.n_cols > k_target_global) {
        V_global = V_global.cols(0, k_target_global - 1);
    }

    NumericVector svd_times(std::max(0, n_list - 1));

    if (n_list > 1) {
        for (int i = 1; i < n_list; ++i) {
            mat Vr = as<mat>(V_list[i]);
            if (Vr.n_cols == 0 || Vr.n_rows != V_global.n_rows) continue;

            mat temp = join_rows(V_global, Vr);
            mat Q, R;
            qr_econ(Q, R, temp);
            mat P_union = Q;

            mat M_proj = P_union.t() * (V_global * V_global.t() + Vr * Vr.t()) * P_union;
            int k_for_svd = std::min(k_target_global, (int)std::min(M_proj.n_rows, M_proj.n_cols));
            if (k_for_svd <= 0) continue;

            mat U;
            vec s;
            double t1 = clock();
            svd_econ(U, s, M_proj);
            double elapsed = (clock() - t1) / (double)CLOCKS_PER_SEC;
            svd_times[i - 1] = elapsed;

            if ((int)U.n_cols > k_for_svd) {
                U = U.cols(0, k_for_svd - 1);
            }
            V_global = P_union * U;
        }
    }

    if ((int)V_global.n_cols > k_target_global) {
        V_global = V_global.cols(0, k_target_global - 1);
    }

    return List::create(Named("V_global") = V_global,
                        Named("svd_times") = svd_times);
}
