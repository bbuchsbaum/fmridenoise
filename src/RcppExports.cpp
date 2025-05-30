#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// Declaration of C++ function
List grassmann_merge_iterative_cpp(List V_list, int k_target_global);
Rcpp::NumericMatrix apply_ar_filter_matrix_cpp(const Rcpp::NumericMatrix& M, const Rcpp::NumericVector& ar_coeffs);
Rcpp::NumericMatrix apply_ar_filter_voxelwise_cpp(const Rcpp::NumericMatrix& Y, const Rcpp::NumericMatrix& coeffs);

// Wrapper for R
extern "C" SEXP _ndx_grassmann_merge_iterative_cpp(SEXP V_listSEXP, SEXP k_target_globalSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< List >::type V_list(V_listSEXP);
    Rcpp::traits::input_parameter< int >::type k_target_global(k_target_globalSEXP);
    return Rcpp::wrap(grassmann_merge_iterative_cpp(V_list, k_target_global));
END_RCPP
}

extern "C" SEXP _ndx_apply_ar_filter_matrix_cpp(SEXP MSEXP, SEXP ar_coeffsSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type ar_coeffs(ar_coeffsSEXP);
    return Rcpp::wrap(apply_ar_filter_matrix_cpp(M, ar_coeffs));
END_RCPP
}

extern "C" SEXP _ndx_apply_ar_filter_voxelwise_cpp(SEXP YSEXP, SEXP coeffsSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type coeffs(coeffsSEXP);
    return Rcpp::wrap(apply_ar_filter_voxelwise_cpp(Y, coeffs));
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ndx_grassmann_merge_iterative_cpp", (DL_FUNC) &_ndx_grassmann_merge_iterative_cpp, 2},
    {"_ndx_apply_ar_filter_matrix_cpp", (DL_FUNC) &_ndx_apply_ar_filter_matrix_cpp, 2},
    {"_ndx_apply_ar_filter_voxelwise_cpp", (DL_FUNC) &_ndx_apply_ar_filter_voxelwise_cpp, 2},
    {NULL, NULL, 0}
};

extern "C" void R_init_ndx(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
