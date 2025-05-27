#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// Declaration of C++ function
List grassmann_merge_iterative_cpp(List V_list, int k_target_global);

// Wrapper for R
extern "C" SEXP _ndx_grassmann_merge_iterative_cpp(SEXP V_listSEXP, SEXP k_target_globalSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< List >::type V_list(V_listSEXP);
    Rcpp::traits::input_parameter< int >::type k_target_global(k_target_globalSEXP);
    return Rcpp::wrap(grassmann_merge_iterative_cpp(V_list, k_target_global));
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ndx_grassmann_merge_iterative_cpp", (DL_FUNC) &_ndx_grassmann_merge_iterative_cpp, 2},
    {NULL, NULL, 0}
};

extern "C" void R_init_ndx(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
