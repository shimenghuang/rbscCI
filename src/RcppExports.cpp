// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// bbb_pvalue
double bbb_pvalue(int n1_tot, int n2_tot, double gamma, int n1_suc, int n2_suc, int p_step);
RcppExport SEXP _fastCI_bbb_pvalue(SEXP n1_totSEXP, SEXP n2_totSEXP, SEXP gammaSEXP, SEXP n1_sucSEXP, SEXP n2_sucSEXP, SEXP p_stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n1_tot(n1_totSEXP);
    Rcpp::traits::input_parameter< int >::type n2_tot(n2_totSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type n1_suc(n1_sucSEXP);
    Rcpp::traits::input_parameter< int >::type n2_suc(n2_sucSEXP);
    Rcpp::traits::input_parameter< int >::type p_step(p_stepSEXP);
    rcpp_result_gen = Rcpp::wrap(bbb_pvalue(n1_tot, n2_tot, gamma, n1_suc, n2_suc, p_step));
    return rcpp_result_gen;
END_RCPP
}
// bscCI
NumericVector bscCI(int n_tot, int n_suc, double conf);
RcppExport SEXP _fastCI_bscCI(SEXP n_totSEXP, SEXP n_sucSEXP, SEXP confSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_tot(n_totSEXP);
    Rcpp::traits::input_parameter< int >::type n_suc(n_sucSEXP);
    Rcpp::traits::input_parameter< double >::type conf(confSEXP);
    rcpp_result_gen = Rcpp::wrap(bscCI(n_tot, n_suc, conf));
    return rcpp_result_gen;
END_RCPP
}
// cpCI
NumericVector cpCI(int n_tot, int n_suc, double conf);
RcppExport SEXP _fastCI_cpCI(SEXP n_totSEXP, SEXP n_sucSEXP, SEXP confSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_tot(n_totSEXP);
    Rcpp::traits::input_parameter< int >::type n_suc(n_sucSEXP);
    Rcpp::traits::input_parameter< double >::type conf(confSEXP);
    rcpp_result_gen = Rcpp::wrap(cpCI(n_tot, n_suc, conf));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fastCI_bbb_pvalue", (DL_FUNC) &_fastCI_bbb_pvalue, 6},
    {"_fastCI_bscCI", (DL_FUNC) &_fastCI_bscCI, 3},
    {"_fastCI_cpCI", (DL_FUNC) &_fastCI_cpCI, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastCI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
