#include <Rcpp.h>
using namespace Rcpp;

#include "bbb.h"

// [[Rcpp::export(".bbb_pvalue")]]
double bbb_pvalue(unsigned int n1_tot, unsigned int n2_tot, double gamma, 
                  unsigned int n1_suc, unsigned int n2_suc, double p_step) {
  
  Barnard bbb(n1_tot, n2_tot, gamma);
  
  return bbb.p_value(n1_suc, n2_suc, p_step);
}

// [[Rcpp::export(".bbb_fast_pvalue")]]
double bbb_fast_pvalue(unsigned int n1_tot, unsigned int n2_tot, double gamma, 
                       unsigned int n1_suc, unsigned int n2_suc, int p_slot, double alpha) {
  
  BarnardFast bbbf(n1_tot, n2_tot, gamma, p_slot, alpha);
  
  return(bbbf.p_value(n1_suc, n2_suc));
}