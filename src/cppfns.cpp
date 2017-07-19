#include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
NumericVector collectC(NumericVector xs, NumericVector ys, int n) {
  int ln = ys.size();
  NumericVector out(n);
  
  for(int i = 0; i < ln; ++i) {
    out[xs[i]] += ys[i];
  }
  return out;
}


// [[Rcpp::export]]
NumericVector takefirst(NumericVector xs, int pos) {
  NumericVector out(pos);
  
  for(int i = 0; i < pos; ++i) {
    out[i] = xs[i];
  }
  return out;
}

// [[Rcpp::export]]
NumericVector takelast(NumericVector xs, int pos, int n) {
  
  NumericVector out(n-pos);
  
  for(int i = pos; i < n; ++i) {
    out[i-pos] = xs[i];
  }
  return out;
}



