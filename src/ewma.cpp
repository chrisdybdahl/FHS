#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ewma(NumericVector data, double init, double lambda) {
  int n = data.size();
  NumericVector ewma_values(n);
  ewma_values[0] = init;
  for (int i = 1; i < n; ++i) {
    ewma_values[i] = lambda * ewma_values[i - 1] + (1 - lambda) * data[i - 1] * data[i - 1];
  }
  return ewma_values;
}
