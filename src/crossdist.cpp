#include <Rcpp.h>

//' Computes  the Euclidean distance between rows of two matrices
//' @name crossdist
//' @param x First matrix
//' @param y Second matrix
//' @return Euclidean distance between rows of \code{x} and \code{y}
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix crossdist(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y) {
  int n1 = x.nrow(), n2 = y.nrow(), ncol = x.ncol(), i, j, k;

  if (ncol != y.ncol()) {
    throw std::runtime_error("Different column number");
  }

  Rcpp::NumericMatrix out(n1, n2);

  for (i = 0; i < n1; i++)
    for (j = 0; j < n2; j++) {
      double sum = 0;
      for (k = 0; k < ncol; k++)
        sum += pow(x(i, k) - y(j, k), 2);
      out(i, j) = sqrt(sum);
    }

  return out;
}
