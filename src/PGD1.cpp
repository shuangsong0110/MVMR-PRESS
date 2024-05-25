// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
NumericVector matrixVectorMult(NumericMatrix mat, NumericVector vec) {
  Map<MatrixXd> X(Rcpp::as<Map<MatrixXd> >(mat));
  Map<VectorXd> Y(Rcpp::as<Map<VectorXd> >(vec));
  VectorXd result = X * Y;
  return Rcpp::wrap(result);
}


NumericVector prox(NumericMatrix A, NumericVector theta, NumericVector b, double tt, double lambda) {
  NumericVector temp = matrixVectorMult(A,theta) ;
  //Rcout <<  temp   << "\n";
  NumericVector s = theta - tt * (temp - b);
  //Rcout <<  s   << "\n";
  int n = s.size();
  for (int i = 0; i < n; i++) {
    if (s[i] > lambda * tt) {
      s[i] = s[i] - lambda * tt;
    } else if (s[i] < -lambda * tt) {
      s[i] = s[i] + lambda * tt;
    } else {
      s[i] = 0;
    }
  }
  return s;
}

// [[Rcpp::export]]
double g(NumericVector theta, NumericMatrix A, NumericVector b) {
  double result = 0;
  int n = theta.size();
  for (int i = 0; i < n; i++) {
    double tmp = 0;
    for (int j = 0; j < n; j++) {
      tmp += A(i, j) * theta[j];
    }
    result += 0.5 * theta[i] * tmp;
  }
  result -= sum(b * theta);
  return result;
}

// [[Rcpp::export]]
NumericVector Gr(NumericVector theta, NumericMatrix A, NumericVector b, double tt, double lambda) {
  return (theta - prox(A, theta, b, tt, lambda)) / tt;
}

// [[Rcpp::export]]
NumericVector PGD(NumericMatrix A, NumericVector b, NumericVector theta, double lambda, int iter = 1000) {
  double beta = 0.2;
  int n = theta.size();
  for (int i = 0; i < iter; i++) {
    double tt = 1;
    NumericVector deg = matrixVectorMult(A, theta) - b;
    //NumericVector temp = A * theta;
    //Rcout <<  deg   << "\n";
    while (g(theta - tt * Gr(theta, A, b, tt, lambda), A, b) > (g(theta, A, b) - tt * sum(deg * Gr(theta, A, b, tt, lambda)) + 0.5 * tt * sum(pow(Gr(theta, A, b, tt, lambda), 2)))) {
      tt = beta * tt;
    }
    theta = prox(A, theta, b, tt, lambda);
  }
  return theta;
}
