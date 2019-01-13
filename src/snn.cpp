// snn.cpp
//
// Based on the implementation in Seurat.

#include <RcppEigen.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::Triplet<double> triplet_t ;

// compute k-nearest neighbors
std::vector<triplet_t> compute_knn_impl(Eigen::MatrixXd x) {
  std::vector<triplet_t> KNN;
  KNN.reserve(x.rows() * x.cols());

  for (int j=0; j <  x.cols(); ++j) {
    for (int i=0; i < x.rows(); ++i) {
      KNN.push_back(triplet_t(i, x(i, j) - 1, 1));
    }
  }

  return KNN ;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> compute_snn_impl(Eigen::MatrixXd x, double prune) {

  int k = x.cols();

  std::vector<triplet_t> KNN = compute_knn_impl(x) ;

  Eigen::SparseMatrix<double> SNN(x.rows(), x.rows()) ;
  SNN.setFromTriplets(KNN.begin(), KNN.end()) ;

  SNN = SNN * (SNN.transpose());

  for (int i=0; i < SNN.outerSize(); ++i) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(SNN, i); it; ++it) {
      it.valueRef() = it.value() / (k + (k - it.value()));
      if (it.value() < prune) {
        it.valueRef() = 0;
      }
    }
  }

  SNN.prune(0.0) ;
  return SNN ;
}
