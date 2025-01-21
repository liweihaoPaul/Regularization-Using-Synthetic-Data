// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
//using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export]]
double estimate_VELOOCV_cpp(const Eigen::MatrixXd & X, 
                               const Eigen::VectorXd & y1, 
                               const Eigen::MatrixXd & Xstar, 
                               const Eigen::VectorXd & ystar, 
                               double tau_0, 
                               const Eigen::VectorXd & betahat_tau0) {
  
  auto sigmoid = [](double t) {
    return 1.0 / (1.0 + exp(-t));
  };
  
  auto deri_sigmoid = [&sigmoid](double t) {
    double s = sigmoid(t);
    return s * (1 - s);
  };
  
  int n = X.rows();
  int M = Xstar.rows();
  int p = X.cols();
  
  Eigen::MatrixXd H(p, p);
  Eigen::MatrixXd Hstar(p, p);
  H.setZero();
  Hstar.setZero();
  
  for(int j = 0; j < n; j++) {
    double tmp = X.row(j).dot(betahat_tau0);
    H.noalias() -= deri_sigmoid(tmp) * X.row(j).transpose() * X.row(j);
  }
  
  for(int js = 0; js < M; js++) {
    double tmp = Xstar.row(js).dot(betahat_tau0);
    Hstar.noalias() -= (tau_0 * n / M) * deri_sigmoid(tmp) * Xstar.row(js).transpose() * Xstar.row(js);
  }

  Eigen::MatrixXd AA = (Hstar + H).inverse();
  
  Eigen::VectorXd betahat_minus_i_T_Xi_seq(n);
  for(int i = 0; i < n; i++) {
    double tmp = X.row(i).dot(betahat_tau0);
    double deritmp=deri_sigmoid(tmp);
    Eigen::VectorXd AAXrowiT= AA *  X.row(i).transpose() ;
    // Adjust the multiplication order and the division to make sure dimensions align
    Eigen::MatrixXd tmp_matrix = AAXrowiT * AAXrowiT.transpose()* deritmp;
    Eigen::MatrixXd A_minus_i = AA - tmp_matrix/ (1.0 + deritmp*  X.row(i) * AAXrowiT);

    // Adjusted the dot product to get qi.
    double qi = tmp + (y1[i] - sigmoid(tmp)) * X.row(i) * A_minus_i * X.row(i).transpose();
    betahat_minus_i_T_Xi_seq[i] = qi;
  }
  


  double  VE = 0.0;
  for(int i = 0; i < n; i++) {
    VE += y1[i]*betahat_minus_i_T_Xi_seq[i]-log(1+exp(betahat_minus_i_T_Xi_seq[i])  );
  }
  

  return -VE;
}
