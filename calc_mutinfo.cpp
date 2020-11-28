#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


double xlogx (double x) {
  if (x > 0) {
    return (x * log(x));
  } else {
    return 0.0;
  }
}



// returns 0 if x <= 0
double zero_log (double x) {
  if (x > 0) {
    return log(x);
  } else {
    return 0.0;
  }
}


double shannon(double prob) {
  return xlogx(prob) + xlogx(1-prob);
}




// [[Rcpp::export]]
NumericVector calc_minfo_normalized_C(
    NumericMatrix prob_mat, 
    arma::vec wt_vec, 
    bool normalize = true) {
  unsigned int n = prob_mat.nrow(), K = prob_mat.ncol();
  double p_joint_1, p_joint_0, p_row_marginal_1, p_row_marginal_0, tmp;
  
  arma::vec log_wt_vec = log(wt_vec);
  
  double sqrt_H_wt = sqrt(- arma::sum(log_wt_vec % wt_vec));
  double sqrt_H_row_marginal;
  
  // Rcpp::Rcout << "sqrt_H_wt: " << sqrt_H_wt << std::endl;
  
  NumericVector minfo_vec(n); 
  double tmp2;
  
  for(unsigned int i = 0; i < n; i++) {
    p_row_marginal_0 = 0;
    p_row_marginal_1 = 0;
    tmp = 0;
    for(unsigned int j = 0; j < K; j++) {
      p_joint_1 = prob_mat(i, j) * wt_vec[j];
      p_joint_0 = (1-prob_mat(i, j)) * wt_vec[j];
      
      p_row_marginal_1 += p_joint_1;
      p_row_marginal_0 += p_joint_0;
      
      tmp += xlogx(p_joint_1) - p_joint_1 * log_wt_vec[j] +
        xlogx(p_joint_0) - p_joint_0 * log_wt_vec[j];
    }
    
    // Rcpp::Rcout << "p_row_marginal_1: " << p_row_marginal_1 << std::endl;
    // Rcpp::Rcout << "xlogx(p_row_marginal_1): " << xlogx(p_row_marginal_1) << std::endl;
    // 
    // Rcpp::Rcout << "p_row_marginal_0: " << p_row_marginal_0 << std::endl;
    // Rcpp::Rcout << "xlogx(p_row_marginal_0): " << xlogx(p_row_marginal_0) << std::endl;
    
    sqrt_H_row_marginal = sqrt(
      - (xlogx(p_row_marginal_0) + xlogx(p_row_marginal_1))
    );
    
    // Rcpp::Rcout << "sqrt_H_row_marginal: " << sqrt_H_row_marginal << std::endl;
    
    tmp2 = tmp - xlogx(p_row_marginal_1) - xlogx(p_row_marginal_0);
    
    // Rcpp::Rcout << "MI: " << tmp2 << std::endl;
    
    minfo_vec[i] = tmp2;
    
    if (normalize) {
      minfo_vec[i] /= (sqrt_H_row_marginal * sqrt_H_wt); 
    }
    
    
    
    // Rcpp::Rcout << "NMI: " << minfo_vec[i] << std::endl;
  }
  
  
  return minfo_vec;
}



// [[Rcpp::export]]
arma::mat calc_minfo_binarycat_normalized_C(
    arma::mat prob_mat, 
    arma::vec wt_vec,
    bool normalized = true) {
  
  unsigned int n = prob_mat.n_rows, K = prob_mat.n_cols;
  
  arma::vec sqrt_H_wt(K);
  for (unsigned int j = 0; j < K; j++) {
    sqrt_H_wt[j] = sqrt(- shannon(wt_vec[j]));
  }
  
  arma::mat out(n, K);
  double p_row_marginal_1, p_row_marginal_0, 
  p1c, p0c, p1cbar, p0cbar, pc, pcbar,
  sqrt_H_row_marginal;
  
  
  for (unsigned int i = 0; i < n; i++) {
    p_row_marginal_0 = 0;
    p_row_marginal_1 = 0;
    
    for (unsigned int j = 0; j < K; j++) {
      p1c = prob_mat(i, j) * wt_vec[j];
      p0c = (1-prob_mat(i, j)) * wt_vec[j];
      
      p_row_marginal_1 += p1c;
      p_row_marginal_0 += p0c;
    }
    
    sqrt_H_row_marginal = sqrt(
      - (xlogx(p_row_marginal_0) + xlogx(p_row_marginal_1)
      )
    );
    
    
    for (unsigned int j = 0; j < K; j++) {
      p1c = prob_mat(i, j) * wt_vec[j];
      p1cbar = p_row_marginal_1 - p1c;
      
      p0c = (1-prob_mat(i, j)) * wt_vec[j];
      p0cbar = p_row_marginal_0 - p0c;
      
      pc = wt_vec[j];
      pcbar = 1 - wt_vec[j];
      
      out(i, j) = 
        p1c * (zero_log(p1c) - zero_log(p_row_marginal_1) - zero_log(pc)) +
        p0c * (zero_log(p0c) - zero_log(p_row_marginal_0) - zero_log(pc)) +
        p1cbar * (zero_log(p1cbar) - zero_log(p_row_marginal_1) - zero_log(pcbar)) +
        p0cbar * (zero_log(p0cbar) - zero_log(p_row_marginal_0) - zero_log(pcbar));
      
      if (normalized) {
        out(i, j) /= (sqrt_H_row_marginal * sqrt_H_wt[j]);
      }
    }
    
  }
  
  return (out);
  
}