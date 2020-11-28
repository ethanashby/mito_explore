Rcpp::sourceCpp("calc_mutinfo.cpp")

calc_minfo <- function(prob_mat, cancer_prob = 1/ncol(prob_mat), 
                       normalize = TRUE, binary_minfo = FALSE) {
  
  if (any(is.null(names(cancer_prob)))) {
    warning("'cancer_prob' has no names. Matching them in the column order of prob_mat")
  }
  
  if (any(is.null(colnames(prob_mat)))) {
    warning("'cancer_prob' has no column names. Matching columns with the order of cancer_prob")
  }
  
  cancer_prob <- cancer_prob/sum(cancer_prob)
  
  
  if (all(!is.null(colnames(prob_mat)), !is.null(names(cancer_prob)))) {
    cancer_names <- colnames(prob_mat)
    cancer_prob <- cancer_prob[cancer_names]
  }
  
  
  if (!binary_minfo) {
    minfo <- calc_minfo_normalized_C(
      data.matrix(prob_mat), 
      cancer_prob, 
      normalize
    )
    names(minfo) <- rownames(prob_mat)
  } else {
    minfo <- calc_minfo_binarycat_normalized_C(
      data.matrix(prob_mat), 
      cancer_prob, 
      normalize
    )
    rownames(minfo) <- rownames(prob_mat)
    colnames(minfo) <- colnames(prob_mat)
  }
  
  minfo
}