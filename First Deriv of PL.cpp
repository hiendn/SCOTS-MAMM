#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

// Copyright Hien Duy Nguyen - University of Queensland 2016/06/02
arma::mat PSEUDO_FUN(SEXP CC, SEXP NEIGH, SEXP PARA, int NN, int GG) {
  
  // Load Primary Variables 
  
  Rcpp::NumericMatrix CC_C(CC);
  arma::mat CC_A(CC_C.begin(),NN,GG,false);
  
  Rcpp::NumericMatrix NEIGH_C(NEIGH);
  arma::mat NEIGH_A(NEIGH_C.begin(),NN,GG,false);
  
  Rcpp::NumericMatrix PARA_C(PARA);
  arma::mat PARA_A(PARA_C.begin(),2,GG,false);  
  
  // Initiate Likelihood Value
  double LL = 0;
  
  // Initiate Loop
  
  for (int nn = 0; nn < NN; nn++) {
    
    // Initiate Dummy Variable
    arma::mat INNER_1.zeros(2,1);
    arma::mat INNER_1.zeros(2,1);

    for (int gg = 0; gg < GG; gg++) {
      
      // Make New Matrix
      arma::mat DAT(1,2);
      DAT(0,0) = 1;
      DAT(0,1) = NEIGH_A(nn,gg);
    
      
      // Combine Dummies
      INNER_1 = INNER_1 + exp(as_scalar(DAT*PARA_A.col(gg)));
      INNER_2 = INNER_2 + CC_A(nn,gg)*as_scalar(DAT*PARA_A.col(gg));
    }
    
    // Add to Log Likelihood
    LL = LL + log(INNER_1) + INNER_2;
  }
  
  return LL;
}
