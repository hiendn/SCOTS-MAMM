#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

double PSEUDO_FUN(SEXP CC, SEXP ETA, SEXP PARA, int NN, int GG) {
  
  // Load Primary Variables 
  
  Rcpp::NumericMatrix CC_C(CC);
  arma::mat CC_A(CC_C.begin(),NN,GG,false);
  
  Rcpp::NumericMatrix ETA_C(ETA);
  arma::mat ETA_A(ETA_C.begin(),NN,GG,false);
  
  Rcpp::NumericMatrix PARA_C(PARA);
  arma::mat PARA_A(PARA_C.begin(),2,GG,false);  
  
  // Initiate Likelihood Value
  double LL = 0;
  
  // Initiate Loop
  
  for (int nn = 0; nn < NN; nn++) {
    
    // Initiate Dummy Variable
    double INNER_1 = 0;
    double INNER_2 = 0;

    for (int gg = 0; gg < GG; gg++) {
      
      // Make New Matrix
      arma::mat DAT(1,2);
      DAT(0,0) = 1;
      DAT(0,1) = ETA_A(nn,gg);
    
      
      // Combine Dummies
      INNER_1 = INNER_1 + exp(as_scalar(DAT*PARA_A.col(gg)));
      INNER_2 = INNER_2 + CC_A(nn,gg)*as_scalar(DAT*PARA_A.col(gg));
    }
    
    // Add to Log Likelihood
    LL = LL + INNER_2 - log(INNER_1);
  }
  
  return LL;
}