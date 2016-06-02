#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

// Copyright Hien Duy Nguyen - University of Queensland 2016/06/02
arma::colvec VAR_FUN(SEXP TAU, SEXP XX, SEXP YY, SEXP BETA, SEXP PII, SEXP VAR, int MM, int NN, int GG, int PP, int MAXPP) {
  
  // Load Primary Variables 
  Rcpp::List XX_R(XX);
  Rcpp::List YY_R(YY);
  int MM_LEFT = MM-PP;

  SEXP TAU_S(TAU);
  Rcpp::NumericMatrix TAU_C(TAU_S);
  arma::mat TAU_A(TAU_C.begin(),NN,GG,false);
  
  Rcpp::NumericVector VAR_C(VAR);
  arma::colvec VAR_A(VAR_C.begin(),GG,false);

  // Load BETA
  Rcpp::NumericMatrix BETA_C(BETA);
  arma::mat BETA_A(BETA_C.begin(),PP+1,GG,false);

  // Initiate Loop
  for (int gg = 0; gg < GG; gg++) {
    
    // Initiate New Dummy Variables
    double INNER_1 = 0;
    double INNER_2 = 0;
    arma::colvec BETA_IN = BETA_A.col(gg);
    
    for (int nn = 0; nn < NN; nn++) {
      
      // Load Variables
      SEXP XX_S(XX_R[nn]);
      Rcpp::NumericMatrix XX_C(XX_S);
      arma::mat XX_A(XX_C.begin(),MM_LEFT,PP+1,false);
      
      SEXP YY_S(YY_R[nn]);
      Rcpp::NumericVector YY_C(YY_S);
      arma::colvec YY_A(YY_C.begin(),MM_LEFT,false);
             
      
      for (int mm = MAXPP-PP; mm < MM_LEFT; mm++) {
        
        INNER_1 = INNER_1 + TAU_A(nn,gg)*pow(YY_A(mm)-as_scalar(XX_A.row(mm)*BETA_IN),2);
        INNER_2 = INNER_2 + TAU_A(nn,gg);
      }
      
    }
    
    // VAR_A(gg) = 0.4*INNER_1/INNER_2 + 0.6*VAR_A(gg);
    VAR_A(gg) = INNER_1/INNER_2;
    
  }
  
  return VAR_A;
}
