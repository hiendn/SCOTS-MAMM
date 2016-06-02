#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

// Copyright Hien Duy Nguyen - University of Queensland 2016/06/02
arma::mat TAU_FUN(SEXP TAU, SEXP XX, SEXP YY, SEXP BETA, SEXP PII, SEXP VAR, int MM, int NN, int GG, int PP, int MAXPP) {
  
  // Load Primary Variables 
  Rcpp::List XX_R(XX);
  Rcpp::List YY_R(YY);
  int MM_LEFT = MM-PP;
  
  SEXP TAU_S(TAU);
  Rcpp::NumericMatrix TAU_C(TAU_S);
  arma::mat TAU_A(TAU_C.begin(),NN,GG,false);
  
  // Initiate Loop
  
  for (int nn = 0; nn < NN; nn++) {
    
    // Initiate Dummy Variable
    double INNER_1 = 0;
    
    
    // Load Variables
    SEXP XX_S(XX_R[nn]);
    Rcpp::NumericMatrix XX_C(XX_S);
    arma::mat XX_A(XX_C.begin(),MM_LEFT,PP+1,false);
    
    SEXP YY_S(YY_R[nn]);
    Rcpp::NumericVector YY_C(YY_S);
    arma::colvec YY_A(YY_C.begin(),MM_LEFT,false);

    Rcpp::NumericVector PI_C(PII);
    arma::colvec PI_A(PI_C.begin(),GG,false);
    
    Rcpp::NumericVector VAR_C(VAR);
    arma::colvec VAR_A(VAR_C.begin(),GG,false);

    // Load BETA
    Rcpp::NumericMatrix BETA_C(BETA);
    arma::mat BETA_A(BETA_C.begin(),PP+1,GG,false);

    // Constant Dummy
    double IN_PLUSD = 100000000;
    for (int gg = 0; gg < GG; gg++) {
      
      // Initiate New Dummy Variables
      double IN_PLUS = 0;
      double VAR_IN = VAR_C(gg);
      double SD_C = sqrt(VAR_IN);
      arma::colvec BETA_IN = BETA_A.col(gg);
      
      
      for (int mm = MAXPP-PP; mm < MM_LEFT; mm++) {
        
        IN_PLUS = IN_PLUS + R::dnorm(YY_A(mm), as_scalar(XX_A.row(mm)*BETA_IN), SD_C, TRUE);
      }
      
      if (IN_PLUS < IN_PLUSD) {
        IN_PLUSD = IN_PLUS;
      }
    }

    for (int gg = 0; gg < GG; gg++) {
      
      // Initiate New Dummy Variables
      double INNER_2 = PI_A(gg);
      double IN_PLUS = 0;
      double VAR_IN = VAR_C(gg);
      double SD_C = sqrt(VAR_IN);
      arma::colvec BETA_IN = BETA_A.col(gg);
      
      
      for (int mm = 0; mm < MM_LEFT; mm++) {
        
        IN_PLUS = IN_PLUS + R::dnorm(YY_A(mm), as_scalar(XX_A.row(mm)*BETA_IN), SD_C, TRUE);
      }
      
      // Combine Dummies
      INNER_1 = INNER_1 + INNER_2*(exp(IN_PLUS-IN_PLUSD)); 
    }
    
    for (int gg = 0; gg < GG; gg++) {
      

      
      // Initiate New Dummy Variables
      double INNER_2 = PI_A(gg);
      double IN_PLUS = 0;
      double VAR_IN = VAR_C(gg);
      double SD_C = sqrt(VAR_IN);
      arma::colvec BETA_IN = BETA_A.col(gg);
      
      
      for (int mm = 0; mm < MM_LEFT; mm++) {
        
        IN_PLUS = IN_PLUS + R::dnorm(YY_A(mm), as_scalar(XX_A.row(mm)*BETA_IN), SD_C, TRUE);
      }
      
      // Combine Dummies
      TAU_A(nn,gg) = INNER_2*(exp(IN_PLUS-IN_PLUSD))/INNER_1; 
    }
    
  }
  
  return TAU_A;
}
