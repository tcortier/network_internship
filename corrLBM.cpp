
#include "Rcpp.h"
#include "utils.h"

using namespace Rcpp;

long fact(int n){

     return (n==0) || (n==1) ? 1 : n* fact(n-1);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix corr_update_tau1(
    Rcpp::NumericMatrix Y,
    Rcpp::NumericMatrix sampl,
    Rcpp::NumericVector alpha1,
    Rcpp::NumericMatrix pi,
    Rcpp::NumericMatrix tau2) {

  int N1 = Y.nrow();
  int N2 = Y.ncol();
  int Q1 = pi.nrow();
  int Q2 = pi.ncol();
  double acc;
  NumericMatrix new_tau1(N1,Q1);

  for(int i=0; i < N1; i++) {
    for(int q=0; q < Q1; q++){

      acc = 0;
      for (int j=0; j < N2; j++) {
        for (int l=0; l < Q2; l++) {
          acc = acc + tau2(j,l) * Y(i,j) *  std::log(pi(q,l) * sampl(i,j)) + tau2(j,l) * (1-Y(i,j))*log(1-sampl(i,j)*pi(q,l)) ;
        }
      }
      new_tau1(i,q) = alpha1[q] * std::exp(acc);
    }
    new_tau1(i,_) = new_tau1(i,_)/sum(new_tau1(i,_));
  }

  return Rcpp::wrap(new_tau1);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix corr_update_tau1_poisson(
    Rcpp::NumericMatrix Y,
    Rcpp::NumericMatrix sampl,
    Rcpp::NumericVector alpha1,
    Rcpp::NumericMatrix lambda,
    Rcpp::NumericMatrix tau2) {

  int N1 = Y.nrow();
  int N2 = Y.ncol();
  int Q1 = lambda.nrow();
  int Q2 = lambda.ncol();
  double acc;
  NumericMatrix new_tau1(N1,Q1);

  for(int i=0; i < N1; i++) {
    for(int q=0; q < Q1; q++){

      acc = 0;
      for (int j=0; j < N2; j++) {
        for (int l=0; l < Q2; l++) {
          acc = acc + tau2(j,l) * Y(i,j) *  std::log(lambda(q,l) * sampl(i,j)) - tau2(j,l)* sampl(i,j)*lambda(q,l)- tau2(j,l)*std::log(fact(Y(i,j)));
        }
      }
      new_tau1(i,q) = alpha1[q] * std::exp(acc);
    }
    new_tau1(i,_) = new_tau1(i,_)/sum(new_tau1(i,_));
  }

  return Rcpp::wrap(new_tau1);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix corr_update_tau2(
    Rcpp::NumericMatrix Y,
    Rcpp::NumericMatrix sampl,
    Rcpp::NumericVector alpha2,
    Rcpp::NumericMatrix pi,
    Rcpp::NumericMatrix tau1) {

  int N1 = Y.nrow();
  int N2 = Y.ncol();
  int Q1 = pi.nrow();
  int Q2 = pi.ncol();
  double acc;
  NumericMatrix new_tau2(N2,Q2);

  for(int j=0; j < N2; j++) {
    for(int l=0; l < Q2; l++){

      acc = 0;
      for (int i=0; i < N1; i++) {
        for (int q=0; q < Q1; q++) {
          acc = acc + tau1(i,q) * Y(i,j) *  std::log(pi(q,l) * sampl(i,j)) + tau1(i,q) * (1-Y(i,j))*log(1-sampl(i,j)*pi(q,l)) ;
        }
      }
      new_tau2(j,l) = alpha2[l] * std::exp(acc);
    }
    new_tau2(j,_) = new_tau2(j,_)/sum(new_tau2(j,_));
  }

  return Rcpp::wrap(new_tau2);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix corr_update_tau2_poisson(
    Rcpp::NumericMatrix Y,
    Rcpp::NumericMatrix sampl,
    Rcpp::NumericVector alpha2,
    Rcpp::NumericMatrix lambda,
    Rcpp::NumericMatrix tau1) {

  int N1 = Y.nrow();
  int N2 = Y.ncol();
  int Q1 = lambda.nrow();
  int Q2 = lambda.ncol();
  double acc;
  NumericMatrix new_tau2(N2,Q2);

  for(int j=0; j < N2; j++) {
    for(int l=0; l < Q2; l++){

      acc = 0;
      for (int i=0; i < N1; i++) {
        for (int q=0; q < Q1; q++) {
          acc = acc + tau1(i,q) * Y(i,j) *  std::log(lambda(q,l) * sampl(i,j)) - tau1(i,q) * sampl(i,j)*lambda(q,l)- tau1(i,q)*std::log(fact(Y(i,j)));
        }
      }
      new_tau2(j,l) = alpha2[l] * std::exp(acc);
    }
    new_tau2(j,_) = new_tau2(j,_)/sum(new_tau2(j,_));
  }

  return Rcpp::wrap(new_tau2);
}


// [[Rcpp::export]]
double corr_compute_logL(
    Rcpp::NumericMatrix members1,
    Rcpp::NumericMatrix members2,
    Rcpp::NumericMatrix pi,
    Rcpp::NumericMatrix Y,
    Rcpp::NumericMatrix sampl
    ){
    
  int Q1 = pi.nrow();
  int N1 = Y.nrow();
  int N2 = Y.ncol();
  int Q2 = pi.ncol();
  double LL;
  LL=0;
  for(int i=0; i < N1; i++) {
    for(int q=0; q < Q1; q++){
      for (int j=0; j < N2; j++) {
        for (int l=0; l < Q2; l++) {
          LL = LL + members1(i,q) * members2(j,l) * Y(i,j) * std::log(sampl(i,j) * pi(q,l)) + members1(i,q) * members2(j,l) * (1-Y(i,j)) * std::log(1-sampl(i,j) * pi(q,l));
        }
      }
    }
  }
  return LL;
}


// [[Rcpp::export]]
double function_essai(
    Rcpp::NumericMatrix Y
    ) {

  int N1 = Y.nrow();
  int N2 = Y.ncol();
  double acc;
  acc = 0;
  for (int i=0; i < N1; i++) {
    for (int j=0; j < N2; j++) {
      acc = acc + Y(i,j) ;
    }
  }    
  return acc;
}




// [[Rcpp::export]]
double corr_compute_logL_poisson(
    Rcpp::NumericMatrix members1,
    Rcpp::NumericMatrix members2,
    Rcpp::NumericMatrix lambda,
    Rcpp::NumericMatrix Y,
    Rcpp::NumericMatrix sampl
    ){
    

  int N1 = Y.nrow();
  int N2 = Y.ncol();
  int Q1 = lambda.nrow();
  int Q2 = lambda.ncol();
  double LL;
  LL=0;
  for(int i=0; i < N1; i++) {
    for(int q=0; q < Q1; q++){
      for (int j=0; j < N2; j++) {
        for (int l=0; l < Q2; l++) {
          LL = LL + members1(i,q) * members2(j,l) * Y(i,j) * std::log(sampl(i,j) * lambda(q,l)) - members1(i,q) * members2(j,l)*sampl(i,j)*lambda(q,l)-members1(i,q)*members2(j,l)*std::log(fact(Y(i,j)));
        }
      }
    }
  }
  return LL;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix essai(
    Rcpp::NumericMatrix members1,
    Rcpp::NumericMatrix members2,
    Rcpp::NumericMatrix lambda,
    Rcpp::NumericMatrix Y,
    Rcpp::NumericMatrix sampl
    ){
    

  int N1 = Y.nrow();
  int N2 = Y.ncol();
  NumericMatrix facY(N1,N2);
  

  for(int i=0; i < N1; i++) {
    for (int j=0; j < N2; j++) {
      facY(i,j) = std::log(fact(Y(i,j)));
    }
  }

  return facY;
}




