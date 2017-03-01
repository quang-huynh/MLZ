#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double MLe(NumericVector stpar, NumericVector Lbar, NumericVector ss,
           NumericVector eff, NumericVector LH, double Lc, double eff_init,
           int n_age, int n_season, int obs_season, double timing, int logpar) {
  double Linf = LH[0];
  double K = LH[1];
  double a0 = LH[2];
  double ac = a0 - log(1 - Lc/Linf)/K;

  int n_yr = Lbar.size();
  int y;
  int a;
  int k;
  double ndata = 0.;
  int astep = n_age*n_season;

  double q;
  double M;

  if(logpar == 0) {
    q = stpar[0];
    M = stpar[1];
  }
  if(logpar == 1) {
    q = exp(stpar[0]);
    M = exp(stpar[1]);
  }

  double Z_init = q * eff_init + M;
  NumericVector Z(n_yr);
  arma::cube N(n_yr,astep,n_season);
  NumericMatrix Nobs(n_yr,astep);
  NumericVector La(astep);
  NumericVector age(astep);

  NumericVector Lpred(n_yr);
  double sum_square = 0.;
  double sigma;
  double nLL;

  for(a=0;a<astep;a++) {
    double ageD = a;
    double seasD = n_season;
    age[a] = ac + ageD/seasD;
    La[a] = Linf * (1 - exp(-K*(age[a] - a0)));
  }
  for(y=0;y<n_yr;y++) Z(y) = q * eff(y) + M;

  N(0,0,0) = 1.;
  for(a=1;a<astep;a++) N(0,a,0) = N(0,a-1,0) * exp(-Z_init/n_season);
  if(n_season>1) {
    for(k=1;k<n_season;k++) {
      N(0,0,k) = 1.;
      for(a=1;a<astep;a++) N(0,a,k) = N(0,a-1,k-1) * exp(-Z(0)/n_season);
    }
  }
  for(y=1;y<n_yr;y++) {
    N(y,0,0) = 1.;
    for(a=1;a<astep;a++) N(y,a,0) = N(y-1,a-1,n_season-1) * exp(-Z(y)/n_season);
    if(n_season>1) {
      for(k=1;k<n_season;k++) {
        N(y,0,k) = 1.;
        for(a=1;a<astep;a++) N(y,a,k) = N(y,a-1,k-1) * exp(-Z(y)/n_season);
      }
    }
  }
  for(y=0;y<n_yr;y++) {
    for(a=0;a<astep;a++) Nobs(y,a) = N(y,a,obs_season-1) * exp(-Z(y) * timing/n_season);
    Lpred(y) = sum(Nobs(y,_)*La)/sum(Nobs(y,_));

    if(ss[y]>0) {
      ndata += 1.;
      sum_square += ss[y] * pow(Lbar[y]- Lpred[y], 2);
    }
  }
  sigma = pow(sum_square/ndata, 0.5);
  nLL = -ndata * log(sigma) - 0.5 * sum_square/(sigma*sigma);
  nLL *= -1;
  return nLL;
}


// [[Rcpp::export]]
List MLepred(NumericVector stpar, NumericVector Lbar, NumericVector ss,
             NumericVector eff, NumericVector LH, double Lc, double eff_init,
             int n_age, int n_season, int obs_season, double timing, int logpar) {
  double Linf = LH[0];
  double K = LH[1];
  double a0 = LH[2];
  double ac = a0 - log(1 - Lc/Linf)/K;

  int n_yr = Lbar.size();
  int y;
  int a;
  int k;
  double ndata = 0.;
  int astep = n_age*n_season;

  double q;
  double M;

  if(logpar == 0) {
    q = stpar[0];
    M = stpar[1];
  }
  if(logpar == 1) {
    q = exp(stpar[0]);
    M = exp(stpar[1]);
  }

  double Z_init = q * eff_init + M;
  NumericVector Z(n_yr);
  arma::cube N(n_yr,astep,n_season);
  NumericMatrix Nobs(n_yr,astep);
  NumericVector La(astep);
  NumericVector age(astep);

  NumericVector Lpred(n_yr);
  double sum_square = 0.;
  double sigma;

  for(a=0;a<astep;a++) {
    double ageD = a;
    double seasD = n_season;
    age[a] = ac + ageD/seasD;
    La[a] = Linf * (1 - exp(-K*(age[a] - a0)));
  }
  for(y=0;y<n_yr;y++) Z(y) = q * eff(y) + M;

  N(0,0,0) = 1.;
  for(a=1;a<astep;a++) N(0,a,0) = N(0,a-1,0) * exp(-Z_init/n_season);
  if(n_season>1) {
    for(k=1;k<n_season;k++) {
      N(0,0,k) = 1.;
      for(a=1;a<astep;a++) N(0,a,k) = N(0,a-1,k-1) * exp(-Z(0)/n_season);
    }
  }
  for(y=1;y<n_yr;y++) {
    N(y,0,0) = 1.;
    for(a=1;a<astep;a++) N(y,a,0) = N(y-1,a-1,n_season-1) * exp(-Z(y)/n_season);
    if(n_season>1) {
      for(k=1;k<n_season;k++) {
        N(y,0,k) = 1.;
        for(a=1;a<astep;a++) N(y,a,k) = N(y,a-1,k-1) * exp(-Z(y)/n_season);
      }
    }
  }
  for(y=0;y<n_yr;y++) {
    for(a=0;a<astep;a++) Nobs(y,a) = N(y,a,obs_season-1) * exp(-Z(y) * timing/n_season);
    Lpred(y) = sum(Nobs(y,_)*La)/sum(Nobs(y,_));

    if(ss[y]>0) {
      ndata += 1.;
      sum_square += ss[y] * pow(Lbar[y]- Lpred[y], 2);
    }
  }
  sigma = pow(sum_square/ndata, 0.5);
  return List::create(Lpred, sigma);
}


// [[Rcpp::export]]
double MLefullnegLL(NumericVector stpar, NumericVector Lbar, NumericVector ss,
                    NumericVector eff, NumericVector LH, double Lc, double eff_init,
                    int n_age, int n_season, int obs_season, double timing, int logpar) {

  double Linf = LH[0];
  double K = LH[1];
  double a0 = LH[2];
  double ac = a0 - log(1 - Lc/Linf)/K;

  int n_yr = Lbar.size();
  int y;
  int a;
  int k;
  double ndata = 0.;
  int astep = n_age*n_season;

  double q;
  double M;
  double sigma = stpar[2];

  if(logpar == 0) {
    q = stpar[0];
    M = stpar[1];
  }
  if(logpar == 1) {
    q = exp(stpar[0]);
    M = exp(stpar[1]);
  }

  double Z_init = q * eff_init + M;
  NumericVector Z(n_yr);
  arma::cube N(n_yr,astep,n_season);
  NumericMatrix Nobs(n_yr,astep);
  NumericVector La(astep);
  NumericVector age(astep);

  NumericVector Lpred(n_yr);
  double nLL = 0.;

  for(a=0;a<astep;a++) {
    double ageD = a;
    double seasD = n_season;
    age[a] = ac + ageD/seasD;
    La[a] = Linf * (1 - exp(-K*(age[a] - a0)));
  }
  for(y=0;y<n_yr;y++) Z(y) = q * eff(y) + M;

  N(0,0,0) = 1.;
  for(a=1;a<astep;a++) N(0,a,0) = N(0,a-1,0) * exp(-Z_init/n_season);
  if(n_season>1) {
    for(k=1;k<n_season;k++) {
      N(0,0,k) = 1.;
      for(a=1;a<astep;a++) N(0,a,k) = N(0,a-1,k-1) * exp(-Z(0)/n_season);
    }
  }
  for(y=1;y<n_yr;y++) {
    N(y,0,0) = 1.;
    for(a=1;a<astep;a++) N(y,a,0) = N(y-1,a-1,n_season-1) * exp(-Z(y)/n_season);
    if(n_season>1) {
      for(k=1;k<n_season;k++) {
        N(y,0,k) = 1.;
        for(a=1;a<astep;a++) N(y,a,k) = N(y,a-1,k-1) * exp(-Z(y)/n_season);
      }
    }
  }
  for(y=0;y<n_yr;y++) {
    for(a=0;a<astep;a++) Nobs(y,a) = N(y,a,obs_season-1) * exp(-Z(y) * timing/n_season);
    Lpred(y) = sum(Nobs(y,_)*La)/sum(Nobs(y,_));

    nLL += -log(sigma) - 0.5 * ss[y] * pow(Lbar[y] - Lpred[y], 2)/(sigma*sigma);
  }

  nLL *= -1;
  return nLL;
}


// [[Rcpp::export]]
double MLefixM(NumericVector stpar, NumericVector Lbar, NumericVector ss,
               NumericVector eff, NumericVector LH, double Lc, double eff_init,
               int n_age, int n_season, int obs_season, double timing, int logpar) {
  double Linf = LH[0];
  double K = LH[1];
  double a0 = LH[2];
  double ac = a0 - log(1 - Lc/Linf)/K;

  int n_yr = Lbar.size();
  int y;
  int a;
  int k;
  double ndata = 0.;
  int astep = n_age*n_season;

  double q;
  double M = LH[3];

  if(logpar == 0) q = stpar[0];
  if(logpar == 1) q = exp(stpar[0]);

  double Z_init = q * eff_init + M;
  NumericVector Z(n_yr);
  arma::cube N(n_yr,astep,n_season);
  NumericMatrix Nobs(n_yr,astep);
  NumericVector La(astep);
  NumericVector age(astep);

  NumericVector Lpred(n_yr);
  double sum_square = 0.;
  double sigma;
  double nLL;

  for(a=0;a<astep;a++) {
    double ageD = a;
    double seasD = n_season;
    age[a] = ac + ageD/seasD;
    La[a] = Linf * (1 - exp(-K*(age[a] - a0)));
  }
  for(y=0;y<n_yr;y++) Z(y) = q * eff(y) + M;

  N(0,0,0) = 1.;
  for(a=1;a<astep;a++) N(0,a,0) = N(0,a-1,0) * exp(-Z_init/n_season);
  if(n_season>1) {
    for(k=1;k<n_season;k++) {
      N(0,0,k) = 1.;
      for(a=1;a<astep;a++) N(0,a,k) = N(0,a-1,k-1) * exp(-Z(0)/n_season);
    }
  }
  for(y=1;y<n_yr;y++) {
    N(y,0,0) = 1.;
    for(a=1;a<astep;a++) N(y,a,0) = N(y-1,a-1,n_season-1) * exp(-Z(y)/n_season);
    if(n_season>1) {
      for(k=1;k<n_season;k++) {
        N(y,0,k) = 1.;
        for(a=1;a<astep;a++) N(y,a,k) = N(y,a-1,k-1) * exp(-Z(y)/n_season);
      }
    }
  }
  for(y=0;y<n_yr;y++) {
    for(a=0;a<astep;a++) Nobs(y,a) = N(y,a,obs_season-1) * exp(-Z(y) * timing/n_season);
    Lpred(y) = sum(Nobs(y,_)*La)/sum(Nobs(y,_));

    if(ss[y]>0) {
      ndata += 1.;
      sum_square += ss[y] * pow(Lbar[y]- Lpred[y], 2);
    }
  }
  sigma = pow(sum_square/ndata, 0.5);
  nLL = -ndata * log(sigma) - 0.5 * sum_square/(sigma*sigma);
  nLL *= -1;
  return nLL;
}


// [[Rcpp::export]]
double MLefixMfullnegLL(NumericVector stpar, NumericVector Lbar, NumericVector ss,
                        NumericVector eff, NumericVector LH, double Lc, double eff_init,
                        int n_age, int n_season, int obs_season, double timing, int logpar) {

  double Linf = LH[0];
  double K = LH[1];
  double a0 = LH[2];
  double ac = a0 - log(1 - Lc/Linf)/K;

  int n_yr = Lbar.size();
  int y;
  int a;
  int k;
  double ndata = 0.;
  int astep = n_age*n_season;

  double q;
  double M = LH[3];
  double sigma = stpar[1];

  if(logpar == 0) q = stpar[0];
  if(logpar == 1) q = exp(stpar[0]);

  double Z_init = q * eff_init + M;
  NumericVector Z(n_yr);
  arma::cube N(n_yr,astep,n_season);
  NumericMatrix Nobs(n_yr,astep);
  NumericVector La(astep);
  NumericVector age(astep);

  NumericVector Lpred(n_yr);
  double nLL = 0.;

  for(a=0;a<astep;a++) {
    double ageD = a;
    double seasD = n_season;
    age[a] = ac + ageD/seasD;
    La[a] = Linf * (1 - exp(-K*(age[a] - a0)));
  }
  for(y=0;y<n_yr;y++) Z(y) = q * eff(y) + M;

  N(0,0,0) = 1.;
  for(a=1;a<astep;a++) N(0,a,0) = N(0,a-1,0) * exp(-Z_init/n_season);
  if(n_season>1) {
    for(k=1;k<n_season;k++) {
      N(0,0,k) = 1.;
      for(a=1;a<astep;a++) N(0,a,k) = N(0,a-1,k-1) * exp(-Z(0)/n_season);
    }
  }
  for(y=1;y<n_yr;y++) {
    N(y,0,0) = 1.;
    for(a=1;a<astep;a++) N(y,a,0) = N(y-1,a-1,n_season-1) * exp(-Z(y)/n_season);
    if(n_season>1) {
      for(k=1;k<n_season;k++) {
        N(y,0,k) = 1.;
        for(a=1;a<astep;a++) N(y,a,k) = N(y,a-1,k-1) * exp(-Z(y)/n_season);
      }
    }
  }
  for(y=0;y<n_yr;y++) {
    for(a=0;a<astep;a++) Nobs(y,a) = N(y,a,obs_season-1) * exp(-Z(y) * timing/n_season);
    Lpred(y) = sum(Nobs(y,_)*La)/sum(Nobs(y,_));

    nLL += -log(sigma) - 0.5 * ss[y] * pow(Lbar[y] - Lpred[y], 2)/(sigma*sigma);
  }

  nLL *= -1;
  return nLL;
}
