#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MLeqnegLL(double Z, NumericVector Lbar, NumericVector ss,
                 NumericVector LH, double Lc, int spCont) {

  int m;
  int count = Lbar.size();

  double Linf = LH[0];
  double K = LH[1];

  double Lpred;
  
  if(spCont == 1) Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
  if(spCont == 0) Lpred = Linf * (1 - ((1 - exp(-Z))/(1 - exp(-Z-K))) * (1 - Lc/Linf));

  double sum_square_Lpred = 0.;
  double nyear = 0.;
  double sigma;
  double nLL;

  for(m=0;m<count;m++) {
    if(ss[m]>0) {
      sum_square_Lpred += ss[m] * pow(Lbar[m]-Lpred,2.);
      nyear += 1.;
    }
  }
  sigma = sqrt(sum_square_Lpred/nyear);
  nLL = -nyear * log(sigma) - 0.5 * sum_square_Lpred/(sigma*sigma);
  nLL *= -1;
  return nLL;
}

// [[Rcpp::export]]
List MLeqpred(double Z, NumericVector Lbar, NumericVector ss,
              NumericVector LH, double Lc, int spCont) {
  int m;
  int count = Lbar.size();

  double Linf = LH[0];
  double K = LH[1];
  
  double Lpred;
  
  if(spCont == 1) Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
  if(spCont == 0) Lpred = Linf * (1 - ((1 - exp(-Z))/(1 - exp(-Z-K))) * (1 - Lc/Linf));

  double sum_square_Lpred = 0.;
  double nyear = 0.;
  double sigma;
  double nLL;

  for(m=0;m<count;m++) {
    if(ss[m]>0) {
      sum_square_Lpred += ss[m] * pow(Lbar[m]-Lpred,2.);
      nyear += 1.;
    }
  }
  sigma = sqrt(sum_square_Lpred/nyear);

  return List::create(Lpred, sigma);
}



// [[Rcpp::export]]
double MLeqfullnegLL(NumericVector stpar, NumericVector Lbar, NumericVector ss,
                     NumericVector LH, double Lc, int spCont) {
  int m;
  int count = Lbar.size();

  double Linf = LH[0];
  double K = LH[1];

  double Z = stpar[0];
  double sigma = stpar[1];
  double Lpred;
  
  if(spCont == 1) Lpred = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
  if(spCont == 0) Lpred = Linf * (1 - ((1 - exp(-Z))/(1 - exp(-Z-K))) * (1 - Lc/Linf));
  
  double nLL = 0.;

  for(m=0;m<count;m++) {
    if(ss[m]>0) {
      nLL += -log(sigma) - 0.5 * ss[m] * pow(Lbar[m] - Lpred, 2)/(sigma*sigma);
    }
  }
  nLL *= -1;
  return nLL;
}
