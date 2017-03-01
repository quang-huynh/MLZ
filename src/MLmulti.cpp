#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double SSM_MSM1_negLL(NumericVector stpar, NumericMatrix Lbar, NumericMatrix ss,
                NumericMatrix LH, NumericVector Lc, int nbreaks, int isMSM1) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericMatrix yearZ(nspec,nbreaks);
  arma::cube dy(nspec,nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix denom(nspec,count);
  NumericMatrix numsum(nspec,count);
  NumericMatrix num(nspec,count);
  NumericMatrix Lpred(nspec,count);

  NumericVector sum_square_Lpred(nspec);
  NumericVector nyear(nspec);
  NumericVector sigma(nspec);
  double nLL = 0.;

  if(nbr == 0 & isMSM1 == 0) yearZ(_,0) = stpar[Range(nspec*(nbreaks+1),nspec*(nbreaks+1)+nspec-1)];
  for(w=0;w<nspec;w++) {
    Z(w,_) = stpar[Range(w*(nbreaks+1),w*(nbreaks+1)+nbreaks)];
    if(nbr == 0 & isMSM1 == 1) yearZ(w,0) = stpar[nspec*(nbreaks+1)];
    if(nbr>0) {
      if(isMSM1 == 0) yearZ(w,_) = stpar[Range(nspec*(nbreaks+1)+w*nbreaks,nspec*(nbreaks+1)+w*nbreaks+nbreaks-1)];
      if(isMSM1 == 1) yearZ(w,_) = stpar[Range(nspec*(nbreaks+1),nspec*(nbreaks+1)+nbreaks-1)];
    }
    for(i=0;i<=nbr;i++) {
      for(m=1;m<=count;m++) {
        if(yearZ(w,i)>=m) dy(w,i,m-1) = 0.;
        else dy(w,i,m-1) = m-yearZ(w,i);
      }
    }
    if(nbreaks>1) {
      for(i=0;i<nbr;i++){
        for(m=0;m<count;m++) dy(w,i,m) -= dy(w,i+1,m);
      }
    }
  }
  for(w=0;w<nspec;w++) {
    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(w,nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(w,nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(w,nbr-j,m));
          }
        }

        if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(w,nbr-i,m)))/Z(w,nbr+1-i);
        if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);

        numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
      }

      num(w,m) = Linf(w) * (denom(w,m) - (1. - Lc(w)/Linf(w)) * numsum(w,m));
      Lpred(w,m) = num(w,m)/denom(w,m);

      if(ss(w,m)>0) {
        sum_square_Lpred(w) += ss(w,m) * pow(Lbar(w,m) - Lpred(w,m), 2);
        nyear(w) += 1.;
      }
    }
    sigma(w) = sqrt(sum_square_Lpred(w)/nyear(w));
    nLL += -nyear(w) * log(sigma(w)) - 0.5 * sum_square_Lpred(w)/pow(sigma(w), 2);
  }
  nLL *= -1;

  return nLL;
}


// [[Rcpp::export]]
double SSM_MSM1_profile(NumericVector stpar, NumericVector year, NumericMatrix Lbar, NumericMatrix ss,
                        NumericMatrix LH, NumericVector Lc, int nbreaks) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericMatrix yearZ(nspec,nbreaks);
  arma::cube dy(nspec,nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix denom(nspec,count);
  NumericMatrix numsum(nspec,count);
  NumericMatrix num(nspec,count);
  NumericMatrix Lpred(nspec,count);

  NumericVector sum_square_Lpred(nspec);
  NumericVector nyear(nspec);
  NumericVector sigma(nspec);
  double nLL = 0.;

  for(w=0;w<nspec;w++) {
    Z(w,_) = stpar[Range(w*(nbreaks+1),w*(nbreaks+1)+nbreaks)];
    yearZ(w,_) = year;

    for(i=0;i<=nbr;i++) {
      for(m=1;m<=count;m++) {
        if(yearZ(w,i)>=m) dy(w,i,m-1) = 0.;
        else dy(w,i,m-1) = m-yearZ(w,i);
      }
    }
    if(nbreaks>1) {
      for(i=0;i<nbr;i++){
        for(m=0;m<count;m++) dy(w,i,m) -= dy(w,i+1,m);
      }
    }
  }
  for(w=0;w<nspec;w++) {
    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(w,nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(w,nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(w,nbr-j,m));
          }
        }

        if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(w,nbr-i,m)))/Z(w,nbr+1-i);
        if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);

        numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
      }

      num(w,m) = Linf(w) * (denom(w,m) - (1. - Lc(w)/Linf(w)) * numsum(w,m));
      Lpred(w,m) = num(w,m)/denom(w,m);

      if(ss(w,m)>0) {
        sum_square_Lpred(w) += ss(w,m) * pow(Lbar(w,m) - Lpred(w,m), 2);
        nyear(w) += 1.;
      }
    }
    sigma(w) = sqrt(sum_square_Lpred(w)/nyear(w));
    nLL += -nyear(w) * log(sigma(w)) - 0.5 * sum_square_Lpred(w)/pow(sigma(w), 2);
  }
  nLL *= -1;

  return nLL;
}


// [[Rcpp::export]]
List SSM_MSM1_pred(NumericVector stpar, NumericMatrix Lbar, NumericMatrix ss,
                      NumericMatrix LH, NumericVector Lc, int nbreaks, int isMSM1) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericMatrix yearZ(nspec,nbreaks);
  arma::cube dy(nspec,nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix denom(nspec,count);
  NumericMatrix numsum(nspec,count);
  NumericMatrix num(nspec,count);
  NumericMatrix Lpred(nspec,count);

  NumericVector sum_square_Lpred(nspec);
  NumericVector nyear(nspec);
  NumericVector sigma(nspec);
  NumericVector loglike(nspec);

  if(nbr == 0 & isMSM1 == 0) yearZ(_,0) = stpar[Range(nspec*(nbreaks+1),nspec*(nbreaks+1)+nspec-1)];
  for(w=0;w<nspec;w++) {
    Z(w,_) = stpar[Range(w*(nbreaks+1),w*(nbreaks+1)+nbreaks)];
    if(nbr == 0 & isMSM1 == 1) yearZ(w,0) = stpar[nspec*(nbreaks+1)];
    if(nbr>0) {
      if(isMSM1 == 0) yearZ(w,_) = stpar[Range(nspec*(nbreaks+1)+w*nbreaks,nspec*(nbreaks+1)+w*nbreaks+nbreaks-1)];
      if(isMSM1 == 1) yearZ(w,_) = stpar[Range(nspec*(nbreaks+1),nspec*(nbreaks+1)+nbreaks-1)];
    }
    for(i=0;i<=nbr;i++) {
      for(m=1;m<=count;m++) {
        if(yearZ(w,i)>=m) dy(w,i,m-1) = 0.;
        else dy(w,i,m-1) = m-yearZ(w,i);
      }
    }
    if(nbreaks>1) {
      for(i=0;i<nbr;i++){
        for(m=0;m<count;m++) dy(w,i,m) -= dy(w,i+1,m);
      }
    }
  }
  for(w=0;w<nspec;w++) {
    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(w,nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(w,nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(w,nbr-j,m));
          }
        }

        if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(w,nbr-i,m)))/Z(w,nbr+1-i);
        if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);

        numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
      }

      num(w,m) = Linf(w) * (denom(w,m) - (1. - Lc(w)/Linf(w)) * numsum(w,m));
      Lpred(w,m) = num(w,m)/denom(w,m);

      if(ss(w,m)>0) {
        sum_square_Lpred(w) += ss(w,m) * pow(Lbar(w,m) - Lpred(w,m), 2);
        nyear(w) += 1.;
      }
    }
    sigma(w) = sqrt(sum_square_Lpred(w)/nyear(w));
    loglike(w) = -nyear(w) * log(sigma(w)) - 0.5 * sum_square_Lpred(w)/pow(sigma(w), 2);
  }

  return List::create(Lpred, sigma, loglike);
}


// [[Rcpp::export]]
double SSM_MSM1_fullnegLL(NumericVector stpar, NumericMatrix Lbar, NumericMatrix ss,
                          NumericMatrix LH, NumericVector Lc, int nbreaks, int isMSM1) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericMatrix yearZ(nspec,nbreaks);
  arma::cube dy(nspec,nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix denom(nspec,count);
  NumericMatrix numsum(nspec,count);
  NumericMatrix num(nspec,count);
  NumericMatrix Lpred(nspec,count);

  NumericVector sigma(nspec);
  double nLL = 0.;

  if(nbr == 0 & isMSM1 == 0) {
    yearZ(_,0) = stpar[Range(nspec*(nbreaks+1),nspec*(nbreaks+1)+nspec-1)];
    sigma = stpar[Range(nspec*(nbreaks+1)+nspec, nspec*(nbreaks+1)+nspec+nspec-1)];
  }

  for(w=0;w<nspec;w++) {
    Z(w,_) = stpar[Range(w*(nbreaks+1),w*(nbreaks+1)+nbreaks)];
    if(nbr == 0 & isMSM1 == 1) {
      yearZ(w,0) = stpar[nspec*(nbreaks+1)];
      sigma(w) = stpar[nspec*(nbreaks+1)+1+w];
    }
    if(nbr>0) {
      if(isMSM1 == 0) {
        yearZ(w,_) = stpar[Range(nspec*(nbreaks+1)+w*nbreaks,nspec*(nbreaks+1)+w*nbreaks+nbreaks-1)];
        sigma(w) = stpar[nspec*(nbreaks+1)+w*nbreaks+nbreaks+w];
      }
      if(isMSM1 == 1) {
        yearZ(w,_) = stpar[Range(nspec*(nbreaks+1),nspec*(nbreaks+1)+nbreaks-1)];
        sigma(w) = stpar[nspec*(nbreaks+1)+nbreaks+w];
      }
    }
    for(i=0;i<=nbr;i++) {
      for(m=1;m<=count;m++) {
        if(yearZ(w,i)>=m) dy(w,i,m-1) = 0.;
        else dy(w,i,m-1) = m-yearZ(w,i);
      }
    }
    if(nbreaks>1) {
      for(i=0;i<nbr;i++){
        for(m=0;m<count;m++) dy(w,i,m) -= dy(w,i+1,m);
      }
    }
  }

  for(w=0;w<nspec;w++) {
    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(w,nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(w,nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(w,nbr-j,m));
          }
        }

        if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(w,nbr-i,m)))/Z(w,nbr+1-i);
        if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);

        numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
      }

      num(w,m) = Linf(w) * (denom(w,m) - (1. - Lc(w)/Linf(w)) * numsum(w,m));
      Lpred(w,m) = num(w,m)/denom(w,m);

      if(ss(w,m)>0) nLL += -log(sigma(w)) - 0.5 * ss(w,m) * pow(Lbar(w,m) - Lpred(w,m), 2)/pow(sigma(w), 2);
    }
  }
  nLL *= -1;
  return nLL;
}
