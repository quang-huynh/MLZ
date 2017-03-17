#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MSM23_negLL(NumericVector stpar, NumericMatrix Lbar, NumericMatrix ss,
                   NumericMatrix LH, NumericVector Lc, int nbreaks, int isMSM3, int spCont) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);
  NumericVector M = LH(_,2);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericVector yearZ(nbreaks);
  NumericMatrix dy(nbreaks,count);

  NumericVector Z1(nspec);
  NumericVector delta(nbreaks);
  NumericVector epsilon(nspec-1);

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

  Z1 = stpar[Range(0,nspec-1)];
  delta = stpar[Range(nspec,nspec+nbr)];
  if(isMSM3 == 0) {
    epsilon = stpar[Range(nspec+nbr+1,nspec+nbr+1+nspec-2)];
    yearZ = stpar[Range(nspec+nbr+1+nspec-1,nspec+nbr+1+nspec-1+nbr)];
  }
  if(isMSM3 == 1) yearZ = stpar[Range(nspec+nbr+1,nspec+nbr+1+nbr)];

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ(i)>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ(i);
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(w=0;w<nspec;w++) {
    Z(w,0) = Z1(w);
    if(w<nspec-1 & isMSM3 == 1) epsilon(w) = 1.;
    for(i=1;i<=nbr+1;i++) {
      if(w==0) Z(w,i) = delta(i-1) * Z(w,i-1) + (1-delta(i-1)) * M(w);
      if(w>0) Z(w,i) = delta(i-1) * epsilon(w-1) * Z(w,i-1) + (1-delta(i-1)*epsilon(w-1)) * M(w);
    }

    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(nbr-j,m));
          }
        }
        
        if(spCont == 1) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/Z(w,nbr+1-i);
          if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);
          numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
        }
        if(spCont == 0) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/(1 - exp(-Z(w,nbr+1-i)));
          if(i==nbr+1) denom(w,m) += a(i,m)/(1 - exp(-Z(w,nbr+1-i)));
          numsum(w,m) += r(i,m) * s(i,m) / (1 - exp(-(Z(w,nbr+1-i) + K(w))));
        }
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
double MSM23_profile(NumericVector stpar, NumericVector year, NumericMatrix Lbar, NumericMatrix ss,
                     NumericMatrix LH, NumericVector Lc, int nbreaks, int isMSM3, int spCont) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);
  NumericVector M = LH(_,2);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericVector yearZ(nbreaks);
  NumericMatrix dy(nbreaks,count);

  NumericVector Z1(nspec);
  NumericVector delta(nbreaks);
  NumericVector epsilon(nspec-1);

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

  Z1 = stpar[Range(0,nspec-1)];
  delta = stpar[Range(nspec,nspec+nbr)];
  if(isMSM3 == 0) epsilon = stpar[Range(nspec+nbr+1,nspec+nbr+1+nspec-2)];
  yearZ = year;

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ(i)>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ(i);
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(w=0;w<nspec;w++) {
    Z(w,0) = Z1(w);
    if(w<nspec-1 & isMSM3 == 1) epsilon(w) = 1.;
    for(i=1;i<=nbr+1;i++) {
      if(w==0) Z(w,i) = delta(i-1) * Z(w,i-1) + (1-delta(i-1)) * M(w);
      if(w>0) Z(w,i) = delta(i-1) * epsilon(w-1) * Z(w,i-1) + (1-delta(i-1)*epsilon(w-1)) * M(w);
    }

    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(nbr-j,m));
          }
        }
        
        if(spCont == 1) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/Z(w,nbr+1-i);
          if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);
          numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
        }
        if(spCont == 0) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/(1 - exp(-Z(w,nbr+1-i)));
          if(i==nbr+1) denom(w,m) += a(i,m)/(1 - exp(-Z(w,nbr+1-i)));
          numsum(w,m) += r(i,m) * s(i,m) / (1 - exp(-(Z(w,nbr+1-i) + K(w))));
        }
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
List MSM23_pred(NumericVector stpar, NumericMatrix Lbar, NumericMatrix ss,
                NumericMatrix LH, NumericVector Lc, int nbreaks, int isMSM3, int spCont) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);
  NumericVector M = LH(_,2);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericVector yearZ(nbreaks);
  NumericMatrix dy(nbreaks,count);

  NumericVector Z1(nspec);
  NumericVector delta(nbreaks);
  NumericVector epsilon(nspec-1);

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
  NumericVector loglikesp(nspec);

  Z1 = stpar[Range(0,nspec-1)];
  delta = stpar[Range(nspec,nspec+nbr)];
  if(isMSM3 == 0) {
    epsilon = stpar[Range(nspec+nbr+1,nspec+nbr+1+nspec-2)];
    yearZ = stpar[Range(nspec+nbr+1+nspec-1,nspec+nbr+1+nspec-1+nbr)];
  }
  if(isMSM3 == 1) yearZ = stpar[Range(nspec+nbr+1,nspec+nbr+1+nbr)];

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ(i)>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ(i);
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(w=0;w<nspec;w++) {
    Z(w,0) = Z1(w);
    if(w<nspec-1 & isMSM3 == 1) epsilon(w) = 1.;
    for(i=1;i<=nbr+1;i++) {
      if(w==0) Z(w,i) = delta(i-1) * Z(w,i-1) + (1-delta(i-1)) * M(w);
      if(w>0) Z(w,i) = delta(i-1) * epsilon(w-1) * Z(w,i-1) + (1-delta(i-1)*epsilon(w-1)) * M(w);
    }

    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(nbr-j,m));
          }
        }
        
        if(spCont == 1) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/Z(w,nbr+1-i);
          if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);
          numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
        }
        if(spCont == 0) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/(1 - exp(-Z(w,nbr+1-i)));
          if(i==nbr+1) denom(w,m) += a(i,m)/(1 - exp(-Z(w,nbr+1-i)));
          numsum(w,m) += r(i,m) * s(i,m) / (1 - exp(-(Z(w,nbr+1-i) + K(w))));
        }
      }

      num(w,m) = Linf(w) * (denom(w,m) - (1. - Lc(w)/Linf(w)) * numsum(w,m));
      Lpred(w,m) = num(w,m)/denom(w,m);

      if(ss(w,m)>0) {
        sum_square_Lpred(w) += ss(w,m) * pow(Lbar(w,m) - Lpred(w,m), 2);
        nyear(w) += 1.;
      }
    }
    sigma(w) = sqrt(sum_square_Lpred(w)/nyear(w));
    loglikesp(w) = -nyear(w) * log(sigma(w)) - 0.5 * sum_square_Lpred(w)/pow(sigma(w), 2);
  }
  return List::create(Lpred, sigma, loglikesp, Z);
}



// [[Rcpp::export]]
double MSM23_fullnegLL(NumericVector stpar, NumericMatrix Lbar, NumericMatrix ss,
                       NumericMatrix LH, NumericVector Lc, int nbreaks, int isMSM3, int spCont) {

  int i;
  int j;
  int m;
  int w;

  int count = Lbar.ncol();
  int nspec = Lbar.nrow();
  int nbr = nbreaks-1;

  NumericVector Linf = LH(_,0);
  NumericVector K = LH(_,1);
  NumericVector M = LH(_,2);

  NumericMatrix Z(nspec,nbreaks+1);
  NumericVector yearZ(nbreaks);
  NumericMatrix dy(nbreaks,count);

  NumericVector Z1(nspec);
  NumericVector delta(nbreaks);
  NumericVector epsilon(nspec-1);
  NumericVector sigma(nspec);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix denom(nspec,count);
  NumericMatrix numsum(nspec,count);
  NumericMatrix num(nspec,count);
  NumericMatrix Lpred(nspec,count);

  double nLL = 0.;

  Z1 = stpar[Range(0,nspec-1)];
  delta = stpar[Range(nspec,nspec+nbr)];
  if(isMSM3 == 0) {
    epsilon = stpar[Range(nspec+nbr+1,nspec+nbr+1+nspec-2)];
    yearZ = stpar[Range(nspec+nbr+1+nspec-1,nspec+nbr+1+nspec-1+nbr)];
    sigma = stpar[Range(nspec+nbr+1+nspec-1+nbr+1,nspec+nbr+1+nspec-1+nbr+nspec)];
  }
  if(isMSM3 == 1) {
    yearZ = stpar[Range(nspec+nbr+1,nspec+nbr+1+nbr)];
    sigma = stpar[Range(nspec+nbr+1+nbr+1,nspec+nbr+1+nbr+nspec)];
  }
  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ(i)>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ(i);
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(w=0;w<nspec;w++) {
    Z(w,0) = Z1(w);
    if(w<nspec-1 & isMSM3 == 1) epsilon(w) = 1.;
    for(i=1;i<=nbr+1;i++) {
      if(w==0) Z(w,i) = delta(i-1) * Z(w,i-1) + (1-delta(i-1)) * M(w);
      if(w>0) Z(w,i) = delta(i-1) * epsilon(w-1) * Z(w,i-1) + (1-delta(i-1)*epsilon(w-1)) * M(w);
    }

    for(m=0;m<count;m++) {
      denom(w,m) = 0.;
      numsum(w,m) = 0.;

      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;

        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(w,nbr+1-i) + K(w)) * dy(nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;

        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(w,nbr+1-j) * dy(nbr-j,m));
            r(i,m) *= exp(-(Z(w,nbr+1-j) + K(w)) * dy(nbr-j,m));
          }
        }
        
        if(spCont == 1) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/Z(w,nbr+1-i);
          if(i==nbr+1) denom(w,m) += a(i,m)/Z(w,nbr+1-i);
          numsum(w,m) += r(i,m) * s(i,m) / (Z(w,nbr+1-i) + K(w));
        }
        if(spCont == 0) {
          if(i<=nbr) denom(w,m) += a(i,m)*(1. - exp(-Z(w,nbr+1-i) * dy(nbr-i,m)))/(1 - exp(-Z(w,nbr+1-i)));
          if(i==nbr+1) denom(w,m) += a(i,m)/(1 - exp(-Z(w,nbr+1-i)));
          numsum(w,m) += r(i,m) * s(i,m) / (1 - exp(-(Z(w,nbr+1-i) + K(w))));
        }
      }

      num(w,m) = Linf(w) * (denom(w,m) - (1. - Lc(w)/Linf(w)) * numsum(w,m));
      Lpred(w,m) = num(w,m)/denom(w,m);

      if(ss(w,m)>0) nLL += -log(sigma(w)) - 0.5 * ss(w,m) * pow(Lbar(w,m) - Lpred(w,m), 2)/pow(sigma(w), 2);
    }
  }
  nLL *= -1;
  return nLL;
}

