#include <Rcpp.h>
using namespace Rcpp;

double pbeta_inc(double x, double a, double b) {
  double ans = R::beta(a,b);
  ans *= R::pbeta(x,a,b,1,0);
  return ans;
}

// [[Rcpp::export]]
double MLWPUEnegLL(NumericVector stpar, NumericVector Lbar, NumericVector ss, NumericVector CPUE,
                   NumericVector LH, double Lc, int nbreaks, int loglikeCPUE, int spCont) {

  int i;
  int j;
  int m;

  int count = Lbar.size();
  int nbr = nbreaks-1;

  double Linf = LH[0];
  double K = LH[1];
  double b = LH[2];

  NumericVector Z(nbreaks+1);
  NumericVector yearZ(nbreaks);
  NumericMatrix dy(nbreaks,count);
  NumericMatrix sumy(nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix w(nbreaks+1,count);
  NumericMatrix int_upper(nbreaks+1,count);
  NumericMatrix int_lower(nbreaks+1,count);
  NumericVector denom(count);
  NumericVector numsum(count);
  NumericVector num(count);
  NumericVector biomass(count);
  NumericVector Lpred(count);
  NumericVector Ipred(count);

  double q;
  NumericVector sum_square(2);
  NumericVector nyear(2);
  NumericVector sigma(2);
  NumericVector loglike(2);
  double nLL;

  for(i=0;i<=nbr;i++) {
    Z[i] = stpar[i];
    yearZ[i] = stpar[nbr+2+i];
  }
  Z[nbr+1] = stpar[nbr+1];

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ[i]>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ[i];
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(m=1;m<=count;m++) {
    for(i=0;i<=nbr;i++) {
      if(yearZ(i)>=m) sumy(i,m-1) = 0.;
      else sumy(i,m-1) = m - yearZ(i);
    }
    for(i=0;i<=nbr+1;i++) {
      if(i==0) {
        if(yearZ(i)>=m) int_lower(i,m-1) = Lc/Linf;
        else int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
        int_upper(i,m-1) = 1.;
      }
      if(i>0 & i<nbr+1) {
        if(yearZ(i-1)<m & yearZ(i)>=m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i)<m) {
          int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i-1)>=m) {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
      if(i==nbr+1) {
        if(yearZ(i-1)<m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        else {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
    }
  }

  for(m=0;m<count;m++) {
    denom[m] = 0.;
    numsum[m] = 0.;
    biomass[m] = 0.;

    for(i=0;i<=nbr+1;i++) {
      a(i,m) = 1.;
      r(i,m) = 1.;
      w(i,m) = 1.;

      if(i<nbr+1) s(i,m) = 1. - exp(-(Z[nbr+1-i]+K) * dy(nbr-i,m));
      if(i==nbr+1) s(i,m) = 1.;

      if(i>0) {
        for(j=0;j<=i-1;j++) {
          a(i,m) *= exp(-Z[nbr+1-j] * dy(nbr-j,m));
          r(i,m) *= exp(-(Z[nbr+1-j] + K) * dy(nbr-j,m));
          w(i,m) *= exp(Z[nbr+1-i] * dy(nbr-j,m));
        }
      }

      if(i<=nbr) denom[m] += a(i,m)*(1. - exp(-Z[nbr+1-i] * dy(nbr-i,m)))/Z[nbr+1-i];
      if(i==nbr+1) denom[m] += a(i,m)/Z[nbr+1-i];

      numsum[m] += r(i,m) * s(i,m) / (Z[nbr+1-i] + K);
      biomass[m] += a(i,m) * w(i,m) * pow(1 - Lc/Linf, -Z[nbr+1-i]/K) *
        (pbeta_inc(int_upper(nbr+1-i,m), b+1, Z[nbr+1-i]/K)-pbeta_inc(int_lower(nbr+1-i,m), b+1, Z[nbr+1-i]/K));
    }
    num[m] = Linf * (denom[m] - (1. - Lc/Linf) * numsum[m]);
    Lpred[m] = num[m]/denom[m];
  }
  for(m=0;m<count;m++) {
    if(ss[m]>0) nyear(0) += 1.;
    if(CPUE[m]>0) nyear(1) += 1.;
  }

  double sum_q = 0.;
  double sum_q2 = 0.;

  if(loglikeCPUE == 0) {
    for(m=0;m<count;m++)
    {
      if(CPUE[m]>0) sum_q += log(CPUE[m]/biomass[m]);
    }
    q = exp(sum_q/nyear(1));
  }
  if(loglikeCPUE == 1) {
    for(m=0;m<count;m++) {
      if(CPUE[m]>0) {
        sum_q += CPUE[m] * biomass[m];
        sum_q2 += pow(biomass[m],2);
      }
    }
    q = sum_q/sum_q2;
  }
  for(m=0;m<count;m++) {
    Ipred[m] = q * biomass[m];
    if(ss[m]>0) sum_square[0] += ss[m] * pow(Lbar[m]-Lpred[m],2.);
    if(CPUE[m]>0) {
      if(loglikeCPUE == 0) sum_square[1] += pow(log(CPUE[m]/Ipred[m]),2);
      if(loglikeCPUE == 1) sum_square[1] += pow(CPUE[m] - Ipred[m], 2);
    }
  }
  for(i=0;i<2;i++) {
    sigma[i] = pow(sum_square[i]/nyear[i], 0.5);
    loglike[i] = -nyear[i] * log(sigma[i]) - 0.5 * sum_square[i]/pow(sigma[i],2);
  }

  nLL = sum(loglike);
  nLL *= -1;

  return nLL;
}


// [[Rcpp::export]]
double MLWPUEprofile(NumericVector Z, NumericVector yearZ, NumericVector Lbar, NumericVector ss,
                     NumericVector CPUE, NumericVector LH, double Lc,
                     int nbreaks, int loglikeCPUE, int spCont) {

  int i;
  int j;
  int m;

  int count = Lbar.size();
  int nbr = nbreaks-1;

  double Linf = LH[0];
  double K = LH[1];
  double b = LH[2];

  NumericMatrix dy(nbreaks,count);
  NumericMatrix sumy(nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix w(nbreaks+1,count);
  NumericMatrix int_upper(nbreaks+1,count);
  NumericMatrix int_lower(nbreaks+1,count);
  NumericVector denom(count);
  NumericVector numsum(count);
  NumericVector num(count);
  NumericVector biomass(count);
  NumericVector Lpred(count);
  NumericVector Ipred(count);

  double q;
  NumericVector sum_square(2);
  NumericVector nyear(2);
  NumericVector sigma(2);
  NumericVector loglike(2);
  double nLL;

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ[i]>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ[i];
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(m=1;m<=count;m++) {
    for(i=0;i<=nbr;i++) {
      if(yearZ(i)>=m) sumy(i,m-1) = 0.;
      else sumy(i,m-1) = m - yearZ(i);
    }
    for(i=0;i<=nbr+1;i++) {
      if(i==0) {
        if(yearZ(i)>=m) int_lower(i,m-1) = Lc/Linf;
        else int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
        int_upper(i,m-1) = 1.;
      }
      if(i>0 & i<nbr+1) {
        if(yearZ(i-1)<m & yearZ(i)>=m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i)<m) {
          int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i-1)>=m) {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
      if(i==nbr+1) {
        if(yearZ(i-1)<m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        else {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
    }
  }

  for(m=0;m<count;m++) {
    denom[m] = 0.;
    numsum[m] = 0.;
    biomass[m] = 0.;

    for(i=0;i<=nbr+1;i++) {
      a(i,m) = 1.;
      r(i,m) = 1.;
      w(i,m) = 1.;

      if(i<nbr+1) s(i,m) = 1. - exp(-(Z[nbr+1-i]+K) * dy(nbr-i,m));
      if(i==nbr+1) s(i,m) = 1.;

      if(i>0) {
        for(j=0;j<=i-1;j++) {
          a(i,m) *= exp(-Z[nbr+1-j] * dy(nbr-j,m));
          r(i,m) *= exp(-(Z[nbr+1-j] + K) * dy(nbr-j,m));
          w(i,m) *= exp(Z[nbr+1-i] * dy(nbr-j,m));
        }
      }

      if(i<=nbr) denom[m] += a(i,m)*(1. - exp(-Z[nbr+1-i] * dy(nbr-i,m)))/Z[nbr+1-i];
      if(i==nbr+1) denom[m] += a(i,m)/Z[nbr+1-i];

      numsum[m] += r(i,m) * s(i,m) / (Z[nbr+1-i] + K);
      biomass[m] += a(i,m) * w(i,m) * pow(1 - Lc/Linf, -Z[nbr+1-i]/K) *
        (pbeta_inc(int_upper(nbr+1-i,m), b+1, Z[nbr+1-i]/K)-pbeta_inc(int_lower(nbr+1-i,m), b+1, Z[nbr+1-i]/K));
    }
    num[m] = Linf * (denom[m] - (1. - Lc/Linf) * numsum[m]);
    Lpred[m] = num[m]/denom[m];
  }
  for(m=0;m<count;m++) {
    if(ss[m]>0) nyear(0) += 1.;
    if(CPUE[m]>0) nyear(1) += 1.;
  }

  double sum_q = 0.;
  double sum_q2 = 0.;

  if(loglikeCPUE == 0) {
    for(m=0;m<count;m++)
    {
      if(CPUE[m]>0) sum_q += log(CPUE[m]/biomass[m]);
    }
    q = exp(sum_q/nyear(1));
  }
  if(loglikeCPUE == 1) {
    for(m=0;m<count;m++) {
      if(CPUE[m]>0) {
        sum_q += CPUE[m] * biomass[m];
        sum_q2 += pow(biomass[m],2);
      }
    }
    q = sum_q/sum_q2;
  }
  for(m=0;m<count;m++) {
    Ipred[m] = q * biomass[m];
    if(ss[m]>0) sum_square[0] += ss[m] * pow(Lbar[m]-Lpred[m],2.);
    if(CPUE[m]>0) {
      if(loglikeCPUE == 0) sum_square[1] += pow(log(CPUE[m]/Ipred[m]),2);
      if(loglikeCPUE == 1) sum_square[1] += pow(CPUE[m] - Ipred[m], 2);
    }
  }
  for(i=0;i<2;i++) {
    sigma[i] = pow(sum_square[i]/nyear[i], 0.5);
    loglike[i] = -nyear[i] * log(sigma[i]) - 0.5 * sum_square[i]/pow(sigma[i],2);
  }

  nLL = sum(loglike);
  nLL *= -1;

  return nLL;
}

// [[Rcpp::export]]
List MLWPUEpred(NumericVector stpar, NumericVector Lbar, NumericVector ss, NumericVector CPUE,
                NumericVector LH, double Lc, int nbreaks, int loglikeCPUE, int spCont) {

  int i;
  int j;
  int m;

  int count = CPUE.size();
  int nbr = nbreaks-1;

  double Linf = LH[0];
  double K = LH[1];
  double b = LH[2];

  NumericVector Z(nbreaks+1);
  NumericVector yearZ(nbreaks);
  NumericMatrix dy(nbreaks,count);
  NumericMatrix sumy(nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix w(nbreaks+1,count);
  NumericMatrix int_upper(nbreaks+1,count);
  NumericMatrix int_lower(nbreaks+1,count);
  NumericVector denom(count);
  NumericVector numsum(count);
  NumericVector num(count);
  NumericVector biomass(count);
  NumericVector Lpred(count);
  NumericVector Ipred(count);

  double q;
  NumericVector sum_square(2);
  NumericVector nyear(2);
  NumericVector sigma(2);
  NumericVector loglike(2);

  for(i=0;i<=nbr;i++) {
    Z[i] = stpar[i];
    yearZ[i] = stpar[nbr+2+i];
  }
  Z[nbr+1] = stpar[nbr+1];

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ[i]>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ[i];
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(m=1;m<=count;m++) {
    for(i=0;i<=nbr;i++) {
      if(yearZ(i)>=m) sumy(i,m-1) = 0.;
      else sumy(i,m-1) = m - yearZ(i);
    }
    for(i=0;i<=nbr+1;i++) {
      if(i==0) {
        if(yearZ(i)>=m) int_lower(i,m-1) = Lc/Linf;
        else int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
        int_upper(i,m-1) = 1.;
      }
      if(i>0 & i<nbr+1) {
        if(yearZ(i-1)<m & yearZ(i)>=m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i)<m) {
          int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i-1)>=m) {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
      if(i==nbr+1) {
        if(yearZ(i-1)<m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        else {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
    }
  }

  for(m=0;m<count;m++) {
    denom[m] = 0.;
    numsum[m] = 0.;
    biomass[m] = 0.;

    for(i=0;i<=nbr+1;i++) {
      a(i,m) = 1.;
      r(i,m) = 1.;
      w(i,m) = 1.;

      if(i<nbr+1) s(i,m) = 1. - exp(-(Z[nbr+1-i]+K) * dy(nbr-i,m));
      if(i==nbr+1) s(i,m) = 1.;

      if(i>0) {
        for(j=0;j<=i-1;j++) {
          a(i,m) *= exp(-Z[nbr+1-j] * dy(nbr-j,m));
          r(i,m) *= exp(-(Z[nbr+1-j] + K) * dy(nbr-j,m));
          w(i,m) *= exp(Z[nbr+1-i] * dy(nbr-j,m));
        }
      }

      if(i<=nbr) denom[m] += a(i,m)*(1. - exp(-Z[nbr+1-i] * dy(nbr-i,m)))/Z[nbr+1-i];
      if(i==nbr+1) denom[m] += a(i,m)/Z[nbr+1-i];

      numsum[m] += r(i,m) * s(i,m) / (Z[nbr+1-i] + K);
      biomass[m] += a(i,m) * w(i,m) * pow(1 - Lc/Linf, -Z[nbr+1-i]/K) *
        (pbeta_inc(int_upper(nbr+1-i,m), b+1, Z[nbr+1-i]/K)-pbeta_inc(int_lower(nbr+1-i,m), b+1, Z[nbr+1-i]/K));
    }
    num[m] = Linf * (denom[m] - (1. - Lc/Linf) * numsum[m]);
    Lpred[m] = num[m]/denom[m];
  }
  for(m=0;m<count;m++) {
    if(ss[m]>0) nyear(0) += 1.;
    if(CPUE[m]>0) nyear(1) += 1.;
  }

  double sum_q = 0.;
  double sum_q2 = 0.;

  if(loglikeCPUE == 0) {
    for(m=0;m<count;m++)
    {
      if(CPUE[m]>0) sum_q += log(CPUE[m]/biomass[m]);
    }
    q = exp(sum_q/nyear(1));
  }
  if(loglikeCPUE == 1) {
    for(m=0;m<count;m++) {
      if(CPUE[m]>0) {
        sum_q += CPUE[m] * biomass[m];
        sum_q2 += pow(biomass[m],2);
      }
    }
    q = sum_q/sum_q2;
  }
  for(m=0;m<count;m++) {
    Ipred[m] = q * biomass[m];
    if(ss[m]>0) sum_square[0] += ss[m] * pow(Lbar[m]-Lpred[m],2.);
    if(CPUE[m]>0) {
      if(loglikeCPUE == 0) sum_square[1] += pow(log(CPUE[m]/Ipred[m]),2);
      if(loglikeCPUE == 1) sum_square[1] += pow(CPUE[m] - Ipred[m], 2);
    }
  }
  for(i=0;i<2;i++) {
    sigma[i] = pow(sum_square[i]/nyear[i], 0.5);
    loglike[i] = -nyear[i] * log(sigma[i]) - 0.5 * sum_square[i]/pow(sigma[i],2);
  }

  return List::create(Lpred, Ipred, q, sigma, loglike);
}


// [[Rcpp::export]]
double MLWPUEfullnegLL(NumericVector stpar, NumericVector Lbar, NumericVector ss, NumericVector CPUE,
                       NumericVector LH, double Lc, int nbreaks, int loglikeCPUE, int spCont) {

  int i;
  int j;
  int m;

  int count = Lbar.size();
  int nbr = nbreaks-1;

  double Linf = LH[0];
  double K = LH[1];
  double b = LH[2];

  NumericVector Z(nbreaks+1);
  NumericVector yearZ(nbreaks);
  NumericMatrix dy(nbreaks,count);
  NumericMatrix sumy(nbreaks,count);

  NumericMatrix a(nbreaks+1,count);
  NumericMatrix s(nbreaks+1,count);
  NumericMatrix r(nbreaks+1,count);
  NumericMatrix w(nbreaks+1,count);
  NumericMatrix int_upper(nbreaks+1,count);
  NumericMatrix int_lower(nbreaks+1,count);
  NumericVector denom(count);
  NumericVector numsum(count);
  NumericVector num(count);
  NumericVector biomass(count);
  NumericVector Lpred(count);
  NumericVector Ipred(count);

  NumericVector sigma(2);
  double q;
  double nLL = 0.;

  for(i=0;i<=nbr;i++) {
    Z[i] = stpar[i];
    yearZ[i] = stpar[nbr+2+i];
  }
  Z[nbr+1] = stpar[nbr+1];
  q = stpar[2*nbr+3];
  sigma[0] = stpar[2*nbr+4];
  sigma[1] = stpar[2*nbr+5];

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      if(yearZ[i]>=m) dy(i,m-1) = 0.;
      else dy(i,m-1) = m-yearZ[i];
    }
  }
  if(nbreaks>1) {
    for(i=0;i<nbr;i++){
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(m=1;m<=count;m++) {
    for(i=0;i<=nbr;i++) {
      if(yearZ(i)>=m) sumy(i,m-1) = 0.;
      else sumy(i,m-1) = m - yearZ(i);
    }
    for(i=0;i<=nbr+1;i++) {
      if(i==0) {
        if(yearZ(i)>=m) int_lower(i,m-1) = Lc/Linf;
        else int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
        int_upper(i,m-1) = 1.;
      }
      if(i>0 & i<nbr+1) {
        if(yearZ(i-1)<m & yearZ(i)>=m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i)<m) {
          int_lower(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1));
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        if(yearZ(i-1)>=m) {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
      if(i==nbr+1) {
        if(yearZ(i-1)<m) {
          int_lower(i,m-1) = Lc/Linf;
          int_upper(i,m-1) = 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1));
        }
        else {
          int_lower(i,m-1) = 0.;
          int_upper(i,m-1) = 0.;
        }
      }
    }
  }

  for(m=0;m<count;m++) {
    denom[m] = 0.;
    numsum[m] = 0.;
    biomass[m] = 0.;

    for(i=0;i<=nbr+1;i++) {
      a(i,m) = 1.;
      r(i,m) = 1.;
      w(i,m) = 1.;

      if(i<nbr+1) s(i,m) = 1. - exp(-(Z[nbr+1-i]+K) * dy(nbr-i,m));
      if(i==nbr+1) s(i,m) = 1.;

      if(i>0) {
        for(j=0;j<=i-1;j++) {
          a(i,m) *= exp(-Z[nbr+1-j] * dy(nbr-j,m));
          r(i,m) *= exp(-(Z[nbr+1-j] + K) * dy(nbr-j,m));
          w(i,m) *= exp(Z[nbr+1-i] * dy(nbr-j,m));
        }
      }

      if(i<=nbr) denom[m] += a(i,m)*(1. - exp(-Z[nbr+1-i] * dy(nbr-i,m)))/Z[nbr+1-i];
      if(i==nbr+1) denom[m] += a(i,m)/Z[nbr+1-i];

      numsum[m] += r(i,m) * s(i,m) / (Z[nbr+1-i] + K);
      biomass[m] += a(i,m) * w(i,m) * pow(1 - Lc/Linf, -Z[nbr+1-i]/K) *
        (pbeta_inc(int_upper(nbr+1-i,m), b+1, Z[nbr+1-i]/K)-pbeta_inc(int_lower(nbr+1-i,m), b+1, Z[nbr+1-i]/K));
    }
    num[m] = Linf * (denom[m] - (1. - Lc/Linf) * numsum[m]);
    Lpred[m] = num[m]/denom[m];
    Ipred[m] = q * biomass[m];

    if(ss[m]>0) nLL += -log(sigma[0]) - 0.5 * ss[m] * pow(Lbar[m]-Lpred[m],2)/pow(sigma[0],2);
    if(CPUE[m]>0) {
      if(loglikeCPUE == 0) nLL += -log(sigma[1]) - 0.5 * pow(log(CPUE[m]/Ipred[m]),2)/pow(sigma[1],2);
      if(loglikeCPUE == 1) nLL += -log(sigma[1]) - 0.5 * pow(CPUE[m] - Ipred[m],2)/pow(sigma[1],2);
    }
  }
  nLL *= -1;
  return nLL;
}
