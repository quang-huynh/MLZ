#include <TMB.hpp>
#include <tiny_ad/beta/pbeta.hpp>
#include <cppad/cppad.hpp>

template<class Type>
Type square(Type x){return x*x;}

template<class Type>
Type pbeta_inc(Type x, Type alpha, Type beta) {
  Type answer = exp(lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta));
  answer *= pbeta(x, alpha, beta);
  return answer;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(Linf);
  DATA_SCALAR(K);
  DATA_SCALAR(Lc);
  DATA_SCALAR(b);
  DATA_INTEGER(nbreaks);
  DATA_VECTOR(Lbar);
  DATA_VECTOR(ss);
  DATA_VECTOR(CPUE);
  DATA_INTEGER(loglikeCPUE);
  DATA_INTEGER(spCont);

  PARAMETER_VECTOR(Z);
  PARAMETER_VECTOR(yearZ);

  int i;
  int j;
  int m;

  int count = Lbar.size();
  int nbr = nbreaks - 1;

  vector<Type> loglike(2);
  vector<Type> sigma(2);
  vector<Type> nyrs(2);

  vector<Type> sum_square(2);

  Type q;
  Type nll;
  
  Type sum_q = 0.;
  Type sum_q2 = 0.;
  
  matrix<Type> dy(nbreaks,count);
  matrix<Type> a(nbreaks+1,count);
  matrix<Type> s(nbreaks+1,count);
  matrix<Type> r(nbreaks+1,count);
  matrix<Type> w(nbreaks+1,count);
  matrix<Type> sumy(nbreaks,count);
  matrix<Type> int_lower(nbreaks+1,count);
  matrix<Type> int_upper(nbreaks+1,count);

  vector<Type> denom(count);
  vector<Type> numsum(count);
  vector<Type> num(count);
  vector<Type> biomass(count);

  vector<Type> Lpred(count);
  vector<Type> Ipred(count);

  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
      Type mm = m;
      dy(i,m-1) = CppAD::CondExpGe(yearZ(i), mm, Type(0), mm-yearZ(i));
	  }
  }

  if(nbreaks>1) {
    for(i=0;i<nbr;i++) {
      for(m=0;m<count;m++) dy(i,m) -= dy(i+1,m);
    }
  }

  for(m=1;m<=count;m++) {
    Type mm = m;
    for(i=0;i<=nbr;i++) {
      sumy(i,m-1) = CppAD::CondExpGe(yearZ(i), mm, Type(0), mm-yearZ(i));
    }
    for(i=0;i<=nbr+1;i++) {
      if(i==0) {
        int_lower(i,m-1) = CppAD::CondExpGe(yearZ(i), mm, Lc/Linf, 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1)));
        int_upper(i,m-1) = 1.;
      }
      if(i>0 & i<nbr+1) {
        //Type lower_test = CppAD::CondExpLe(mm, yearZ(i-1), Type(0), Type(1));
        //if(lower_test == 1) lower_test = CppAD::CondExpLe(mm, yearZ(i), Lc/Linf, 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1)));
        //int_lower(i,m-1) = lower_test;
        int_lower(i,m-1) = CppAD::CondExpLe(mm, yearZ(i-1), Type(0), CppAD::CondExpLe(mm, yearZ(i), Lc/Linf, 
                                            1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1))));
        int_upper(i,m-1) = CppAD::CondExpLe(mm, yearZ(i-1), Type(0), 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1)));
      }
      if(i==nbr+1) {
        int_lower(i,m-1) = CppAD::CondExpLt(yearZ(i-1), mm, Lc/Linf, Type(0));
        int_upper(i,m-1) = CppAD::CondExpLt(yearZ(i-1), mm, 1 - (1 - Lc/Linf) * exp(-K * sumy(i-1,m-1)), Type(0));
      }
    }
  }

  for(m=0;m<count;m++) {
    denom(m) = 0.;
    numsum(m) = 0.;
    biomass(m) = 0.;

    for(i=0;i<=nbr+1;i++) {
      a(i,m) = 1.;
      r(i,m) = 1.;
      w(i,m) = 1.;
      if(i<nbr+1) s(i,m) = 1. - exp(-(Z(nbr+1-i)+K) * dy(nbr-i,m));
      if(i==nbr+1) s(i,m) = 1.;

	    if(i>0) {
	      for(j=0;j<=i-1;j++) {
          a(i,m) *= exp(-Z(nbr+1-j) * dy(nbr-j,m));
          r(i,m) *= exp(-(Z(nbr+1-j) + K) * dy(nbr-j,m));
          w(i,m) *= exp(Z(nbr+1-i) * dy(nbr-j,m));
        }
      }

      if(i<=nbr) denom(m) += a(i,m) * (1. - exp(-Z(nbr+1-i) * dy(nbr-i,m)))/Z(nbr+1-i);
      if(i==nbr+1) denom(m) += a(i,m)/Z(nbr+1-i);

      numsum(m) += r(i,m) * s(i,m) / (Z(nbr+1-i) + K);
      biomass(m) += a(i,m) * w(i,m) * pow(1 - Lc/Linf, -Z(nbr+1-i)/K) *
        (pbeta_inc(int_upper(nbr+1-i,m), b+1, Z(nbr+1-i)/K) - pbeta_inc(int_lower(nbr+1-i,m), b+1, Z(nbr+1-i)/K));
    }
    num(m) = Linf * (denom(m) - (1. - Lc/Linf)*numsum(m));
    Lpred(m) = num(m)/denom(m);
    
    if(ss(m)>0) {
      nyrs(0) += 1.;
      sum_square(0) += ss(m) * square(Lbar(m)-Lpred(m));
    }
    if(CPUE(m)>0) {
      nyrs(1) += 1.;
      if(loglikeCPUE == 0) sum_q += log(CPUE(m)/biomass(m));
      if(loglikeCPUE == 1) {
        sum_q += CPUE(m) * biomass(m);
        sum_q2 += square(biomass(m));
      }
    }
  }

  REPORT(Lpred);

  if(loglikeCPUE == 0) q = exp(sum_q/nyrs(1));
  if(loglikeCPUE == 1) q = sum_q/sum_q2;

  ADREPORT(q);

  for(m=0;m<count;m++) {
    Ipred(m) = q * biomass(m);
	if(CPUE(m)>0) {
	  if(loglikeCPUE == 0) sum_square(1) += square(log(CPUE(m)/Ipred(m)));
	  if(loglikeCPUE == 1) sum_square(1) += square(CPUE(m) - Ipred(m));
	}
  }

  REPORT(Ipred);

  for(i=0;i<2;i++) {
    sigma(i) = sqrt(sum_square(i)/nyrs(i));
    loglike(i) = -nyrs(i) * log(sigma(i)) - 0.5 * sum_square(i)/square(sigma(i));
  }

  ADREPORT(sigma);

  nll = loglike.sum();
  nll *= -1;

  REPORT(loglike);
  REPORT(nll);
  REPORT(Z);
  REPORT(yearZ);
  REPORT(q);
  REPORT(int_lower);
  REPORT(int_upper);

  return nll;
}


