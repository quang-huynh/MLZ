

template<class Type> 
Type square(Type x) {return x*x;}

template<class Type>
Type pbeta_inc(Type x, Type alpha, Type beta) {
  Type answer = exp(lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta));
  answer *= pbeta(x, alpha, beta);
  return answer;
}

template<class Type>
struct modelOutput{
  matrix<Type> dy;
  matrix<Type> sumy;
  matrix<Type> a;
  matrix<Type> v;
  matrix<Type> s;
  matrix<Type> r;
  vector<Type> denom;
  vector<Type> numsum;
  vector<Type> num;
  vector<Type> biomass;
  vector<Type> Lpred;  
};

template<class Type>
modelOutput<Type> model_output(vector<Type> Z, vector<Type> yearZ, Type Linf, Type K, Type Lc, Type b, int nbreaks, int count) {
  int i, j, m;
  // Adjust for index start at 0 for matrices and vectors in Cpp
  int nbr = nbreaks - 1;
  matrix<Type> dy(nbreaks,count);
  matrix<Type> sumy(nbreaks,count);
  matrix<Type> int_lower(nbreaks+1,count);
  matrix<Type> int_upper(nbreaks+1,count);
  
  matrix<Type> a(nbreaks+1,count);
  matrix<Type> v(nbreaks+1,count);
  matrix<Type> s(nbreaks+1,count);
  matrix<Type> r(nbreaks+1,count);
  matrix<Type> w(nbreaks+1,count);
  
  vector<Type> denom(count);
  vector<Type> numsum(count);
  vector<Type> num(count);
  vector<Type> biomass(count);
  
  denom.setZero();
  numsum.setZero();
  biomass.setZero();
  
  vector<Type> Lpred(count);
  Type mm;  
  for(i=0;i<=nbr;i++) {
    for(m=1;m<=count;m++) {
	    mm = m;
      dy(i,m-1) = CppAD::CondExpGe(yearZ(i), mm, Type(0.), mm-yearZ(i));
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
      sumy(i,m-1) = CppAD::CondExpGe(yearZ(i), mm, Type(0.), mm-yearZ(i));
    }
    for(i=0;i<=nbr+1;i++) {
      if(i==0) {
        int_lower(i,m-1) = CppAD::CondExpGe(yearZ(i), mm, Lc/Linf, 1 - (1 - Lc/Linf) * exp(-K * sumy(i,m-1)));
        int_upper(i,m-1) = 1.;
      }
      if((i>0) & (i<nbr+1)) {
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
    for(i=0;i<=nbr+1;i++) {
      a(i,m) = 1.;
      r(i,m) = 1.;
      w(i,m) = 1.;
	  
      if(i>0) {
        for(j=0;j<=i-1;j++) {
          a(i,m) *= exp(-Z(nbr+1-j) * dy(nbr-j,m));
          r(i,m) *= exp(-(Z(nbr+1-j) + K) * dy(nbr-j,m));
          w(i,m) *= exp(Z(nbr+1-i) * dy(nbr-j,m));
        }
      }
      if(i<=nbr) {
        v(i,m) = 1. - exp(-Z(nbr+1-i) * dy(nbr-i,m));
        s(i,m) = 1. - exp(-(Z(nbr+1-i)+K) * dy(nbr-i,m));
      }
      if(i==nbr+1) {
        v(i,m) = 1.;
        s(i,m) = 1.;
      }
      denom(m) += a(i,m) * v(i,m)/Z(nbr+1-i);
      numsum(m) += r(i,m) * s(i,m) /(Z(nbr+1-i) + K);
      biomass(m) += a(i,m) * w(i,m) * pow(1 - Lc/Linf, -Z(nbr+1-i)/K) *
        (pbeta_inc(int_upper(nbr+1-i,m), b+1, Z(nbr+1-i)/K) - pbeta_inc(int_lower(nbr+1-i,m), b+1, Z(nbr+1-i)/K));
    }
    num(m) = Linf * (denom(m) - (1. - Lc/Linf) * numsum(m));
    Lpred(m) = num(m)/denom(m);
  }
  
  modelOutput<Type> output;
  output.dy = dy;
  output.a = a;
  output.v = v;
  output.r = r;
  output.s = s;
  output.numsum = numsum;
  output.num = num;
  output.denom = denom;
  output.biomass = biomass;
  output.Lpred = Lpred;
  
  return output;
}


template<class Type>
modelOutput<Type> model_output_eq(Type Z, Type Linf, Type K, Type Lc, Type b, int count) {
  int m;
  vector<Type> denom(count);
  vector<Type> biomass(count);
  vector<Type> Lpred(count);
  Type Zeq = Z;
  for(m=0;m<count;m++) {
    denom(m) = 1/Zeq;
    biomass(m) = pow(1 - Lc/Linf, -Zeq/K) * (pbeta_inc(Type(1), b+1, Zeq/K) - pbeta_inc(Lc/Linf, b+1, Zeq/K));
    Lpred(m) = Linf * (1 - (Zeq/(Zeq+K)) * (1 - Lc/Linf));
  }
  
  modelOutput<Type> output;
  output.denom = denom;
  output.biomass = biomass;
  output.Lpred = Lpred;
  
  return output;
}

template<class Type>
Type estimate_sigmaL(vector<Type> Lbar, vector<Type> Lpred, vector<Type> ss, int count) {
  Type sum_square = 0.;
  Type nyrs = 0.;
  for(int m=0;m<count;m++) {
    if(ss(m)>0) {
      sum_square += ss(m) * square(Lbar(m)-Lpred(m));
      nyrs += 1.;
    }
  }  
  Type sigmaL = sqrt(sum_square/nyrs);
  return sigmaL;
}

//template<class Type>
//vector<Type> estimate_sigmaL(matrix<Type> Lbar, matrix<Type> Lpred, matrix<Type> ss, int nspec, int count) {
//  vector<Type> sum_square(nspec);
// vector<Type> nyrs(nspec);
//  vector<Type> sigmaL(nspec);
//  sum_square.setZero();
//  nyrs.setZero();
//  for(int sp=0;sp<nspec;sp++) {
//    for(int m=0;m<count;m++) {
//      if(ss(m,sp)>0) {
//        sum_square(sp) += ss(m,sp) * square(Lbar(m,sp)-Lpred(m,sp));
//        nyrs(sp) += 1.;
//      }
//    }
//    sigmaL(sp) = sqrt(sum_square(sp)/nyrs(sp));
// }
//  return sigmaL;
//}

template<class Type>
Type estimate_q(vector<Type> CPUE, vector<Type> population, int loglikeCPUE, int count) {
  Type sum_q = 0.;
  Type sum_q2 = 0.;
  Type nyrs = 0.;
  Type q = 0.;
  for(int m=0;m<count;m++) {
    if(CPUE(m)>0) {
      if(loglikeCPUE == 0) { 
	      nyrs += 1.;	
	      sum_q += log(CPUE(m)/population(m));
      } else {
        sum_q += CPUE(m) * population(m);
        sum_q2 += square(population(m));
      }
    }
  }
  if(loglikeCPUE == 0) {
    q = exp(sum_q/nyrs);
  } else {
    q = sum_q/sum_q2;
  }
  return q;
}

template<class Type>
vector<Type> calculate_Ipred(vector<Type> population, Type q, int count) {
  vector<Type> Ipred(count);
  for(int m=0;m<count;m++) Ipred(m) = q * population(m);
  return Ipred;
}

template<class Type>
Type estimate_sigmaI(vector<Type> CPUE, vector<Type> Ipred, int count, int loglikeCPUE) {
  Type sum_square = 0.;
  Type nyrs = 0.;
  for(int m=0;m<count;m++) {
    if(CPUE(m)>0) {
      nyrs += 1.;
      if(loglikeCPUE == 0) {
        sum_square += square(log(CPUE(m)/Ipred(m)));
      } else {
        if(loglikeCPUE == 1) sum_square += square(CPUE(m) - Ipred(m));
      }
    }
  }
  Type sigmaI = sqrt(sum_square/nyrs);
  return sigmaI;
}

template<class Type>
Type nll_Lbar(vector<Type> Lbar, vector<Type> Lpred, vector<Type> ss, Type sigma, int count) {
  Type nll = 0.;
  for(int m=0;m<count;m++) {
    if(ss(m)>0) nll -= dnorm(Lbar(m), Lpred(m), sigma/sqrt(ss(m)), true);
  }
  return nll;
}

//template<class Type>
//vector<Type> nll_Lbar(matrix<double> Lbar, matrix<Type> Lpred, matrix<Type> ss, vector<Type> sigma, int nspec, int count) {
//  vector<Type> nll(nspec);
//  nll.setZero();
//  for(int w=0;w<nspec;w++) {
//    for(int m=0;m<count;m++) {
//      if(ss(m,w)>0) nll(w) -= dnorm(Lbar(m,w), Lpred(m,w), sigma(w)/sqrt(ss(m,w)), true);
//    }
//  }
//  return nll;
//}

template<class Type>
Type nll_CPUE(vector<Type> CPUE, vector<Type> Ipred, Type sigma, int count, int loglikeCPUE) {
  Type nll = 0.;
  for(int m=0;m<count;m++) {
    if(CPUE(m)>0) {
      if(loglikeCPUE == 0) {
        nll -= dnorm(log(CPUE(m)), log(Ipred(m)), sigma, true);
      } else {
        nll -= dnorm(CPUE(m), Ipred(m), sigma, true);
      }
    }   
  }
  return nll;
}


template<class Type>
matrix<Type> ML_effort_Neq(int astep, int n_season, Type Z_init, Type Z, Type seasD) {
  matrix<Type> N(astep,n_season);
  int a,k;
  N(0,0) = 1.;
  for(a=1;a<astep;a++) N(a,0) = N(a-1,0) * exp(-Z_init/seasD);
  if(n_season>1) {
    for(k=1;k<n_season;k++) {
      N(0,k) = 1.;
      for(a=1;a<astep;a++) N(a,k) = N(a-1,k-1) * exp(-Z/seasD);
    }
  }
  return N;
}

template<class Type>
matrix<Type> ML_effort_N(int astep, int n_season, Type Zprev, Type Z, Type seasD, matrix<Type> Nprev) {
  matrix<Type> N(astep,n_season);
  int a,k;
  N(0,0) = 1.;
  for(a=1;a<astep;a++) N(a,0) = Nprev(a-1,n_season-1) * exp(-Zprev/seasD);
  if(n_season>1) {
    for(k=1;k<n_season;k++) {
      N(0,k) = 1.;
      for(a=1;a<astep;a++) N(a,k) = N(a-1,k-1) * exp(-Z/seasD);
    }
  }
  return N;
}



#include "ML.hpp"
#include "MLCR.hpp"
#include "MSM1S.hpp" 
#include "MSM23.hpp"
#include "MLeffort.hpp" 

