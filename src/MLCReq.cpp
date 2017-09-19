#include <TMB.hpp>
#include <tiny_ad/beta/pbeta.hpp>

template <class Type> 
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
  DATA_VECTOR(Lbar);
  DATA_VECTOR(ss);
  DATA_VECTOR(CPUE);
  DATA_INTEGER(loglikeCPUE);
  DATA_INTEGER(isNPUE); // WPUE or NPUE
  
  PARAMETER(Z);
  
  int i;
  int m;
  
  vector<Type> loglike(2);
  vector<Type> sigma(2);
  vector<Type> nyrs(2);
  
  vector<Type> sum_square(2);
  sum_square.setZero();
  
  Type q;
  Type nll;
  
  Type sum_q = 0.;
  Type sum_q2 = 0.;
  
  vector<Type> Lpred(Lbar.size());
  vector<Type> Ipred(Lbar.size());
  Type denom = 0.;
  
  for(m=0;m<Lbar.size();m++) {
    Lpred(m) = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
    if(isNPUE == 1) denom = 1/Z;
    if(isNPUE == 0) denom = pow(1 - Lc/Linf, -Z/K) * 
      pbeta_inc(Type(1), b+1, Z/K) - pbeta_inc(Lc/Linf, b+1, Z/K);
    
    if(ss(m)>0) {
      nyrs(0) += 1.;
      sum_square(0) += ss(m) * square(Lbar(m)-Lpred(m));
    }
    if(CPUE(m)>0) {
      nyrs(1) += 1.;
      if(loglikeCPUE == 0) sum_q += log(CPUE(m)/denom);
      if(loglikeCPUE == 1) {
        sum_q += CPUE(m) * denom;
        sum_q2 += square(denom);
      }
    }
  }  
  REPORT(Lpred);
  
  if(loglikeCPUE == 0) q = exp(sum_q/nyrs(1));
  if(loglikeCPUE == 1) q = sum_q/sum_q2;
  
  ADREPORT(q);
  
  for(m=0;m<Lbar.size();m++) {
    Ipred(m) = q * denom;
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
  REPORT(q);
  
  return nll;
}  
