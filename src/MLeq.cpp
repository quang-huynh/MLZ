#include <TMB.hpp>

template <class Type> 
Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(Linf);
  DATA_SCALAR(K);
  DATA_SCALAR(Lc);
  DATA_VECTOR(Lbar);
  DATA_VECTOR(ss);
  DATA_INTEGER(spCont)
  
  PARAMETER(Z);
  
  int m;
  
  Type nll;
  Type sigma;
  Type nyrs = 0.;
  
  Type sum_square_Lpred = 0.;
  
  vector<Type> Lpred(Lbar.size());
  
  for(m=0;m<Lbar.size();m++) {
    if(spCont == 1) Lpred(m) = Linf * (1 - (Z/(Z+K)) * (1 - Lc/Linf));
    if(spCont == 0) Lpred(m) = Linf * (1 - ((1 - exp(-Z))/(1 - exp(-(Z+K)))) * (1 - Lc/Linf));
    if(ss(m)>0) {
      sum_square_Lpred += ss(m) * square(Lbar(m)-Lpred(m));
      nyrs += 1.;
    }
  }  
  REPORT(Lpred);
  
  sigma = sqrt(sum_square_Lpred/nyrs);
  ADREPORT(sigma);
  
  nll = -nyrs * log(sigma) - 0.5 * sum_square_Lpred/square(sigma);
  nll *= -1;
  REPORT(nll);
  REPORT(Z);
  return nll;
}  
