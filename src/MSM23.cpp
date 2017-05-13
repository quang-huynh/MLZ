#include <TMB.hpp>
#include <cppad/cppad.hpp>

template<class Type> 
  Type square(Type x){return x*x;}

template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Linf);
  DATA_VECTOR(K);
  DATA_VECTOR(Lc);
  DATA_VECTOR(M);
  DATA_INTEGER(nbreaks);
  DATA_INTEGER(nspec);
  DATA_MATRIX(Lbar);
  DATA_MATRIX(ss);
  DATA_INTEGER(spCont);
  
  PARAMETER_VECTOR(Z1);
  PARAMETER_VECTOR(delta);
  PARAMETER_VECTOR(epsilon);
  PARAMETER_VECTOR(yearZ);
  
  int i;
  int j;
  int m;
  int w;
  
  int count = Lbar.col(0).size();
  int nbr = nbreaks - 1;
    
  vector<Type> loglikesp(nspec);
  vector<Type> sigma(nspec);
  vector<Type> nyrs(nspec);
  
  Type nll;
  
  vector<Type> sum_square_Lpred(nspec);
  
  matrix<Type> Z(nbreaks+1,nspec);
  
  matrix<Type> dy(nbreaks,count);
  matrix<Type> a(nbreaks+1,count);
  matrix<Type> s(nbreaks+1,count);
  matrix<Type> r(nbreaks+1,count);
  
  matrix<Type> denom(count,nspec);
  matrix<Type> numsum(count,nspec);
  matrix<Type> num(count,nspec);
  
  matrix<Type> Lpred(count,nspec);
  
  for(w=0;w<nspec;w++) {
    Z(0,w) = Z1(w);
    for(i=1;i<=nbr+1;i++) {
      if(w==0) Z(i,w) = delta(i-1) * Z(i-1,w) + (1-delta(i-1)) * M(w);
      if(w>0) Z(i,w) = delta(i-1) * epsilon(w-1) * Z(i-1,w) + (1-delta(i-1)*epsilon(w-1)) * M(w);
    }
  }
  ADREPORT(Z);
  
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
  
  for(w=0;w<nspec;w++) {
    sum_square_Lpred(w) = 0.;
    nyrs(w) = 0.;
    
    for(m=0;m<count;m++) {
      denom(m,w) = 0.;
      numsum(m,w) = 0.;
      
      for(i=0;i<=nbr+1;i++) {
        a(i,m) = 1.;
        r(i,m) = 1.;
        if(i<nbr+1) s(i,m) = 1. - exp(-(Z(nbr+1-i,w)+K(w)) * dy(nbr-i,m));
        if(i==nbr+1) s(i,m) = 1.;
        
        if(i>0) {
          for(j=0;j<=i-1;j++) {
            a(i,m) *= exp(-Z(nbr+1-j,w) * dy(nbr-j,m));
            r(i,m) *= exp(-(Z(nbr+1-j,w) + K(w)) * dy(nbr-j,m));
          }
        }
        
        if(spCont == 1) {
          if(i<=nbr) denom(m,w) += a(i,m)*(1. - exp(-Z(nbr+1-i,w) * dy(nbr-i,m)))/Z(nbr+1-i,w);
          if(i==nbr+1) denom(m,w) += a(i,m)/Z(nbr+1-i,w);
          numsum(m,w) += r(i,m) * s(i,m) / (Z(nbr+1-i,w) + K(w));
        }
        if(spCont == 0) {
          if(i<=nbr) denom(m,w) += a(i,m)*(1. - exp(-Z(nbr+1-i,w) * dy(nbr-i,m)))/(1 - exp(-Z(nbr+1-i,w)));
          if(i==nbr+1) denom(m,w) += a(i,m)/(1 - exp(-Z(nbr+1-i,w)));
          numsum(m,w) += r(i,m) * s(i,m) / (1 - exp(-(Z(nbr+1-i,w) + K(w))));
        }
      }
      num(m,w) = Linf(w) * (denom(m,w) - (1. - Lc(w)/Linf(w))*numsum(m,w));
      Lpred(m,w) = num(m,w)/denom(m,w);
      
      if(ss(m,w)>0) {
        sum_square_Lpred(w) += ss(m,w) * square(Lbar(m,w)-Lpred(m,w));
        nyrs(w) += 1.;
      }
    }
    sigma(w) = sqrt(sum_square_Lpred(w)/nyrs(w));
    loglikesp(w) = -nyrs(w) * log(sigma(w)) - 0.5 * sum_square_Lpred(w)/square(sigma(w));
  }
  
  REPORT(Lpred);
  ADREPORT(sigma);
  
  nll = loglikesp.sum();
  nll *= -1;
  REPORT(loglikesp);
  REPORT(nll);
  REPORT(Z);
  REPORT(yearZ);
  return nll;
}
