
#ifndef MSM23_hpp
#define MSM23_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type MSM23(objective_function<Type> *obj) {
  
  DATA_VECTOR(Linf);
  DATA_VECTOR(K);
  DATA_VECTOR(Lc);
  DATA_VECTOR(M);
  DATA_INTEGER(nbreaks);
  DATA_INTEGER(nspec);
  DATA_MATRIX(Lbar);
  DATA_MATRIX(ss); // Matrix of sample sizes for mean lengths
  
  PARAMETER_VECTOR(Z1);
  PARAMETER_VECTOR(delta);
  PARAMETER_VECTOR(epsilon);
  PARAMETER_VECTOR(yearZ);
  
  int count = Lbar.col(0).size();
  int i, sp;
  
  // Derive time series of Z
  matrix<Type> Z(nbreaks+1,nspec);
  int nbr = nbreaks - 1;
  for(sp=0;sp<nspec;sp++) {
    Z(0,sp) = Z1(sp);
    for(i=1;i<=nbr+1;i++) {
      if(sp==0) Z(i,sp) = delta(i-1) * Z(i-1,sp) + (1-delta(i-1)) * M(sp);
      if(sp>0) Z(i,sp) = delta(i-1) * epsilon(sp-1) * Z(i-1,sp) + (1-delta(i-1)*epsilon(sp-1)) * M(sp);
    }
  }
  ADREPORT(Z);
  
  matrix<Type> Lpred(count,nspec);
  vector<Type> sigmaL(nspec);
  vector<Type> nllc(nspec);
  
  for(sp=0;sp<nspec;sp++) {
    // Calculate a, v, r, s, denom, num, numsum, Lpred
    modelOutput<Type> output;
    vector<Type> Ztmp = Z.col(sp);
	
	  output = model_output(Ztmp, yearZ, Linf(sp), K(sp), Lc(sp), Type(1e-4), nbreaks, count);
	  Lpred.col(sp) = output.Lpred;
    
    // Analytical solution for sigmaL
    vector<Type> Lbartmp = Lbar.col(sp);
    vector<Type> sstmp = ss.col(sp);
    
    sigmaL(sp) = estimate_sigmaL(Lbartmp, output.Lpred, sstmp, count);
    
    // Negative log-likelihood of mean lengths
    nllc(sp) = nll_Lbar(Lbartmp, output.Lpred, sstmp, sigmaL(sp), count);
  }
  REPORT(Lpred);
  ADREPORT(sigmaL);
  
  Type nll = nllc.sum();
  
  REPORT(Z);
  REPORT(yearZ);
  REPORT(nllc);
  REPORT(nll);
  return nll;
  
}
  
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
  
#endif

