
#ifndef MSM1S_hpp
#define MSM1S_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type MSM1S(objective_function<Type> *obj) {
  
  DATA_VECTOR(Linf);
  DATA_VECTOR(K);
  DATA_VECTOR(Lc);
  DATA_INTEGER(nbreaks);
  DATA_INTEGER(nspec);
  DATA_MATRIX(Lbar);
  DATA_MATRIX(ss); // Matrix of sample sizes for mean lengths
  
  PARAMETER_MATRIX(Z);
  PARAMETER_MATRIX(yearZ);
  
  int count = Lbar.col(0).size();
  int sp;
  
  matrix<Type> Lpred(count,nspec);
  vector<Type> sigmaL(nspec);
  vector<Type> nllc(nspec);
    
  for(sp=0;sp<nspec;sp++) {
    // Calculate a, v, r, s, denom, num, numsum, Lpred
    modelOutput<Type> output;
    vector<Type> Ztmp = Z.col(sp);
    vector<Type> yrZtmp = yearZ.col(sp);
    
    if(nbreaks == 0) {
      output = model_output_eq(Ztmp(0), Linf(sp), K(sp), Lc(sp), Type(1e-4), count);
    } else {
      output = model_output(Ztmp, yrZtmp, Linf(sp), K(sp), Lc(sp), Type(1e-4), nbreaks, count);
    }  
    
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


