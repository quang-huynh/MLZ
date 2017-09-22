//template<class Type>
//Type objective_function<Type>::operator() ()
//{
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
    vector<Type> Ztmp(nbreaks+1);
    Ztmp = Z.col(sp);
    Type Linftemp = Linf(sp);
	  Type Ktemp = K(sp);
	  Type Lctemp = Lc(sp);
	  Type btemp = 1e-4;	
	
	  output = model_output(Ztmp, yearZ, Linftemp, Ktemp, Lctemp, btemp, nbreaks, count);
	  Lpred.col(sp) = output.Lpred;
    
    // Analytical solution for sigmaL
    vector<Type> Lbartmp;
    Lbartmp = Lbar.col(sp);
    vector<Type> sstmp;
    sstmp = ss.col(sp);
    
    Type sigma = estimate_sigmaL(Lbartmp, output.Lpred, sstmp, count);
    sigmaL(sp) = sigma;
    
    // Negative log-likelihood of mean lengths
    nllc(sp) = nll_Lbar(Lbartmp, output.Lpred, sstmp, sigma, count);
  }
  REPORT(Lpred);
  ADREPORT(sigmaL);  
  
  Type nll = nllc.sum();
  
  REPORT(Z);
  REPORT(yearZ);
  REPORT(nllc);
  REPORT(nll);
  return nll;
//}


