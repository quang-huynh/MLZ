//template<class Type>
//Type objective_function<Type>::operator() ()
//{
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
    vector<Type> Ztmp(nbreaks+1);
	Ztmp = Z.col(sp);
	vector<Type> yrZtmp(nbreaks);
	yrZtmp = yearZ.col(sp);
	Type Linftemp = Linf(sp);
	Type Ktemp = K(sp);
	Type Lctemp = Lc(sp);
	Type btemp = 1e-4;	
	
	if(nbreaks == 0) output = model_output_eq(Ztmp(0), Linftemp, Ktemp, Lctemp, btemp, count);
    if(nbreaks > 0) output = model_output(Ztmp, yrZtmp, Linftemp, Ktemp, Lctemp, btemp, nbreaks, count);
    Lpred.col(sp) = output.Lpred;
    
    vector<Type> Lbartmp;
    Lbartmp = Lbar.col(sp);
    vector<Type> sstmp;
    sstmp = ss.col(sp);
    
    // Analytical solution for sigmaL
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


