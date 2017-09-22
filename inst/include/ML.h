
//template<class Type>
//Type objective_function<Type>::operator() ()
//{
  DATA_SCALAR(Linf);
  DATA_SCALAR(K);
  DATA_SCALAR(Lc);
  DATA_INTEGER(nbreaks);
  DATA_VECTOR(Lbar);
  DATA_VECTOR(ss); // Vector of sample sizes for mean lengths
  
  PARAMETER_VECTOR(Z);
  PARAMETER_VECTOR(yearZ);
  
  int count = Lbar.size();
  
  // Calculate a, v, r, s, denom, num, numsum, Lpred
  modelOutput<Type> output;
  if(nbreaks == 0) output = model_output_eq(Z(0), Linf, K, Lc, Type(1e-3), count);
  if(nbreaks > 0) output = model_output(Z, yearZ, Linf, K, Lc, Type(1e-3), nbreaks, count);
  
  vector<Type> Lpred(count);
  Lpred = output.Lpred;
  REPORT(Lpred);
  
  // Analytical solution for sigmaL
  Type sigma = estimate_sigmaL(Lbar, Lpred, ss, count);  
  ADREPORT(sigma);
  
  // Negative log-likelihood of mean lengths
  Type nll = nll_Lbar(Lbar, Lpred, ss, sigma, count);

  REPORT(Z);
  REPORT(yearZ);
  REPORT(nll);
  return nll;
//}
