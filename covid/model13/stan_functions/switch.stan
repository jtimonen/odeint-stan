// Sigmoidal switch
real switch_eta(real t, real t1, real eta, real nu, real xi) {
  return(eta+(1-eta)/(1+exp(xi*(t-t1-nu))));
}