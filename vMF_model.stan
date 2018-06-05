/**
 * Fit data generated from the vMF
 * distribution to infer mu and kappa.
 */

functions {

  /**
   * Define the vMF PDF.
   * NB: cannot be vectorised.
   */
  real vMF_lpdf(vector v, vector mu, real kappa) {

    real lprob = kappa * dot_product(v, mu) + log(kappa) - log(4 * pi() * sinh(kappa));
    return lprob;
    
  }
  
}

data {

  int N;
  unit_vector[3] data[N];
  
}

parameters {

  unit_vector[3] mu;
  real kappa;
  
}

model {

  for (i in 1:N) {
    data[i] ~ vMF(mu, kappa);
  }

}
