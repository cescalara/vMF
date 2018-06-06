/**
 * Fit data generated from the vMF
 * distribution to infer mu and kappa.
 */

functions {

  /**
   * Define the vMF PDF.
   * NB: Cannot be vectorised.
   * Uses sinh(kappa) ~ exp(kappa)/2 
   * approximation for kappa > 100.
   */
  real vMF_lpdf(vector v, vector mu, real kappa) {

    real lprob;
    if (kappa > 100) {
      lprob = kappa * dot_product(v, mu) + log(kappa) - log(4 * pi()) - kappa + log(2);
    }
    else {
      lprob = kappa * dot_product(v, mu) + log(kappa) - log(4 * pi() * sinh(kappa));
    }
    return lprob;
    
  }
    
}

data {

  int<lower=0> N;
  unit_vector[3] d[N];
  
}

parameters {

  unit_vector[3] mu;
  real<lower=0> kappa;
  
}

model {

  for (i in 1:N) {
    d[i] ~ vMF(mu, kappa);
  }

}
