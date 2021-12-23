data {
  int<lower=1> nreps;
  int<lower=1> notus;
  int datamatrix[nreps, notus];
}

parameters {
  real<lower=0> theta;
  simplex[notus] pi;
  simplex[notus] p[nreps];
}

model {
    target += exponential_lpdf(theta | 0.001);
    target += dirichlet_lpdf(pi | rep_vector(0.000000001, notus));
    for(j in 1:nreps){
      target += dirichlet_lpdf(p[j] | theta*pi);
      target += multinomial_lpmf(datamatrix[j,] | p[j]);
  }
}
