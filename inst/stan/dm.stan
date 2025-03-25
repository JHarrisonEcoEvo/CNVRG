data {
  int<lower=1> N;
  int<lower=1> nreps;
  int<lower=1> notus;

  array[N] int<lower=1> start;
  array[N] int<lower=1> end;

  array[nreps, notus] int datamatrix;
}

parameters {
  array[N] real<lower=0> theta;
  array[N] simplex[notus] pi;
  array[nreps] simplex[notus] p;
}

model {
  for(i in 1:N){
    target += exponential_lpdf(theta[i] | 0.001);
    target += dirichlet_lpdf(pi[i] | rep_vector(0.0000001, notus));
    for(j in start[i]:end[i]){
      target += dirichlet_lpdf(p[j] | theta[i]*pi[i]);
      target += multinomial_lpmf(datamatrix[j,] | p[j]);
    }
  }
}
