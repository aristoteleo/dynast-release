// This model was adapted from the NASC-seq pipeline
// https://github.com/sandberg-lab/NASC-seq/blob/master/data/NASCseqModel.stan

data {
  int<lower=0> N;
  int<lower=0> contents[N];
  int<lower=0> conversions[N];
  int<lower=0> counts[N];
  real<lower=0, upper=1> p_c;
  real<lower=0, upper=1> p_e;
}

parameters {
  real<lower=-3, upper=3> log_alpha;
  real<lower=-3, upper=3> log_beta;
  real<lower=0, upper=1> pi_g;
}

transformed parameters {
  real alpha = exp(log_alpha);
  real beta = exp(log_beta);
  real log_pi_g = log(pi_g);
  real log_pi_g_inv = log(1 - pi_g);
}

model {
  pi_g ~ beta(alpha, beta);
  for (i in 1:N) {
    target += counts[i] * log_sum_exp(
      binomial_lpmf(conversions[i] | contents[i], p_c) + log_pi_g,
      binomial_lpmf(conversions[i] | contents[i], p_e) + log_pi_g_inv
    );
  }
}
