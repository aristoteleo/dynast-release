data {
  int<lower=0> N;
  int<lower=0> contents[N];
  int<lower=0> conversions[N];
  real<lower=0, upper=1> p_c;
  real<lower=0, upper=1> p_e;
}

parameters {
  real<lower=-4, upper=4> log_alpha;
  real<lower=-4, upper=4> log_beta;
  real<lower=0, upper=1> pi_g;
}

transformed parameters{
    real alpha;
    real beta;
    alpha = exp(log_alpha);
    beta = exp(log_beta);
}

model {
    pi_g ~ beta(alpha, beta);
    for (i in 1:N){
        target += log_sum_exp(binomial_lpmf(conversions[i] | contents[i], p_c) + bernoulli_lpmf(1|pi_g),binomial_lpmf(conversions[i] | contents[i],p_e) + bernoulli_lpmf(0|pi_g));
    }
}
