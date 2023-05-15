
data {

   int<lower=2> nSamples_t ;
   int<lower=2> nSamples_h ;

   int logE_t[nSamples_t];  // observed tx abundance
   int logE_h[nSamples_h];  // observed tx abundance

}


parameters {

   positive_ordered[2]        mu ;
   real<lower=0.001,upper=3>  phi ;
   real<lower=0,upper=1>      pi_on ;    // on/off mixing weight

}

model {

   mu ~ normal(7,5) ;

   for (s in 1:nSamples_t) {
      target += log_mix(pi_on,
                        neg_binomial_2_lpmf(logE_t[s] | mu[1], phi),
                        neg_binomial_2_lpmf(logE_t[s] | mu[2], phi));
   }

}


generated quantities {

   vector [nSamples_t] log_lik_t ;
   vector [nSamples_h] log_lik_h ;

   for (s in 1:nSamples_t) {
      log_lik_t[s] = log_mix(pi_on,
                        neg_binomial_2_lpmf(logE_t[s] | mu[1], phi),
                        neg_binomial_2_lpmf(logE_t[s] | mu[2], phi)); }

   for (s in 1:nSamples_h) {
      log_lik_h[s] = log_mix(pi_on,
                        neg_binomial_2_lpmf(logE_h[s] | mu[1], phi),
                        neg_binomial_2_lpmf(logE_h[s] | mu[2], phi)); }


}

