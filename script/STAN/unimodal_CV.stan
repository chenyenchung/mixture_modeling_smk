/*
 * Models log(O_gs) as arising from a single gaussian component
 *
 * 1. P(E_g) in bimodal-on v bimodal-off: mix of 2 normals
 *    1. 'on' nuclei that express gene
 *    2. 'off' nuclei that do not express gene

 * given GS knowns: O_gs, observed expression matrix (genes x samples)
 *
 * infer 5 unknowns:
 * - P(E_g|on) mu and sigma
 * - P(E_g|off) mu and sigma
 * - pi_on mixture weight
 *
 */


data {

   int<lower=2> nSamples_t ;
   int<lower=2> nSamples_h ;

   int<lower=0> logE_t[nSamples_t];  // observed tx abundance
   int<lower=0> logE_h[nSamples_h];  // observed tx abundance

}


parameters {

// gene-specific estimates of P(E_g) component 1 and 2
// GOAL: assign gene,samples  to on/off components
//   real<lower=-10,upper=17>     mu1 ;
   real<lower=0,upper=14>     mu1 ;
   real<lower=0.001,upper=3>  phi ;

}


model {

   target += neg_binomial_2_lpmf(logE_t | mu1, phi) ;
   // for (s in 1:nSamples_t) {
   //    target += neg_binomial_2_lpmf(logE_t[s] | mu1, phi) ;
   // }

}


generated quantities {

   vector [nSamples_t] log_lik_t ;
   vector [nSamples_h] log_lik_h ;

   for (s in 1:nSamples_t) {
      log_lik_t[s] = neg_binomial_2_lpmf(logE_t[s] | mu1, phi) ; }

   for (s in 1:nSamples_h) {
      log_lik_h[s] = neg_binomial_2_lpmf(logE_h[s] | mu1, phi) ; }

}
