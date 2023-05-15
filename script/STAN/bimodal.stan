
data {

   int<lower=2> nSamples ;
   int<lower=1> nDrivers ;
   int<lower=1> nCells ;

   int logE[nSamples];       // observed tx abundance
   int                         driver[nSamples] ;     // driver
   int                         cell[nSamples] ;       // cell type

}


parameters {

// gene-specific estimates of P(E_g) component 1 and 2
// GOAL: assign gene,samples  to on/off components

   positive_ordered[2]        mu ;       // component means
   real<lower=0.001,upper=3>  phi ;      // equal overdispersion
   real<lower=0,upper=1>      pi_on ;    // on/off mixing weight

}

transformed parameters {

   real mu_on ;
   real mu_off ;

   mu_on  = mu[2] ;
   mu_off = mu[1] ;

}


model {

   mu ~ normal(7,5) ;

   for (s in 1:nSamples) {
      target += log_mix(pi_on,
                        neg_binomial_2_lpmf(logE[s] | mu_on, phi),
                        neg_binomial_2_lpmf(logE[s] | mu_off, phi));
   }

}


generated quantities {

   vector [nSamples]    pon_gs ;   // posterior likelihood of expression
   vector [nDrivers]    pon_gd ;   // cell-level likelihood of expression
   vector [nCells]      pon_gc ;   // cell-level likelihood of expression

   vector [nDrivers]    pon_gd_onmass ;
   vector [nDrivers]    pon_gd_offmass ;

   vector [nCells]      pon_gc_onmass ;
   vector [nCells]      pon_gc_offmass ;


   for (d in 1:nDrivers) {
      pon_gd_onmass[d] = 0 ;
      pon_gd_offmass[d] = 0 ;
   }

   for (c in 1:nCells) {
      pon_gc_onmass[c] = 0 ;
      pon_gc_offmass[c] = 0 ;
   }

   for (s in 1:nSamples) {

      real ps[2];  // temp for log component P
      real sumps;  // temp for sum of log component P

      ps[1] = neg_binomial_2_lpmf(logE[s] | mu_on, phi) ;
      ps[2] = neg_binomial_2_lpmf(logE[s] | mu_off, phi) ;

      sumps = log_mix(pi_on, ps[1], ps[2]) ;

      pon_gs[s] = log(pi_on) + ps[1] - sumps ;

      pon_gd_onmass[driver[s]]  = pon_gd_onmass[driver[s]] + ps[1] ;
      pon_gd_offmass[driver[s]] = pon_gd_offmass[driver[s]] + ps[2] ;

      pon_gc_onmass[cell[s]]  = pon_gc_onmass[cell[s]] + ps[1] ;
      pon_gc_offmass[cell[s]] = pon_gc_offmass[cell[s]] + ps[2] ;

   }

   for (c in 1:nCells) {
      pon_gc[c] = log(pi_on) + pon_gc_onmass[c] -
                  log_mix(pi_on, pon_gc_onmass[c], pon_gc_offmass[c]) ;
   }

   for (d in 1:nDrivers) {
      pon_gd[d] = log(pi_on) + pon_gd_onmass[d] -
                  log_mix(pi_on, pon_gd_onmass[d], pon_gd_offmass[d]) ;
   }

}

