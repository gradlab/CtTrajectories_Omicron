functions {
  // mufun returns dCt, the delta Ct below the LOD:
  real mufun(real t, real tp, real wp, real wr, real dp){
    // Viral load rises between onset and peak: 
    if(t>(tp-wp) && t<=tp)
      return(dp/wp*(t-(tp-wp)));
    // Viral load falls between peak and recovery: 
    else if(t>tp && t<=(tp+wr))
      return(dp - dp/wr*(t-tp));
    // Ct = lod after recovery: 
    else 
      return(0);
  }
}

data {
  int<lower=0> N;             // Number of concatenated data points
  int<lower=0> n_id;          // Number of individuals 
  real<lower=0> lod;          // Limit of detection
  int<lower=0> id[N];         // Vector marking which datum belongs to which id
  int<lower=0> special[n_id];    // Vector marking special ids
  real t[N];                  // Vector marking the time for each data point 
  real<lower=0, upper=lod> y[N];  // Concatenated data vector 
  real<lower=0> tpsd;         // Prior sd for the onset time (days)
  real<lower=0> dpmin;        // Min peak Ct (delta from lod)
  real<lower=0> dpmean_prior; // Prior mean peak Ct (delta from lod)
  real<lower=0> dpsd_prior;   // Prior sd peak Ct (delta from lod)
  real<lower=0> wpmin;        // Min proliferation duration
  real<lower=0> wpmax;        // Max proliferation duration
  real<lower=0> wpmean_prior; // Prior mean proliferation duration
  real<lower=0> wpsd_prior;   // Prior sd proliferation duration
  real<lower=0> wrmin;        // Min clearance duration
  real<lower=0> wrmax;        // Max clearance duration
  real<lower=0> wrmean_prior; // Prior mean clearance duration 
  real<lower=0> wrsd_prior;   // Prior sd clearance duration 
  real<lower=0> sigma_max;    // Max allowable value for observation noise
  real<lower=0> sigma_prior_scale;  // Prior observation noise Cauchy scale
  real<lower=0, upper=1> lambda; // Mixture probability (~1-sensitivity)
  real<lower=0> fpmean;        // False positive mean Ct
}

transformed data {
  real<lower=0, upper=lod> ydrop[N];  // Concatenated deviation from LOD 

  real loglambda;
  real log1mlambda;

  real dpcauchypriorscale;
  real wpcauchypriorscale;
  real wrcauchypriorscale;

  // real sigma;

  for(i in 1:N){
    ydrop[i] = lod-y[i];
  }

  loglambda = log(lambda);
  log1mlambda = log1m(lambda);

  // Define cauchy prior scales so that 90% of the half-distribution mass lies below the max cutoff for that parameter. 
  dpcauchypriorscale = lod/tan(pi()*(0.95-0.5));
  wpcauchypriorscale = wpmax/tan(pi()*(0.95-0.5));
  wrcauchypriorscale = wrmax/tan(pi()*(0.95-0.5));

}

parameters {

  real<lower=dpmin, upper=lod> dpmeanB;   // Poplation peak Ct drop mean
  real<lower=wpmin, upper=wpmax> wpmeanB; // Population onset-to-peak time mean
  real<lower=wrmin, upper=wrmax> wrmeanB; // Population peak-to-recovery time mean 

  real<lower=dpmin, upper=lod> dpmeanW;   // Poplation peak Ct drop mean
  real<lower=wpmin, upper=wpmax> wpmeanW; // Population onset-to-peak time mean
  real<lower=wrmin, upper=wrmax> wrmeanW; // Population peak-to-recovery time mean 

  real<lower=0> dpsd;          // Poplation peak Ct drop sd
  real<lower=0> wpsd;          // Population onset-to-peak time sd
  real<lower=0> wrsd;          // Population peak-to-recovery time sd

  real tp[n_id];                        // Peak time
  real<lower=dpmin, upper=lod> dp[n_id];    // Peak Ct drop
  real<lower=wpmin, upper=wpmax> wp[n_id];  // Onset-to-peak time
  real<lower=wrmin, upper=wrmax> wr[n_id];  // Peak-to-recovery time 

  real<lower=0, upper=sigma_max> sigma;    // Process noise during infection
}

transformed parameters {
  real process_sd[N];
  real mu[N];
  
  for(i in 1:N){
    mu[i]=mufun(t[i], tp[id[i]], wp[id[i]], wr[id[i]], dp[id[i]]);
    process_sd[i]=sigma;
  };
}

model {

  // Hierarchical priors:
  dpmeanB ~ normal(dpmean_prior,dpsd_prior) T[dpmin,lod];
  dpmeanW ~ normal(dpmean_prior,dpsd_prior) T[dpmin,lod];
  
  wpmeanB ~ normal(wpmean_prior, wpsd_prior) T[wpmin,wpmax];
  wpmeanW ~ normal(wpmean_prior, wpsd_prior) T[wpmin,wpmax];

  wrmeanB ~ normal(wrmean_prior, wrsd_prior) T[wrmin,wrmax];
  wrmeanW ~ normal(wrmean_prior, wrsd_prior) T[wrmin,wrmax];

  dpsd ~ cauchy(0,dpcauchypriorscale) T[0,];
  wpsd ~ cauchy(0,wpcauchypriorscale) T[0,];
  wrsd ~ cauchy(0,wrcauchypriorscale) T[0,];  

  sigma ~ cauchy(0,sigma_prior_scale) T[0,10];

  // Individual parameter specifications:
  tp ~ normal(0,tpsd);
  for(i in 1:n_id){
    if(special[i] == 1){
      dp[i] ~ normal(dpmeanB, dpsd) T[dpmin, lod];
      wp[i] ~ normal(wpmeanB, wpsd) T[wpmin, wpmax];
      wr[i] ~ normal(wrmeanB, wrsd) T[wrmin, wrmax];
      } else {
      dp[i] ~ normal(dpmeanW, dpsd) T[dpmin, lod];
      wp[i] ~ normal(wpmeanW, wpsd) T[wpmin, wpmax];
      wr[i] ~ normal(wrmeanW, wrsd) T[wrmin, wrmax];
      }
  } 

  // Main model specification: 
  for(i in 1:N){

    target += log_sum_exp(
      log1mlambda + normal_lpdf(ydrop[i] | mu[i], process_sd[i]),
      loglambda + exponential_lpdf(ydrop[i] | 1/fpmean));

    if (ydrop[i] < 0 || ydrop[i] > lod)
      target += negative_infinity();
    else
      target += -log_diff_exp(
        log1mlambda + normal_lcdf(lod | mu[i], process_sd[i]),
        log1mlambda + normal_lcdf(0 | mu[i], process_sd[i]));
    // see https://mc-stan.org/docs/2_18/reference-manual/sampling-statements-section.html

    }

}


