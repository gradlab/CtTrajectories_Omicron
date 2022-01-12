# ============================================================================
if(run_pars$trapfit==TRUE){ # Trapezoid fit -----------------------------------
# ============================================================================

ct_model <- stan_model("code/fit_posteriors_special_trap.stan") 

fit_startq <- Sys.time()
ct_fit <- sampling(ct_model, 
	data=list(
		N=nrow(indiv_data), 
		n_id=length(unique(indiv_data$id_clean)),
		lod=global_pars[["lod"]], 
		id=indiv_data$id_clean,
		special=as.list(prior_pars)$special,
		t=indiv_data$t, 
		y=indiv_data$y, 
		tpsd=as.list(prior_pars)$tpsd,
		dpmin=as.list(prior_pars)$dpmin,
		dpmean_prior=as.list(prior_pars)$dpmean_prior,
		dpsd_prior=as.list(prior_pars)$dpsd_prior,
		wpmin=as.list(prior_pars)$wpmin,
		wpmax=as.list(prior_pars)$wpmax,
		wpmean_prior=as.list(prior_pars)$wpmean_prior,
		wpsd_prior=as.list(prior_pars)$wpsd_prior,
		wxmin=as.list(prior_pars)$wxmin,
		wxmax=as.list(prior_pars)$wxmax,
		wxmean_prior=as.list(prior_pars)$wxmean_prior,
		wxsd_prior=as.list(prior_pars)$wxsd_prior,
		wrmin=as.list(prior_pars)$wrmin,
		wrmax=as.list(prior_pars)$wrmax,
		wrmean_prior=as.list(prior_pars)$wrmean_prior,
		wrsd_prior=as.list(prior_pars)$wrsd_prior,
		sigma_max=as.list(prior_pars)$sigma_max,
		sigma_prior_scale=as.list(prior_pars)$sigma_prior_scale,
		lambda=as.list(prior_pars)$lambda,
		fpmean=as.list(prior_pars)$fpmean), 
	iter=1000, chains=4) # iter=2000
# , control = list(adapt_delta=0.85)
# , control = list(adapt_delta=0.99)
# control = list(adapt_delta=0.95, max_treedepth=15)
fit_endq <- Sys.time()
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))

# launch_shinystan_nonblocking(ct_fit)

# ============================================================================
} else { # --- Regular fit ---------------------------------------------------
# ============================================================================

ct_model <- stan_model("code/fit_posteriors_special.stan") 
# ct_model <- stan_model("code/fit_posteriors_special_lognorm.stan") 

fit_startq <- Sys.time()
ct_fit <- sampling(ct_model, 
	data=list(
		N=nrow(indiv_data), 
		n_id=length(unique(indiv_data$id_clean)),
		lod=global_pars[["lod"]], 
		id=indiv_data$id_clean,
		special=as.list(prior_pars)$special,
		t=indiv_data$t, 
		y=indiv_data$y, 
		tpsd=as.list(prior_pars)$tpsd,
		dpmin=as.list(prior_pars)$dpmin,
		dpmean_prior=as.list(prior_pars)$dpmean_prior,
		dpsd_prior=as.list(prior_pars)$dpsd_prior,
		wpmin=as.list(prior_pars)$wpmin,
		wpmax=as.list(prior_pars)$wpmax,
		wpmean_prior=as.list(prior_pars)$wpmean_prior,
		wpsd_prior=as.list(prior_pars)$wpsd_prior,
		wrmin=as.list(prior_pars)$wrmin,
		wrmax=as.list(prior_pars)$wrmax,
		wrmean_prior=as.list(prior_pars)$wrmean_prior,
		wrsd_prior=as.list(prior_pars)$wrsd_prior,
		sigma_max=as.list(prior_pars)$sigma_max,
		sigma_prior_scale=as.list(prior_pars)$sigma_prior_scale,
		lambda=as.list(prior_pars)$lambda,
		fpmean=as.list(prior_pars)$fpmean), 
	iter=1000, chains=4) # iter=2000
# , control = list(adapt_delta=0.85)
# , control = list(adapt_delta=0.99)
# control = list(adapt_delta=0.95, max_treedepth=15)
fit_endq <- Sys.time()
print(paste0("Fit time: ",difftime(fit_endq, fit_startq, units="min")," mins"))

# launch_shinystan_nonblocking(ct_fit)

}

