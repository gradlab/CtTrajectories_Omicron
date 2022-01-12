if(run_pars$trapfit==TRUE){
	dist_summary <- shared_params_df %>% 
		select(-dpmeanW, -dpmeanB, -wpmeanW, -wpmeanB, -wxmeanW, -wxmeanB, -wrmeanW, -wrmeanB) %>%
			rename(dpmeanW=dpmeanW_trans,
			dpmeanB=dpmeanB_trans,
			wpmeanW=wpmeanW_trans,
			wpmeanB=wpmeanB_trans,
			wxmeanW=wxmeanW_trans,
			wxmeanB=wxmeanB_trans,
			wrmeanW=wrmeanW_trans,
			wrmeanB=wrmeanB_trans) %>% 
		summarise(
			peak.ct.NonVariant_mean=mean(global_pars[["lod"]]-dpmeanW),
			peak.ct.NonVariant_lwr=quantile(global_pars[["lod"]]-dpmeanW,0.025),
			peak.ct.NonVariant_upr=quantile(global_pars[["lod"]]-dpmeanW,0.975),
			peak.geml.NonVariant_mean=mean(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanW)),
			peak.geml.NonVariant_lwr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanW),0.025),
			peak.geml.NonVariant_upr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanW),0.975),
			proliferation.time.NonVariant_mean=mean(wpmeanW),
			proliferation.time.NonVariant_lwr=quantile(wpmeanW,0.025),
			proliferation.time.NonVariant_upr=quantile(wpmeanW,0.975),
			plateau.time.NonVariant_mean=mean(wxmeanW),
			plateau.time.NonVariant_lwr=quantile(wxmeanW,0.025),
			plateau.time.NonVariant_upr=quantile(wxmeanW,0.975),
			clearance.time.NonVariant_mean=mean(wrmeanW),
			clearance.time.NonVariant_lwr=quantile(wrmeanW,0.025),
			clearance.time.NonVariant_upr=quantile(wrmeanW,0.975),
			total.duration.NonVariant_mean=mean(wpmeanW+wxmeanW+wrmeanW),
			total.duration.NonVariant_lwr=quantile(wpmeanW+wxmeanW+wrmeanW,0.025),
			total.duration.NonVariant_upr=quantile(wpmeanW+wxmeanW+wrmeanW,0.975),
			peak.ct.Variant_mean=mean(global_pars[["lod"]]-dpmeanB),
			peak.ct.Variant_lwr=quantile(global_pars[["lod"]]-dpmeanB,0.025),
			peak.ct.Variant_upr=quantile(global_pars[["lod"]]-dpmeanB,0.975),
			peak.geml.Variant_mean=mean(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB)),
			peak.geml.Variant_lwr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB),0.025),
			peak.geml.Variant_upr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB),0.975),
			proliferation.time.Variant_mean=mean(wpmeanB),
			proliferation.time.Variant_lwr=quantile(wpmeanB,0.025),
			proliferation.time.Variant_upr=quantile(wpmeanB,0.975),
			plateau.time.Variant_mean=mean(wxmeanB),
			plateau.time.Variant_lwr=quantile(wxmeanB,0.025),
			plateau.time.Variant_upr=quantile(wxmeanB,0.975),
			clearance.time.Variant_mean=mean(wrmeanB),
			clearance.time.Variant_lwr=quantile(wrmeanB,0.025),
			clearance.time.Variant_upr=quantile(wrmeanB,0.975),
			total.duration.Variant_mean=mean(wpmeanB+wxmeanB+wrmeanB),
			total.duration.Variant_lwr=quantile(wpmeanB+wxmeanB+wrmeanB,0.025),
			total.duration.Variant_upr=quantile(wpmeanB+wxmeanB+wrmeanB,0.975)
			) %>%
		pivot_longer(everything()) %>%
		separate(name, c("parameter", "statistic"), sep="_") %>%
		pivot_wider(names_from=statistic, values_from=value) %>%
		arrange(parameter)
	} else {
	dist_summary <- shared_params_df %>% 
		select(-dpmeanW, -dpmeanB, -wpmeanW, -wpmeanB, -wrmeanW, -wrmeanB) %>%
			rename(dpmeanW=dpmeanW_trans,
			dpmeanB=dpmeanB_trans,
			wpmeanW=wpmeanW_trans,
			wpmeanB=wpmeanB_trans,
			wrmeanW=wrmeanW_trans,
			wrmeanB=wrmeanB_trans) %>% 
		summarise(
			peak.ct.NonVariant_mean=mean(global_pars[["lod"]]-dpmeanW),
			peak.ct.NonVariant_lwr=quantile(global_pars[["lod"]]-dpmeanW,0.025),
			peak.ct.NonVariant_upr=quantile(global_pars[["lod"]]-dpmeanW,0.975),
			peak.geml.NonVariant_mean=mean(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanW)),
			peak.geml.NonVariant_lwr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanW),0.025),
			peak.geml.NonVariant_upr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanW),0.975),
			proliferation.time.NonVariant_mean=mean(wpmeanW),
			proliferation.time.NonVariant_lwr=quantile(wpmeanW,0.025),
			proliferation.time.NonVariant_upr=quantile(wpmeanW,0.975),
			clearance.time.NonVariant_mean=mean(wrmeanW),
			clearance.time.NonVariant_lwr=quantile(wrmeanW,0.025),
			clearance.time.NonVariant_upr=quantile(wrmeanW,0.975),
			total.duration.NonVariant_mean=mean(wpmeanW+wrmeanW),
			total.duration.NonVariant_lwr=quantile(wpmeanW+wrmeanW,0.025),
			total.duration.NonVariant_upr=quantile(wpmeanW+wrmeanW,0.975),
			peak.ct.Variant_mean=mean(global_pars[["lod"]]-dpmeanB),
			peak.ct.Variant_lwr=quantile(global_pars[["lod"]]-dpmeanB,0.025),
			peak.ct.Variant_upr=quantile(global_pars[["lod"]]-dpmeanB,0.975),
			peak.geml.Variant_mean=mean(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB)),
			peak.geml.Variant_lwr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB),0.025),
			peak.geml.Variant_upr=quantile(convert_Ct_logGEML(global_pars[["lod"]]-dpmeanB),0.975),
			proliferation.time.Variant_mean=mean(wpmeanB),
			proliferation.time.Variant_lwr=quantile(wpmeanB,0.025),
			proliferation.time.Variant_upr=quantile(wpmeanB,0.975),
			clearance.time.Variant_mean=mean(wrmeanB),
			clearance.time.Variant_lwr=quantile(wrmeanB,0.025),
			clearance.time.Variant_upr=quantile(wrmeanB,0.975),
			total.duration.Variant_mean=mean(wpmeanB+wrmeanB),
			total.duration.Variant_lwr=quantile(wpmeanB+wrmeanB,0.025),
			total.duration.Variant_upr=quantile(wpmeanB+wrmeanB,0.975)
			) %>%
		pivot_longer(everything()) %>%
		separate(name, c("parameter", "statistic"), sep="_") %>%
		pivot_wider(names_from=statistic, values_from=value) %>%
		arrange(parameter)
	}

dist_summary$mean <- round(dist_summary$mean, 2)
dist_summary$lwr <- round(dist_summary$lwr, 2)
dist_summary$upr <- round(dist_summary$upr, 2)

