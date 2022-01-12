if(run_pars$trapfit==TRUE){
	params <- rstan::extract(ct_fit)
	indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wx","wr"), n_indiv) %>% 
		rename(id_clean=id) %>% 
		left_join(id_map, by="id_clean") %>%
		left_join(special_map, by="id")

	shared_params_df <- make_shared_params_df(params, c("dpmeanW","wpmeanW","wxmeanW","wrmeanW","dpmeanB","wpmeanB","wxmeanB","wrmeanB","dpsd","wpsd","wxsd","wrsd")) %>% 
	mutate(dpmeanW_trans=truncnormmean(dpmeanW,dpsd,run_pars$dpmin,global_pars[["lod"]])) %>%
	mutate(dpmeanB_trans=truncnormmean(dpmeanB,dpsd,run_pars$dpmin,global_pars[["lod"]])) %>%
	mutate(wpmeanW_trans=truncnormmean(wpmeanW,wpsd,run_pars$wpmin,run_pars$wpmax)) %>%
	mutate(wpmeanB_trans=truncnormmean(wpmeanB,wpsd,run_pars$wpmin,run_pars$wpmax)) %>%
	mutate(wxmeanW_trans=truncnormmean(wxmeanW,wxsd,run_pars$wxmin,run_pars$wxmax)) %>%
	mutate(wxmeanB_trans=truncnormmean(wxmeanB,wxsd,run_pars$wxmin,run_pars$wxmax)) %>%
	mutate(wrmeanW_trans=truncnormmean(wrmeanW,wrsd,run_pars$wrmin,run_pars$wrmax)) %>%
	mutate(wrmeanB_trans=truncnormmean(wrmeanB,wrsd,run_pars$wrmin,run_pars$wrmax))

	params_df <- indiv_params_df %>% 
		left_join(shared_params_df, by="iteration") %>% 
		select(-iteration) 

	# write.csv(shared_params_df,file='output/shared_params_df.csv',row.names=FALSE)
	# write.csv(indiv_params_df,file='output/indiv_params_df.csv',row.names=FALSE)
} else {
	params <- rstan::extract(ct_fit)
	indiv_params_df <- make_indiv_params_df(params, c("tp","dp","wp","wr"), n_indiv) %>% 
		rename(id_clean=id) %>% 
		left_join(id_map, by="id_clean") %>%
		left_join(special_map, by="id")

	shared_params_df <- make_shared_params_df(params, c("dpmeanW","wpmeanW","wrmeanW","dpmeanB","wpmeanB","wrmeanB","dpsd","wpsd","wrsd")) %>% 
	mutate(dpmeanW_trans=truncnormmean(dpmeanW,dpsd,run_pars$dpmin,global_pars[["lod"]])) %>%
	mutate(dpmeanB_trans=truncnormmean(dpmeanB,dpsd,run_pars$dpmin,global_pars[["lod"]])) %>%
	mutate(wpmeanW_trans=truncnormmean(wpmeanW,wpsd,run_pars$wpmin,run_pars$wpmax)) %>%
	mutate(wpmeanB_trans=truncnormmean(wpmeanB,wpsd,run_pars$wpmin,run_pars$wpmax)) %>%
	mutate(wrmeanW_trans=truncnormmean(wrmeanW,wrsd,run_pars$wrmin,run_pars$wrmax)) %>%
	mutate(wrmeanB_trans=truncnormmean(wrmeanB,wrsd,run_pars$wrmin,run_pars$wrmax))

	params_df <- indiv_params_df %>% 
		left_join(shared_params_df, by="iteration") %>% 
		select(-iteration) 

	# write.csv(shared_params_df,file='output/shared_params_df.csv',row.names=FALSE)
	# write.csv(indiv_params_df,file='output/indiv_params_df.csv',row.names=FALSE)
}