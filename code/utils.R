library(tidyverse) 

clean_person_id <- function(ct_dat){
	id_list <- sort(unique(ct_dat$PersonID))
	id_df <- data.frame(PersonID=id_list, PersonIDClean=1:length(id_list), stringsAsFactors=FALSE)
	out <- left_join(ct_dat, id_df, by="PersonID")
	return(out)
}

clean_infection_events <- function(ct_dat){
	ie_list <- sort(unique(ct_dat$InfectionEvent))
	ie_df <- data.frame(InfectionEvent=ie_list, InfectionEventClean=1:length(ie_list), stringsAsFactors=FALSE)
	out <- left_join(ct_dat, ie_df, by="InfectionEvent")
	return(out)
}

# ct_dat_refined %>% 
# 	clean_infection_events %>% 
# 	select(PersonID, TestDate, CtT1, InfectionEvent, InfectionEventClean) %>% 
# 	print(n=Inf)

save_figlist <- function(obj,dir,name,driver,width=8,height=5){
	mapply(function(x,y) ggsave(x,
		file=paste0(dir,name,"_",y,".",driver), width=width, height=height),
	obj,1:length(obj))
}

get_infection_events <- function(df){
	df <- df %>% 
		split(.$PersonID) %>% 
		map(~ get_infection_events_indiv(.)) %>% 
		bind_rows() %>% 
		arrange(PersonID, TestDate) %>% 
		mutate(InfectionEvent=cumsum(FirstKeep)) %>% 
		select(-FirstKeep)
}

get_infection_events_indiv <- function(df_person, winlen=14){

	df_person <- df_person %>% 
		mutate(Trim=0)

	for(r in 1:nrow(df_person)){
	ctsum <- df_person %>% 
		filter(
			(TestDateIndex >= df_person[r,]$TestDateIndex-winlen) & 
			(TestDateIndex <= df_person[r,]$TestDateIndex+winlen))	 %>% 
			summarise(CtSum=sum(CtT1-global_pars[["lod"]])) %>% 
			pull(CtSum)
	if(ctsum==0){
		df_person[r,]$Trim <- 1
		}
	}

	df_person <- df_person %>% 
		mutate(Keep=1-Trim) %>% 
		mutate(FirstKeep=case_when(
			# (Keep==1 & (lag(Keep)==0 | is.na(lag(Keep))))~1, TRUE~0)) %>% 
			(Keep==1 & (lag(Keep)==0 | is.na(lag(Keep)) | lag(TestDateIndex)<(TestDateIndex-winlen)))~1, TRUE~0)) %>% 
			mutate(InfectionEvent=cumsum(FirstKeep)) %>% 
			filter(Keep==1) %>% 
		select(-Trim,-Keep) 
}

# ct_dat %>% 
# 	filter(Lineage!="None") %>% 
# 	group_by(PersonID) %>% 
# 	summarise(nl=n_distinct(Lineage)) %>% 
# 	arrange(desc(nl))

propagate_lineage <- function(ct_dat){
	# Error control: 
	linreps <- ct_dat %>% 
		filter(Lineage!="None") %>% 
		group_by(PersonID) %>% 
		summarise(nl=n_distinct(Lineage)) %>% 
		arrange(desc(nl))

	maxlinreps <- max(linreps$nl) 
	if(maxlinreps>1){
		print("Error: some people have multiple lineages")
		return(NA)
	}

	# Propagate lineages: 
	linjoin <- ct_dat %>% 
		mutate(Lineage=case_when(Lineage=="None"~"ZZZNone", TRUE~Lineage)) %>% 
		group_by(PersonID) %>% 
		arrange(Lineage) %>% 
		select(PersonID, Lineage) %>% 
		slice(1) %>% 
		mutate(Lineage=case_when(Lineage=="ZZZNone"~"None", TRUE~Lineage)) 

	out <- ct_dat %>% 
		select(-Lineage) %>% 
		left_join(linjoin)

	return(out) 
}

# filter(temp, PersonID==1717) %>% select(PersonID, Lineage)
make_test_date_index <- function(ct_dat){
	minctdf <- ct_dat %>% 
		arrange(PersonID, TestDate) %>%
		group_by(PersonID) %>%
		mutate(index=1:n()) %>%
		filter(CtT1==min(CtT1)) %>%
		slice(1) %>%
		select(PersonID, TestDate) %>%
		rename(MinCtTestDate=TestDate)

	out <- ct_dat %>% 
		left_join(minctdf, by="PersonID") %>%
		mutate(TestDateIndex=as.numeric(difftime(TestDate,MinCtTestDate,units="days"))) %>%
		select(-MinCtTestDate)

	return(out)
}

make_test_date_index_infevent <- function(ct_dat){
	minctdf <- ct_dat %>% 
		arrange(PersonID, TestDate) %>%
		group_by(PersonID,InfectionEvent) %>%
		mutate(index=1:n()) %>%
		filter(CtT1==min(CtT1)) %>%
		slice(1) %>%
		select(PersonID,InfectionEvent,TestDate) %>%
		rename(MinCtTestDate=TestDate)

	out <- ct_dat %>% 
		select(-TestDateIndex) %>% 
		left_join(minctdf, by=c("PersonID","InfectionEvent")) %>%
		mutate(TestDateIndex=as.numeric(difftime(TestDate,MinCtTestDate,units="days"))) %>%
		select(-MinCtTestDate)

	return(out)
}

# trim_negatives <- function(indiv_data, global_pars){
# 	# lod <- 40
# 	out <- indiv_data %>% 
# 		split(.$id) %>% 
# 		map(~ arrange(., t)) %>% 
# 		map(~ mutate(., rowindex=1:n())) %>% 
# 		map(~ mutate(., ispositive=case_when(y<global_pars[["lod"]] ~ 1, TRUE~0))) %>% 
# 		map(~ mutate(., ispositive_lag=lag(ispositive))) %>%
# 		map(~ mutate(., ispositive_lag2=lag(ispositive,2))) %>%
# 		map(~ mutate(., ispositive_lead=lead(ispositive))) %>%
# 		map(~ mutate(., ispositive_lead2=lead(ispositive,2))) %>%
# 		map(~ filter(., ispositive==1 | ispositive_lag==1 | ispositive_lag2==1 | ispositive_lead==1 | ispositive_lead2==1)) %>%
# 		map(~ select(., id, id_clean, t, y, special)) %>%
# 		bind_rows()
# 	return(out)
# }

trim_negatives <- function(indiv_data, global_pars){
	# lod <- 40
	out <- indiv_data %>% 
		split(.$id) %>% 
		map(~ arrange(., t)) %>% 
		map(~ mutate(., rowindex=1:n())) %>% 
		map(~ mutate(., ispositive=case_when(y<global_pars[["lod"]] ~ 1, TRUE~0))) %>% 
		map(~ mutate(., ispositive_lag=lag(ispositive))) %>%
		map(~ mutate(., ispositive_lag2=lag(ispositive,2))) %>%
		map(~ mutate(., ispositive_lead=lead(ispositive))) %>%
		map(~ mutate(., ispositive_lead2=lead(ispositive,2))) %>%
		map(~ filter(., ispositive==1 | ispositive_lag==1 | ispositive_lag2==1 | ispositive_lead==1 | ispositive_lead2==1)) %>%
		map(~ select(., -rowindex, -ispositive, -ispositive_lag, -ispositive_lag2, -ispositive_lead, -ispositive_lead2)) %>%
		bind_rows()
	return(out)
}

launch_shinystan_nonblocking <- function(fit) {
  library(future)
  plan(multisession)
  future(
    launch_shinystan(fit) #You can replace this with any other Shiny app
  )
}

sample_n_groups <- function(grouped_df, size, replace=FALSE, weight=NULL){
	# From https://cmdlinetips.com/2019/07/how-to-randomly-select-groups-in-r-with-dplyr/
	grp_var <- grouped_df %>% 
	    groups %>%
	    unlist %>% 
	    as.character

	random_grp <- grouped_df %>% 
	    summarise() %>% 
	    slice_sample(n=size, replace=replace, weight_by=weight)# %>% 
	    # mutate(unique_id = 1:n())

	out <- grouped_df %>% 
	    right_join(random_grp, by=grp_var)

	return(out)
}

# For generating names for a matrix-turned-data frame: 
makenames <- function(parname, n_indiv){
	unlist(lapply(1:n_indiv, function(x) paste0(parname, "_", x)))
}

# For parsing the 'parameters' output from Stan: 
parseparam <- function(extracted_params, parname, n_indiv){
	as_tibble(setNames(as.data.frame(extracted_params[[parname]]), makenames(parname,n_indiv)))
}

make_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value") %>%
		select(-iteration) 
}

make_indiv_params_df <- function(extracted_params, parnames, n_indiv){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) parseparam(extracted_params,x,n_indiv)), cbind) %>% 
		as_tibble %>% 
		mutate(iteration=1:n()) %>% 
		pivot_longer(-iteration) %>% 
		separate(name, c("param","id"), sep="_") %>%
		pivot_wider(c("id","iteration"), names_from="param", values_from="value")
}

make_shared_params_df <- function(extracted_params, parnames){
	# Use "reduce" here
	out <- reduce(lapply(parnames, function(x) 
		as_tibble(setNames(as.data.frame(extracted_params[[x]]),x))
		), cbind) %>%
		as_tibble() %>%
		mutate(iteration=1:n())
}

plot_ct_fit <- function(params_df, global_pars, indiv_data, ctalpha=0.01,specialcol="red", basecol="blue", xlim=c(NA,NA), vlabel="Variant", nvlabel="Non-Variant",ntraces=100){
	with(as.list(global_pars),{
	params_df %>% 
		mutate(id_clean=as.integer(id_clean)) %>%
		group_by(id_clean) %>%
		sample_n(ntraces) %>% 
		ungroup() %>%
		ggplot() + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
			geom_point(data=indiv_data, aes(x=t, y=y,col=factor(special)), size=0.5) + 
			scale_color_manual(values=c("1"=specialcol,"0"=basecol), labels=c("1"=vlabel,"0"=nvlabel)) + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
			labs(x="Time since min Ct (days)", y="Ct") + 
			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
			scale_x_continuous(limits=xlim) + 
			facet_wrap(~id) # Change to id_clean to obscure identities
			})
}


plot_ct_fit_split <- function(params_df, global_pars, indiv_data, ctalpha=0.01,specialcol="red", basecol="blue", xlim=c(NA,NA), vlabel="Variant", nvlabel="Non-Variant",nsplits=4,shade_ids=c(),ntraces=100){
	with(as.list(global_pars),{

	shade_df <- tibble(id=sort(unique(params_df$id))) %>% 
		mutate(shade=case_when(id %in% shade_ids~"Yes", TRUE~"No"))

	params_df <- params_df %>% 
		mutate(id_clean=as.integer(id_clean)) %>%
		mutate(splitval=ceiling((id_clean)/(max(id_clean)/nsplits)))

	splitmap <- params_df %>% 
		group_by(id) %>% 
		summarise(splitval=first(splitval)) 

	indiv_data_split <- indiv_data %>% 
		left_join(splitmap,by="id") %>% 
		split(.$splitval)

	shade_df_split <- shade_df %>% 
		left_join(splitmap,by="id") %>% 
		split(.$splitval)


	params_df %>% 
		split(.$splitval) %>% 
		map(~ group_by(.,id_clean)) %>%
		map(~ sample_n(.,ntraces)) %>% 
		map(~ ungroup(.)) %>%
		imap(~ ggplot(.x) + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
			geom_point(data=indiv_data_split[[.y]], aes(x=t, y=y,col=factor(special)), size=0.5) + 
			scale_color_manual(values=c("1"=specialcol,"0"=basecol), labels=c("1"=vlabel,"0"=nvlabel)) + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
			labs(x="Time since min Ct (days)", y="Ct") + 
			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
			scale_x_continuous(limits=xlim) + 
			geom_rect(data=shade_df_split[[.y]],aes(fill=shade),xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf, col="white", alpha=0.2) + 
			scale_fill_manual(values=c("Yes"="Gray","No"="White"), guide="none") + 
			facet_wrap(~id)) # Change to id_clean to obscure identities
			})
}


plot_ct_fit_trap <- function(params_df, global_pars, indiv_data, ctalpha=0.01, specialcol="red", xlim=c(NA,NA), ntraces=100){
	with(as.list(global_pars),{
	params_df %>% 
		mutate(id_clean=as.integer(id_clean)) %>%
		group_by(id_clean) %>%
		sample_n(ntraces) %>% 
		ungroup() %>%
		ggplot() + 
			# Plot traces:
			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wx, yend=lod-dp), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wx, y=lod-dp, xend=tp+wx+wr, yend=lod), alpha=ctalpha) + 
			geom_segment(aes(x=tp+wx+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
			# Plot data:
			geom_point(data=indiv_data, aes(x=t, y=y,col=factor(special)), size=0.5) + 
			scale_color_manual(values=c("1"=specialcol,"0"="blue"), labels=c("1"="Variant","0"="Non-variant")) + 
			theme_minimal() + 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
			labs(x="Time since min Ct (days)", y="Ct") + 
			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
			scale_x_continuous(limits=xlim) + 
			facet_wrap(~id) # Change to id_clean to obscure identities
			})
}

# plot_ct_fit_b117_trap_rawid <- function(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=100){
# 	with(as.list(global_pars),{
# 	params_df %>% 
# 		filter(id==3229) %>%
# 		mutate(id=as.integer(id)) %>%
# 		group_by(id) %>%
# 		sample_n(ntraces) %>% 
# 		ungroup() %>%
# 		ggplot() + 
# 			# Plot traces:
# 			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wx, yend=lod-dp), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp+wx, y=lod-dp, xend=tp+wx+wr, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp+wx+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
# 			# Plot data:
# 			geom_point(data=filter(indiv_data, id=="3229"), aes(x=t, y=y, col=factor(b117)), size=0.5) + 
# 			scale_color_manual(values=c("1"="red","0"="blue"), labels=c("1"="B117","0"="non-B117")) + 
# 			theme_minimal() + 
# 			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
# 			labs(x="Time since min Ct (days)", y="Ct") + 
# 			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
# 			# scale_x_continuous(limits=c(-7,14)) + 
# 			facet_wrap(~id)
# 			})
# }

# plot_ct_fit_b117_fulldata <- function(params_df, global_pars, ct_dat_refined, ctalpha=0.01, ntraces=100){
# 	with(as.list(global_pars),{
# 	params_df %>% 
# 		mutate(id_clean=as.integer(id_clean)) %>%
# 		group_by(id_clean) %>%
# 		sample_n(ntraces) %>% 
# 		ungroup() %>%
# 		ggplot() + 
# 			# Plot traces:
# 			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
# 			# Plot data:
# 			geom_point(data=mutate(ct_dat_refined,id_clean=as.integer(PersonIDClean)), aes(x=TestDateIndex, y=CtT1,col=factor(B117Status)), size=0.5) + 
# 			scale_color_manual(values=c("1"="red","0"="blue","Yes"="red","No"="blue"), labels=c("1"="B117","0"="non-B117","Yes"="B117","No"="non-B117")) + 
# 			theme_minimal() + 
# 			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8)) + 
# 			labs(x="Time since min Ct (days)", y="Ct") + 
# 			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
# 			facet_wrap(~id_clean)
# 			})
# }

# plot_ct_fit_b117_noid <- function(params_df, global_pars, indiv_data, ctalpha=0.01, ntraces=100){
# 	with(as.list(global_pars),{
# 	params_df %>% 
# 		mutate(id_clean=as.integer(id_clean)) %>%
# 		group_by(id_clean) %>%
# 		sample_n(ntraces) %>% 
# 		ungroup() %>%
# 		ggplot() + 
# 			# Plot traces:
# 			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
# 			# Plot data:
# 			geom_point(data=indiv_data, aes(x=t, y=y,col=factor(b117)), size=0.5) + 
# 			scale_color_manual(values=c("1"="red","0"="blue"), labels=c("1"="B117","0"="non-B117")) + 
# 			theme_minimal() + 
# 			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8), strip.text.x=element_blank()) + 
# 			labs(x="Time since min Ct (days)", y="Ct") + 
# 			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
# 			facet_wrap(~b117+id_clean)
# 			})
# }

# plot_ct_fit_b117_fulldata_noid <- function(params_df, global_pars, ct_dat_refined, ctalpha=0.01, ntraces=100){
# 	with(as.list(global_pars),{
# 	params_df %>% 
# 		mutate(id_clean=as.integer(id_clean)) %>%
# 		group_by(id_clean) %>%
# 		sample_n(ntraces) %>% 
# 		ungroup() %>%
# 		ggplot() + 
# 			# Plot traces:
# 			geom_segment(aes(x=-Inf, y=lod, xend=tp-wp, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp-wp, y=lod, xend=tp, yend=lod-dp), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp, y=lod-dp, xend=tp+wr, yend=lod), alpha=ctalpha) + 
# 			geom_segment(aes(x=tp+wr, y=lod, xend=Inf, yend=lod), alpha=ctalpha) + 
# 			# Plot data:
# 			geom_point(data=mutate(ct_dat_refined,
# 				id_clean=as.integer(PersonIDClean),
# 				b117=case_when(B117Status=="Yes"~1, TRUE~0)), aes(x=TestDateIndex, y=CtT1,col=factor(B117Status)), size=0.5) + 
# 			scale_color_manual(values=c("1"="red","0"="blue","Yes"="red","No"="blue"), labels=c("1"="B117","0"="non-B117","Yes"="B117","No"="non-B117")) + 
# 			theme_minimal() + 
# 			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), axis.text.y=element_text(size=8), strip.text.x=element_blank()) + 
# 			labs(x="Time since min Ct (days)", y="Ct") + 
# 			scale_y_reverse(breaks=c(40,30,20,10), labels=c("(-)","30","20","10"), sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + 
# 			facet_wrap(~b117+id_clean)
# 			})
# }

grid_off <- list(theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()))

y_ticks_off <- list(theme(axis.ticks.y=element_blank(),axis.text.y=element_blank()))

convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
	out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
	return(out) 
}

srise <- function(x, dp, wp){
	out <- dp*(1+x/wp)
	return(out)
}
sfall <- function(x, dp, wr){
	out <- dp*(1-x/wr)
	return(out)
}

sfall_trap <- function(x, dp, wx, wr){
	out <- -dp/wr*x + dp*(wx/wr + 1)
	out[x <= wx] <- dp[x <= wx]
	# if(x <= wx){
	# 	out <- dp
	# 	} else {
	# 	out <- -dp/wr*x + dp*(wx/wr + 1)
	# 	}
	return(out)
}

make_sample_trajectory_special <- function(shared_params_df, global_pars, siglevel=0.95, ge=FALSE, referencefill="blue", specialfill="red"){
	# For asymptomatic:
	with(as.list(global_pars),{

	shared_params_df <- shared_params_df %>% 
		select(-dpmeanW, -dpmeanB, -wpmeanW, -wpmeanB, -wrmeanW, -wrmeanB) %>%
		rename(dpmeanW=dpmeanW_trans,
		dpmeanB=dpmeanB_trans,
		wpmeanW=wpmeanW_trans,
		wpmeanB=wpmeanB_trans,
		wrmeanW=wrmeanW_trans,
		wrmeanB=wrmeanB_trans)

	wp_mean_W <- mean(shared_params_df$wpmeanW)
	wp_lwr_W <- quantile(shared_params_df$wpmeanW,(1-siglevel)/2)
	wp_upr_W <- quantile(shared_params_df$wpmeanW,1-(1-siglevel)/2)

	wr_mean_W <- mean(shared_params_df$wrmeanW)
	wr_lwr_W <- quantile(shared_params_df$wrmeanW,(1-siglevel)/2)
	wr_upr_W <- quantile(shared_params_df$wrmeanW,1-(1-siglevel)/2)

	dp_mean_W <- mean(shared_params_df$dpmeanW)
	dp_lwr_W <- quantile(shared_params_df$dpmeanW,(1-siglevel)/2)
	dp_upr_W <- quantile(shared_params_df$dpmeanW,1-(1-siglevel)/2)

	xvals_proliferation_W <- seq(from=-wp_upr_W, 0, length.out=500)
	xvals_clearance_W <- seq(from=0, wr_upr_W, length.out=500)

	yvals_upr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),(1-siglevel)/2)))
	yvals_upr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),1-(1-siglevel)/2)))
	yvals_lwr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),(1-siglevel)/2)))

	# For symptomatic:
	wp_mean_B <- mean(shared_params_df$wpmeanB)
	wp_lwr_B <- quantile(shared_params_df$wpmeanB,(1-siglevel)/2)
	wp_upr_B <- quantile(shared_params_df$wpmeanB,1-(1-siglevel)/2)

	wr_mean_B <- mean(shared_params_df$wrmeanB)
	wr_lwr_B <- quantile(shared_params_df$wrmeanB,(1-siglevel)/2)
	wr_upr_B <- quantile(shared_params_df$wrmeanB,1-(1-siglevel)/2)

	dp_mean_B <- mean(shared_params_df$dpmeanB)
	dp_lwr_B <- quantile(shared_params_df$dpmeanB,(1-siglevel)/2)
	dp_upr_B <- quantile(shared_params_df$dpmeanB,1-(1-siglevel)/2)

	xvals_proliferation_B <- seq(from=-wp_upr_B, 0, length.out=500)
	xvals_clearance_B <- seq(from=0, wr_upr_B, length.out=500)

	yvals_upr_proliferation_B <- unlist(lapply(xvals_proliferation_B, 
	function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_B <- unlist(lapply(xvals_proliferation_B, 
	function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),(1-siglevel)/2)))
	yvals_upr_clearance_B <- unlist(lapply(xvals_clearance_B, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),1-(1-siglevel)/2)))
	yvals_lwr_clearance_B <- unlist(lapply(xvals_clearance_B, 
	function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),(1-siglevel)/2)))

	if(ge==FALSE){
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=yvals_lwr_proliferation_W, 
					yvals_upr=yvals_upr_proliferation_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=lod,yend=lod-dp_mean_W),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_B, 
					yvals_lwr=yvals_lwr_proliferation_B, 
					yvals_upr=yvals_upr_proliferation_B),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=-wp_mean_B,xend=0,y=lod,yend=lod-dp_mean_B),col=specialfill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=yvals_lwr_clearance_W, 
					yvals_upr=yvals_upr_clearance_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=0,xend=wr_mean_W,y=lod-dp_mean_W,yend=lod),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_B, 
					yvals_lwr=yvals_lwr_clearance_B, 
					yvals_upr=yvals_upr_clearance_B),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=0,xend=wr_mean_B,y=lod-dp_mean_B,yend=lod),col=specialfill) + 
			coord_cartesian(ylim=c(40,15), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Ct") + 
			scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + # "log"["10"]" RNA copies/ml"
			theme(text=element_text(size=18))
	} else {
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=(yvals_lwr_proliferation_W), 
					yvals_upr=(yvals_upr_proliferation_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_W)),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_B, 
					yvals_lwr=(yvals_lwr_proliferation_B), 
					yvals_upr=(yvals_upr_proliferation_B)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=-wp_mean_B,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_B)),col=specialfill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=(yvals_lwr_clearance_W), 
					yvals_upr=(yvals_upr_clearance_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=0,xend=wr_mean_W,y=10^convert_Ct_logGEML(lod-dp_mean_W),yend=10^convert_Ct_logGEML(lod)),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_B, 
					yvals_lwr=(yvals_lwr_clearance_B), 
					yvals_upr=(yvals_upr_clearance_B)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=0,xend=wr_mean_B,y=10^convert_Ct_logGEML(lod-dp_mean_B),yend=10^convert_Ct_logGEML(lod)),col=specialfill) + 
			coord_cartesian(ylim=c(10^convert_Ct_logGEML(40),10^convert_Ct_logGEML(15)), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="RNA copies per ml") + 
			scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
			theme(text=element_text(size=18))
	}

	return(out)

	})

}


make_sample_trajectory_special_trap <- function(shared_params_df, global_pars, siglevel=0.95, ge=FALSE, referencefill="blue", specialfill="red"){
	# For asymptomatic:
	with(as.list(global_pars),{

	shared_params_df <- shared_params_df %>% 
		select(-dpmeanW, -dpmeanB, -wpmeanW, -wpmeanB, -wxmeanW, -wxmeanB, -wrmeanW, -wrmeanB) %>%
		rename(dpmeanW=dpmeanW_trans,
		dpmeanB=dpmeanB_trans,
		wpmeanW=wpmeanW_trans,
		wpmeanB=wpmeanB_trans,
		wxmeanW=wxmeanW_trans,
		wxmeanB=wxmeanB_trans,
		wrmeanW=wrmeanW_trans,
		wrmeanB=wrmeanB_trans)

	wp_mean_W <- mean(shared_params_df$wpmeanW)
	wp_lwr_W <- quantile(shared_params_df$wpmeanW,(1-siglevel)/2)
	wp_upr_W <- quantile(shared_params_df$wpmeanW,1-(1-siglevel)/2)

	wx_mean_W <- mean(shared_params_df$wxmeanW)
	wx_lwr_W <- quantile(shared_params_df$wxmeanW,(1-siglevel)/2)
	wx_upr_W <- quantile(shared_params_df$wxmeanW,1-(1-siglevel)/2)

	wr_mean_W <- mean(shared_params_df$wrmeanW)
	wr_lwr_W <- quantile(shared_params_df$wrmeanW,(1-siglevel)/2)
	wr_upr_W <- quantile(shared_params_df$wrmeanW,1-(1-siglevel)/2)

	dp_mean_W <- mean(shared_params_df$dpmeanW)
	dp_lwr_W <- quantile(shared_params_df$dpmeanW,(1-siglevel)/2)
	dp_upr_W <- quantile(shared_params_df$dpmeanW,1-(1-siglevel)/2)

	xvals_proliferation_W <- seq(from=-wp_upr_W, 0, length.out=500)
	# xvals_plateau_W <- seq(from=0, wx_upr_W, length.out=500)
	xvals_clearance_W <- seq(from=0, wx_upr_W+wr_upr_W, length.out=500)

	yvals_upr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),(1-siglevel)/2)))
	# yvals_upr_plateau_W <- rep(dp_upr_W, length(xvals_plateau_W))
	# yvals_lwr_plateau_W <- rep(dp_lwr_W, length(xvals_plateau_W))
	yvals_upr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall_trap(x, shared_params_df$dpmeanW, shared_params_df$wxmeanW, shared_params_df$wrmeanW),1-(1-siglevel)/2)))
	yvals_lwr_clearance_W <- unlist(lapply(xvals_clearance_W, 
	function(x) quantile(sfall_trap(x, shared_params_df$dpmeanW, shared_params_df$wxmeanW, shared_params_df$wrmeanW),(1-siglevel)/2)))

	# For symptomatic:
	wp_mean_B <- mean(shared_params_df$wpmeanB)
	wp_lwr_B <- quantile(shared_params_df$wpmeanB,(1-siglevel)/2)
	wp_upr_B <- quantile(shared_params_df$wpmeanB,1-(1-siglevel)/2)

	wx_mean_B <- mean(shared_params_df$wxmeanB)
	wx_lwr_B <- quantile(shared_params_df$wxmeanB,(1-siglevel)/2)
	wx_upr_B <- quantile(shared_params_df$wxmeanB,1-(1-siglevel)/2)

	wr_mean_B <- mean(shared_params_df$wrmeanB)
	wr_lwr_B <- quantile(shared_params_df$wrmeanB,(1-siglevel)/2)
	wr_upr_B <- quantile(shared_params_df$wrmeanB,1-(1-siglevel)/2)

	dp_mean_B <- mean(shared_params_df$dpmeanB)
	dp_lwr_B <- quantile(shared_params_df$dpmeanB,(1-siglevel)/2)
	dp_upr_B <- quantile(shared_params_df$dpmeanB,1-(1-siglevel)/2)

	xvals_proliferation_B <- seq(from=-wp_upr_B, 0, length.out=500)
	# xvals_plateau_B <- seq(from=0, wx_upr_B, length.out=500)
	xvals_clearance_B <- seq(from=0, wx_upr_B+wr_upr_B, length.out=500)

	yvals_upr_proliferation_B <- unlist(lapply(xvals_proliferation_B, 
	function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),1-(1-siglevel)/2)))
	yvals_lwr_proliferation_B <- unlist(lapply(xvals_proliferation_B, 
	function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),(1-siglevel)/2)))
	# yvals_upr_plateau_B <- rep(dp_upr_B, length(xvals_plateau_B))
	# yvals_lwr_plateau_B <- rep(dp_lwr_B, length(xvals_plateau_B))
	yvals_upr_clearance_B <- unlist(lapply(xvals_clearance_B, 
	function(x) quantile(sfall_trap(x, shared_params_df$dpmeanB, shared_params_df$wxmeanB, shared_params_df$wrmeanB),1-(1-siglevel)/2)))
	yvals_lwr_clearance_B <- unlist(lapply(xvals_clearance_B, 
	function(x) quantile(sfall_trap(x, shared_params_df$dpmeanB, shared_params_df$wxmeanB, shared_params_df$wrmeanB),(1-siglevel)/2)))

	if(ge==FALSE){
		out <- ggplot() + 
			# Proliferation: ---
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=yvals_lwr_proliferation_W, 
					yvals_upr=yvals_upr_proliferation_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=lod,yend=lod-dp_mean_W),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_B, 
					yvals_lwr=yvals_lwr_proliferation_B, 
					yvals_upr=yvals_upr_proliferation_B),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=-wp_mean_B,xend=0,y=lod,yend=lod-dp_mean_B),col=specialfill) + 
			# Plateau: ---
			# geom_ribbon(
			# 	data=tibble(
			# 		xvals=xvals_plateau_W, 
			# 		yvals_lwr=yvals_lwr_plateau_W, 
			# 		yvals_upr=yvals_upr_plateau_W),
			# 	aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=0,xend=wx_mean_W,y=lod-dp_mean_W,yend=lod-dp_mean_W),col=referencefill) + 
			# geom_ribbon(
			# 	data=tibble(
			# 		xvals=xvals_plateau_B, 
			# 		yvals_lwr=yvals_lwr_plateau_B, 
			# 		yvals_upr=yvals_upr_plateau_B),
			# 	aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=0,xend=wx_mean_B,y=lod-dp_mean_B,yend=lod-dp_mean_B),col=specialfill) + 
			# Clearance: ---
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=yvals_lwr_clearance_W, 
					yvals_upr=yvals_upr_clearance_W),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=wx_mean_W,xend=wx_mean_W+wr_mean_W,y=lod-dp_mean_W,yend=lod),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_B, 
					yvals_lwr=yvals_lwr_clearance_B, 
					yvals_upr=yvals_upr_clearance_B),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=wx_mean_B,xend=wx_mean_B+wr_mean_B,y=lod-dp_mean_B,yend=lod),col=specialfill) + 
			# Plot options: ---
			coord_cartesian(ylim=c(40,15), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="Ct") + 
			scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + # "log"["10"]" RNA copies/ml"
			theme(text=element_text(size=18))
	} else {
		out <- ggplot() + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_W, 
					yvals_lwr=(yvals_lwr_proliferation_W), 
					yvals_upr=(yvals_upr_proliferation_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=-wp_mean_W,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_W)),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation_B, 
					yvals_lwr=(yvals_lwr_proliferation_B), 
					yvals_upr=(yvals_upr_proliferation_B)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=-wp_mean_B,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_B)),col=specialfill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_W, 
					yvals_lwr=(yvals_lwr_clearance_W), 
					yvals_upr=(yvals_upr_clearance_W)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=referencefill) + 
			geom_segment(aes(x=0,xend=wr_mean_W,y=10^convert_Ct_logGEML(lod-dp_mean_W),yend=10^convert_Ct_logGEML(lod)),col=referencefill) + 
			geom_ribbon(
				data=tibble(
					xvals=xvals_clearance_B, 
					yvals_lwr=(yvals_lwr_clearance_B), 
					yvals_upr=(yvals_upr_clearance_B)),
				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill=specialfill) + 
			geom_segment(aes(x=0,xend=wr_mean_B,y=10^convert_Ct_logGEML(lod-dp_mean_B),yend=10^convert_Ct_logGEML(lod)),col=specialfill) + 
			coord_cartesian(ylim=c(10^convert_Ct_logGEML(40),10^convert_Ct_logGEML(15)), expand=FALSE) + 
			theme_minimal() + 
			labs(x="Days from peak", y="RNA copies per ml") + 
			scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
			theme(text=element_text(size=18))
	}

	return(out)

	})

}

# make_sample_trajectory_single(shared_params_df, global_pars, siglevel=0.95)
make_sample_trajectory_single <- function(shared_params_df, global_pars, siglevel=0.95, linecol="black", special=TRUE){
	# For asymptomatic:
	with(as.list(global_pars),{

	shared_params_df <- shared_params_df %>% 
		select(-dpmeanW, -dpmeanB, -wpmeanW, -wpmeanB, -wrmeanW, -wrmeanB) %>%
		rename(dpmeanW=dpmeanW_trans,
		dpmeanB=dpmeanB_trans,
		wpmeanW=wpmeanW_trans,
		wpmeanB=wpmeanB_trans,
		wrmeanW=wrmeanW_trans,
		wrmeanB=wrmeanB_trans)

	if(special==FALSE){
		wp_mean <- mean(shared_params_df$wpmeanW)
		wp_lwr <- quantile(shared_params_df$wpmeanW,(1-siglevel)/2)
		wp_upr <- quantile(shared_params_df$wpmeanW,1-(1-siglevel)/2)

		wr_mean <- mean(shared_params_df$wrmeanW)
		wr_lwr <- quantile(shared_params_df$wrmeanW,(1-siglevel)/2)
		wr_upr <- quantile(shared_params_df$wrmeanW,1-(1-siglevel)/2)

		dp_mean <- mean(shared_params_df$dpmeanW)
		dp_lwr <- quantile(shared_params_df$dpmeanW,(1-siglevel)/2)
		dp_upr <- quantile(shared_params_df$dpmeanW,1-(1-siglevel)/2)

		xvals_proliferation <- seq(from=-wp_upr, 0, length.out=500)
		xvals_clearance <- seq(from=0, wr_upr, length.out=500)

		yvals_upr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),1-(1-siglevel)/2)))
		yvals_lwr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),(1-siglevel)/2)))
		yvals_upr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),1-(1-siglevel)/2)))
		yvals_lwr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),(1-siglevel)/2)))
	} else {
		wp_mean <- mean(shared_params_df$wpmeanB)
		wp_lwr <- quantile(shared_params_df$wpmeanB,(1-siglevel)/2)
		wp_upr <- quantile(shared_params_df$wpmeanB,1-(1-siglevel)/2)

		wr_mean <- mean(shared_params_df$wrmeanB)
		wr_lwr <- quantile(shared_params_df$wrmeanB,(1-siglevel)/2)
		wr_upr <- quantile(shared_params_df$wrmeanB,1-(1-siglevel)/2)

		dp_mean <- mean(shared_params_df$dpmeanB)
		dp_lwr <- quantile(shared_params_df$dpmeanB,(1-siglevel)/2)
		dp_upr <- quantile(shared_params_df$dpmeanB,1-(1-siglevel)/2)

		xvals_proliferation <- seq(from=-wp_upr, 0, length.out=500)
		xvals_clearance <- seq(from=0, wr_upr, length.out=500)

		yvals_upr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),1-(1-siglevel)/2)))
		yvals_lwr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),(1-siglevel)/2)))
		yvals_upr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),1-(1-siglevel)/2)))
		yvals_lwr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),(1-siglevel)/2)))
	}

	out <- ggplot() + 
		geom_line(data=tibble(
				xvals=xvals_proliferation, 
				yvals=yvals_lwr_proliferation),
			aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed")+
		geom_line(data=tibble(
				xvals=xvals_proliferation, 
				yvals=yvals_upr_proliferation),
			aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed")+
		geom_segment(aes(x=-wp_mean,xend=0,y=lod,yend=lod-dp_mean),col=linecol) + 
		geom_line(
			data=tibble(
				xvals=xvals_clearance, 
				yvals=yvals_lwr_clearance),
			aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed") + 
		geom_line(
			data=tibble(
				xvals=xvals_clearance, 
				yvals=yvals_upr_clearance),
			aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed") + 
		geom_segment(aes(x=0,xend=wr_mean,y=lod-dp_mean,yend=lod),col=linecol) + 
		coord_cartesian(ylim=c(40,15), expand=FALSE) + 
		theme_minimal() + 
		labs(x="Days from peak", y="Ct") + 
		scale_y_reverse(sec.axis=sec_axis(~convert_Ct_logGEML(.), name=expression(log[10]~RNA~copies/ml))) + # "log"["10"]" RNA copies/ml"
		theme(text=element_text(size=18))

	return(out)

	})

}

# temp <- make_sample_trajectory_single_list(shared_params_df, global_pars, linecol="black", alphalevel=0.2, siglevel=0.95)
# fig_nv_raw_points + temp
make_sample_trajectory_single_list <- function(shared_params_df, global_pars, siglevel=0.95, alphalevel=0.4, linecol="black", special=TRUE){
	# For asymptomatic:
	with(as.list(global_pars),{

	shared_params_df <- shared_params_df %>% 
		select(-dpmeanW, -dpmeanB, -wpmeanW, -wpmeanB, -wrmeanW, -wrmeanB) %>%
		rename(dpmeanW=dpmeanW_trans,
		dpmeanB=dpmeanB_trans,
		wpmeanW=wpmeanW_trans,
		wpmeanB=wpmeanB_trans,
		wrmeanW=wrmeanW_trans,
		wrmeanB=wrmeanB_trans)

	if(special==FALSE){
		wp_mean <- mean(shared_params_df$wpmeanW)
		wp_lwr <- quantile(shared_params_df$wpmeanW,(1-siglevel)/2)
		wp_upr <- quantile(shared_params_df$wpmeanW,1-(1-siglevel)/2)

		wr_mean <- mean(shared_params_df$wrmeanW)
		wr_lwr <- quantile(shared_params_df$wrmeanW,(1-siglevel)/2)
		wr_upr <- quantile(shared_params_df$wrmeanW,1-(1-siglevel)/2)

		dp_mean <- mean(shared_params_df$dpmeanW)
		dp_lwr <- quantile(shared_params_df$dpmeanW,(1-siglevel)/2)
		dp_upr <- quantile(shared_params_df$dpmeanW,1-(1-siglevel)/2)

		xvals_proliferation <- seq(from=-wp_upr, 0, length.out=500)
		xvals_clearance <- seq(from=0, wr_upr, length.out=500)

		yvals_upr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),1-(1-siglevel)/2)))
		yvals_lwr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),(1-siglevel)/2)))
		yvals_upr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),1-(1-siglevel)/2)))
		yvals_lwr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),(1-siglevel)/2)))
	} else {
		wp_mean <- mean(shared_params_df$wpmeanB)
		wp_lwr <- quantile(shared_params_df$wpmeanB,(1-siglevel)/2)
		wp_upr <- quantile(shared_params_df$wpmeanB,1-(1-siglevel)/2)

		wr_mean <- mean(shared_params_df$wrmeanB)
		wr_lwr <- quantile(shared_params_df$wrmeanB,(1-siglevel)/2)
		wr_upr <- quantile(shared_params_df$wrmeanB,1-(1-siglevel)/2)

		dp_mean <- mean(shared_params_df$dpmeanB)
		dp_lwr <- quantile(shared_params_df$dpmeanB,(1-siglevel)/2)
		dp_upr <- quantile(shared_params_df$dpmeanB,1-(1-siglevel)/2)

		xvals_proliferation <- seq(from=-wp_upr, 0, length.out=500)
		xvals_clearance <- seq(from=0, wr_upr, length.out=500)

		yvals_upr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),1-(1-siglevel)/2)))
		yvals_lwr_proliferation <- unlist(lapply(xvals_proliferation, 
		function(x) quantile(srise(x, shared_params_df$dpmeanB, shared_params_df$wpmeanB),(1-siglevel)/2)))
		yvals_upr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),1-(1-siglevel)/2)))
		yvals_lwr_clearance <- unlist(lapply(xvals_clearance, 
		function(x) quantile(sfall(x, shared_params_df$dpmeanB, shared_params_df$wrmeanB),(1-siglevel)/2)))
	}

	out <- list(
		# geom_line(data=tibble(
		# 		xvals=xvals_proliferation, 
		# 		yvals=yvals_lwr_proliferation),
		# 	aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed"),
		# geom_line(data=tibble(
		# 		xvals=xvals_proliferation, 
		# 		yvals=yvals_upr_proliferation),
		# 	aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed"),
		geom_ribbon(
				data=tibble(
					xvals=xvals_proliferation, 
					yvals_lwr=yvals_lwr_proliferation, 
					yvals_upr=yvals_upr_proliferation),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=alphalevel, fill=linecol),
		geom_segment(aes(x=-wp_mean,xend=0,y=lod,yend=lod-dp_mean),col=linecol) , 
		# geom_line(
		# 	data=tibble(
		# 		xvals=xvals_clearance, 
		# 		yvals=yvals_lwr_clearance),
		# 	aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed") , 
		# geom_line(
		# 	data=tibble(
		# 		xvals=xvals_clearance, 
		# 		yvals=yvals_upr_clearance),
		# 	aes(x=xvals, y=lod-yvals), alpha=0.6, col=linecol, linetype="dashed") , 
		geom_ribbon(
				data=tibble(
					xvals=xvals_clearance, 
					yvals_lwr=yvals_lwr_clearance, 
					yvals_upr=yvals_upr_clearance),
				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=alphalevel, fill=linecol),
		geom_segment(aes(x=0,xend=wr_mean,y=lod-dp_mean,yend=lod),col=linecol),
		coord_cartesian(ylim=c(40,11), expand=FALSE)
		)

	return(out)

	})

}


# make_sample_trajectory_wtonly <- function(shared_params_df, global_pars, siglevel=0.9, ge=FALSE){
# 	# For asymptomatic:
# 	with(as.list(global_pars),{

# 	wp_mean_W <- mean(shared_params_df$wpmeanW)
# 	wp_lwr_W <- quantile(shared_params_df$wpmeanW,(1-siglevel)/2)
# 	wp_upr_W <- quantile(shared_params_df$wpmeanW,1-(1-siglevel)/2)

# 	wr_mean_W <- mean(shared_params_df$wrmeanW)
# 	wr_lwr_W <- quantile(shared_params_df$wrmeanW,(1-siglevel)/2)
# 	wr_upr_W <- quantile(shared_params_df$wrmeanW,1-(1-siglevel)/2)

# 	dp_mean_W <- mean(shared_params_df$dpmeanW)
# 	dp_lwr_W <- quantile(shared_params_df$dpmeanW,(1-siglevel)/2)
# 	dp_upr_W <- quantile(shared_params_df$dpmeanW,1-(1-siglevel)/2)

# 	xvals_proliferation_W <- seq(from=-wp_upr_W, 0, length.out=500)
# 	xvals_clearance_W <- seq(from=0, wr_upr_W, length.out=500)

# 	yvals_upr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
# 	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),1-(1-siglevel)/2)))
# 	yvals_lwr_proliferation_W <- unlist(lapply(xvals_proliferation_W, 
# 	function(x) quantile(srise(x, shared_params_df$dpmeanW, shared_params_df$wpmeanW),(1-siglevel)/2)))
# 	yvals_upr_clearance_W <- unlist(lapply(xvals_clearance_W, 
# 	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),1-(1-siglevel)/2)))
# 	yvals_lwr_clearance_W <- unlist(lapply(xvals_clearance_W, 
# 	function(x) quantile(sfall(x, shared_params_df$dpmeanW, shared_params_df$wrmeanW),(1-siglevel)/2)))


# 	if(ge==FALSE){
# 		out <- ggplot() + 
# 			geom_ribbon(
# 				data=tibble(
# 					xvals=xvals_proliferation_W, 
# 					yvals_lwr=yvals_lwr_proliferation_W, 
# 					yvals_upr=yvals_upr_proliferation_W),
# 				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="blue") + 
# 			geom_segment(aes(x=-wp_mean_W,xend=0,y=lod,yend=lod-dp_mean_W),col="blue") + 
# 			geom_ribbon(
# 				data=tibble(
# 					xvals=xvals_clearance_W, 
# 					yvals_lwr=yvals_lwr_clearance_W, 
# 					yvals_upr=yvals_upr_clearance_W),
# 				aes(x=xvals, ymin=lod-yvals_lwr, ymax=lod-yvals_upr), alpha=0.2, fill="blue") + 
# 			geom_segment(aes(x=0,xend=wr_mean_W,y=lod-dp_mean_W,yend=lod),col="blue") + 
# 			coord_cartesian(ylim=c(40,15), expand=FALSE) + 
# 			theme_minimal() + 
# 			labs(x="Days from peak", y="Ct") + 
# 			scale_y_reverse() + 
# 			theme(text=element_text(size=18))
# 	} else {
# 		out <- ggplot() + 
# 			geom_ribbon(
# 				data=tibble(
# 					xvals=xvals_proliferation_W, 
# 					yvals_lwr=(yvals_lwr_proliferation_W), 
# 					yvals_upr=(yvals_upr_proliferation_W)),
# 				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="blue") + 
# 			geom_segment(aes(x=-wp_mean_W,xend=0,y=10^convert_Ct_logGEML(lod),yend=10^convert_Ct_logGEML(lod-dp_mean_W)),col="blue") + 
# 			geom_ribbon(
# 				data=tibble(
# 					xvals=xvals_clearance_W, 
# 					yvals_lwr=(yvals_lwr_clearance_W), 
# 					yvals_upr=(yvals_upr_clearance_W)),
# 				aes(x=xvals, ymin=10^convert_Ct_logGEML(lod-yvals_lwr), ymax=10^convert_Ct_logGEML(lod-yvals_upr)), alpha=0.2, fill="blue") + 
# 			geom_segment(aes(x=0,xend=wr_mean_W,y=10^convert_Ct_logGEML(lod-dp_mean_W),yend=10^convert_Ct_logGEML(lod)),col="blue") + 
# 			coord_cartesian(ylim=c(10^convert_Ct_logGEML(40),10^convert_Ct_logGEML(15)), expand=FALSE) + 
# 			theme_minimal() + 
# 			labs(x="Days from peak", y="RNA copies per ml") + 
# 			scale_y_continuous(trans='log10', labels = trans_format("log10", math_format(10^.x))) + 
# 			theme(text=element_text(size=18))
# 	}

# 	return(out)

# 	})

# }

make_normal_prior_df <- function(mean, sd, pmin, pmax, step){
	xmin <- qnorm(pmin, mean=mean, sd=sd)
	xmax <- qnorm(pmax, mean=mean, sd=sd)
	xvals <- seq(from=xmin, to=xmax, by=step)
	out <- as_tibble(data.frame(
		x=xvals, 
		density=dnorm(xvals, mean=mean, sd=sd)))
	return(out)
}


truncnormmean <- function(mu, sigma, T0, T1){
	alpha <- (T0 - mu)/sigma
	beta <- (T1 - mu)/sigma
	out <- mu + 
		(dnorm(alpha, 0, 1) - dnorm(beta, 0, 1))/
		(pnorm(beta, 0, 1) - pnorm(alpha, 0, 1))*sigma
	return(out)
}