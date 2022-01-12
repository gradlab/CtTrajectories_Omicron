# =============================================================================
# Fit model for all parameter sets
# =============================================================================

library(tidyverse) 
library(scales)
library(lubridate)

source('code/utils.R')
source('code/utils_private.R')
source("code/set_global_pars.R")
load("data/ct_dat_refined.RData")
source("code/set_run_pars.R")

for(run_pars_index in 6:6){ # length(run_pars_list)

	run_pars <- run_pars_list[[run_pars_index]]

	source("code/fit_posteriors_preamble.R")
	source("code/fit_posteriors.R")
	source("code/extract_parameters.R")

	source("code/make_figures.R")
	source("code/summarise_dists.R")

	savedir <- paste0("figures/run_pars_",run_pars_index,"/")
	savedir_output <- paste0("output/run_pars_",run_pars_index,"/")

	source("code/save_figures.R")
	write_csv(dist_summary, file=paste0(savedir,"dist_summary.csv"))
	save(ct_fit, file=paste0(extdrive,savedir_output,"ct_fit.RData"))

	source('code/diagnose_fit.R')
	write_csv(warningtab, file=paste0(savedir,"warningtab.csv"))

	print(paste0("Finished parameter set ",run_pars_index,": ",Sys.time()))

}
