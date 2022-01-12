# =============================================================================
# Import 
# =============================================================================

library(tidyverse) 
# Assuming that we've already generated a ct_fit object: 

# =============================================================================
# Generate Rhat diagnostics:  
# =============================================================================

rhatvals <- as.data.frame((summary(ct_fit)$summary)[,"Rhat"])
rhatvals <- cbind(rownames(rhatvals), data.frame(rhatvals, row.names=NULL))
names(rhatvals) <- c("param","Rhat")
rhatvals <- as_tibble(rhatvals) 

rhatwarnings <- rhatvals %>% 
	filter(Rhat>1.1) %>% 
	summarise(rhatabove1.1=n())

# =============================================================================
# Generate divergences summary: 
# =============================================================================

ndivergent <- Reduce(rbind, get_sampler_params(ct_fit, inc_warmup=FALSE)) %>% 
	as_tibble() %>% 
	summarise(divergent=sum(divergent__))

# =============================================================================
# Combine warnings: 
# =============================================================================

warningtab <- cbind(rhatwarnings, ndivergent)