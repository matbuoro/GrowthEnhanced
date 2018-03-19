# POSTERIOR CHECK
library("rstanarm")
library(bayesplot)

  # Load results (mcmc)
  load(paste0("results/",myvar,"_",data.name,".Rdata")) ## load .Rdata from results/ folder 
  mcmc <- as.array(best$stanfit) ## extract mcmc samples
  
  # Precictive values
  y <- best$y # extract data
  yrep <- posterior_predict(best) # simulate posterior prediction
  group <- best$data$Treatment # group (e.g. SHAM vs GH)
  
  # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
  plot1 <- ppc_stat(y, yrep)
  plot2 <- ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
  plot3 <- pp_check(y, yrep, ppc_dens_overlay)
  # # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
  plot4 <- ppc_scatter_avg(y, yrep)
  # ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)
  


