rm(list=ls())   # Clear memory

#--------------- PACKAGES ---------------#
library(readxl)
# see https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library("rstanarm")
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
#library(betareg)
#library(bayesplot)
#library(gridExtra)

SEED <- 1234
set.seed(SEED)

# # Install from Github (development version)
# if (!require(devtools)) {
#   install.packages("devtools")
#   library(devtools)
# }
# install_github("stan-dev/rstanarm", args = "--preclean", build_vignettes = FALSE)

#--------------- FUNCTIONS ---------------#
invlogit<-function(x) {1/(1+exp(-(x)))} # inverse logit
logit<-function(x) {log(x/(1-x))} # logisitic transformation
gf <- function(x){if(median(x)>=0){mean(x<=0)}else{mean(x>0)}} # function to calculate the proportion of posterior with same sign as mean


#--------------- DIRECTORY --------------#
setwd("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/")
#setwd("~/Documents/mbuoro/Ims/")

#-----------------ANALYSIS---------------#
CHAINS  = 3 # number of chains
CORES = CHAINS # number of core
ITER = 10000 # number of iterations
WARM = 2000 # warmup
ndelta = 0.99




###______________________________________________###
###_________________DATA_________________________###
###______________________________________________###


###_________________EXCRETION____________________###
# Load dataset:
excretion <- read_excel("data/Stream_Channels_2016_Mat.xlsx", sheet = "Population excretion")
excretion <- as.data.frame(excretion)
excretion$Sum_SRP <- as.numeric(excretion$Sum_SRP)
excretion$Sum_NH4 <- as.numeric(excretion$Sum_NH4)
excretion$Treatment <- as.factor(excretion$Treatment)
excretion <- within(excretion, Treatment <- relevel(Treatment, ref = 3)) # SHAM as reference
#View(df)


###____________ Sum_SRP    _____________ ###
fit1 <- stan_glmer(
  log(Sum_SRP) ~ Treatment +(1|Block), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = excretion,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Sum_SRP"
  Sum_SRP <- fit1
  ### SAVE
  save(best,file=paste0("results/IMS2016_",model.name,".Rdata"))
}




###____________ Sum_NH4    _____________ ###
fit1 <- stan_glmer(
  log(Sum_NH4) ~ Treatment +(1|Block), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = excretion,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)



best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Sum_NH4"
  Sum_SRP <- best
  ### SAVE
  save(best,file=paste0("results/IMS2016_",model.name,".Rdata"))
}






###______________________________________________###
###_________________NUTRIENT ____________________###
###______________________________________________###


nutrient <- subset(nutrient, Time != "T0" & Position != "Up" )



###____________ NH4    _____________ ###
fit1 <- stan_glmer(
  log(NH4) ~Treatment+Time+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = nutrient,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  log(NH4) ~Treatment*Time +(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = nutrient,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


# fit3 <- stan_glmer(
#   log(NH4) ~Treatment + Time + Position+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = nutrient,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit4 <- stan_glmer(
#   log(NH4) ~Treatment*Position+ Time +(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = nutrient,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
# loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position +(1|Block/Section)
# loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position+Time +(1|Block/Section)
loo <- compare(loo1, loo2)#, loo3, loo4) 

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)#, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "NH4"
  NH4 <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}

mcmc <- as.array(best$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat("Effect of GH: ", median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]\n")
gf(TreatmentGH)

y <- best$y
yrep <- posterior_predict(best)
group <- best$data$Treatment


#par(mfrow=c(2,2))
#par(mar=c(3.5, 4.5, 4.5, 1.5))
# # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
plot1 <- ppc_stat(y, yrep)
plot2 <- ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
#color_scheme_set("teal")
#ppc_stat_grouped(y, yrep, group)
plot3 <- pp_check(y, yrep, ppc_dens_overlay)
# pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
# # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
plot4 <- ppc_scatter_avg(y, yrep)
# ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)

grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)




###____________ PO4    _____________ ###
fit1 <- stan_glmer(
  log(PO4) ~Treatment+Time+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = nutrient,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  log(PO4) ~Treatment*Time +(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = nutrient,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


# fit3 <- stan_glmer(
#   log(PO4+1) ~Treatment + Time + Position+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = nutrient,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit4 <- stan_glmer(
#   log(PO4+1) ~Treatment*Position+ Time +(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = nutrient,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
# loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position +(1|Block/Section)
# loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position+Time +(1|Block/Section)
loo <- compare(loo1, loo2, loo3, loo4) 


# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)#, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "PO4"
  PO4 <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}


mcmc <- as.array(best$stanfit)

TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
TreatmentGHLD <- as.vector(mcmc[,1:CORES,"TreatmentGH_LD"])
TreatmentNF <- as.vector(mcmc[,1:CORES,"TreatmentNF"])

x <- TreatmentGHLD
cat("Effect of GH: ", median(x),"["
    ,quantile(x, probs=0.025),";"
    ,quantile(x, probs=0.975),"]\n")
gf(x)

y <- best$y
yrep <- posterior_predict(best)
group <- best$data$Treatment


#par(mfrow=c(2,2))
#par(mar=c(3.5, 4.5, 4.5, 1.5))
# # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
plot1 <- ppc_stat(y, yrep)
plot2 <- ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
#color_scheme_set("teal")
#ppc_stat_grouped(y, yrep, group)
plot3 <- pp_check(y, yrep, ppc_dens_overlay)
# pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
# # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
plot4 <- ppc_scatter_avg(y, yrep)
# ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)

grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)





###______________________________________________###
###_________________NUTRIENT FLUX____________________###
###______________________________________________###


# Load dataset:
flux <- read_excel("data/Flux_Nutrients_to_Mat.xlsx")
flux <- as.data.frame(flux)
flux$Flux_PO4 <- as.numeric(flux$Flux_PO4)
flux$Flux_NH4 <- as.numeric(flux$Flux_NH4)
flux$Flux_NNO3 <- as.numeric(flux$Flux_NNO3)
flux$Flux_N_total <- as.numeric(flux$Flux_N_total)
flux$Time <- as.factor(flux$Time)
flux$Block <- as.factor(flux$Block)
flux$Section <- as.factor(flux$Stream_channel)
flux$Treatment[flux$Treatment=="GHA"] <- "GH"
flux$Treatment[flux$Treatment=="GHB"] <- "GH_LD"
flux$Treatment <- as.factor(flux$Treatment)
flux <- within(flux, Treatment <- relevel(Treatment, ref = 4)) # SHAM as reference
# #View(df)
# 
# # up <- subset(nutrient,Position == "Up")
# # down <- subset(nutrient,Position == "Down")
# # # DELTA = UP - DOWN
# # deltas.NH4 <- up$NH4 - down$NH4
# # deltas.PO4 <- up$PO4 - down$PO4
# # ratio.NH4 <- (up$NH4/down$NH4) 
# # up$PO4 <- up$PO4 + .1
# # down$PO4 <- down$PO4 + .1
# # ratio.PO4 <- (up$PO4/down$PO4 )
# # delta <- cbind(up,deltas.NH4=deltas.NH4,deltas.PO4=deltas.PO4, ratio.NH4=ratio.NH4, ratio.PO4=ratio.PO4)
# # delta <- subset(delta, Time != "T0")
# 
###____________ DELTA NH4    _____________ ###
fit1 <- stan_glmer(
  Flux_NH4 ~Treatment*Time+(1|Block/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = flux,
  prior = student_t(df = 7),
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Flux_NH4 ~Treatment+Time+(1|Block/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = flux,
  prior = student_t(df = 7),
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
loo <- compare(loo1, loo2)
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher.
# If the difference is be positive then the second model is preferred.

#Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred.
#That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)



## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "FluxNH4"
  FluxNH4 <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}




###____________ DELTA PO4    _____________ ###
fit1 <- stan_glmer(
  Flux_PO4 ~Treatment*Time+(1|Block/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = flux,
  prior = student_t(df = 7),
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Flux_PO4 ~Treatment+Time+(1|Block/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = flux,
  prior = student_t(df = 7),
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
loo <- compare(loo1, loo2)
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher.
# If the difference is be positive then the second model is preferred.

#Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred.
#That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)


## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "FluxPO4"
  FluxPO4 <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}
# 
# 
# 
# ###____________ DELTA Flux_NNO3    _____________ ###
# fit1 <- stan_glmer(
#   Flux_NNO3 ~Treatment+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = subset(flux, Time=="T2"),
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# best <-fit1
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Flux_NNO3"
#   Flux_NNO3 <- best
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
# }
# 
# 
# ###____________ DELTA Flux_N_total    _____________ ###
# fit1 <- stan_glmer(
#   Flux_N_total ~Treatment+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = subset(flux, Time=="T2"),
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# best <-fit1
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Flux_N_total"
#   Flux_N_total <- best
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
# }
# 




###______________________________________________###
###______________Filamentous ____________________###
###______________________________________________###


# Load dataset:
filamentous <- read_excel("data/Stream_Channels_2016_Mat.xlsx", sheet = "Filamentous")
filamentous <- as.data.frame(filamentous)
filamentous$Fila <- as.numeric(filamentous$Fila)
filamentous$Block <- as.factor(filamentous$Block)
filamentous$Section <- as.factor(filamentous$Section)
filamentous$Position <- as.factor(filamentous$Position)
filamentous$Treatment <- as.factor(filamentous$Treatment)
filamentous <- within(filamentous, Treatment <- relevel(Treatment, ref = 4)) # SHAM as reference
#View(df)

# data transformed
filamentous$Fila[filamentous$Fila==1] <- .99 # must be <1
filamentous$logit.Fila <- logit(filamentous$Fila)


###____________ Fila    _____________ ###
## /!\ Does not work with (1|Block/Section)
# Fila <- stan_betareg(
#   Fila ~ Treatment+Position,#+(1|Block/Section),
#   link = "logit", link.phi = "log", 
#   algorithm = "sampling",
#   data = filamentous,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),QR = TRUE,
# )


fit1 <- stan_glmer(
  logit.Fila ~ Treatment*Position+(1|Block/Section), 
  #family = binomial(link=logit),
  na.action = "na.omit",
  data = filamentous,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)



fit2 <- stan_glmer(
  logit.Fila ~Treatment + Position+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = filamentous,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position +(1|Block/Section)
loo <- compare(loo1, loo2) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

#Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
#That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)


## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Fila"
  Fila <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}




###______________________________________________###
###______________Others_______ __________________###
###______________________________________________###

# Load dataset:
sp <- read_excel("data/Stream_Channels_2016_Mat.xlsx", sheet = "Autres")
sp <- as.data.frame(sp)
sp$K_CM <- as.numeric(sp$K_CM)
sp$K_FM <- as.numeric(sp$K_FM)
sp$Cyanobacteria <- as.numeric(sp$Cyanobacteria)
sp$Green_algae <- as.numeric(sp$Green_algae)
sp$Diatoms <- as.numeric(sp$Diatoms)
sp$All_Chloro <- as.numeric(sp$All_Chloro)
sp$Nb_Chironomidae <- as.numeric(sp$Nb_Chironomidae)
sp$Time <- as.factor(sp$Time)
sp$Block <- as.factor(sp$Block)
sp$Section <- as.factor(sp$Section)
sp$Position <- as.factor(sp$Position)
sp$Position_Bis <- as.factor(sp$Position_Bis)
sp$Treatment <- as.factor(sp$Treatment)

sp <- within(sp, Treatment <- relevel(Treatment, ref = 4)) # SHAM as reference
#View(df)



###____________ K_CM    _____________ ###
fit1 <- stan_glmer(
  log(K_CM) ~ Treatment*Time*Position_Bis+(1|Block/Section),
  #family = binomial(link=logit),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  log(K_CM) ~Treatment*Time + Position_Bis+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  log(K_CM) ~Treatment + Time + Position_Bis+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  log(K_CM) ~Treatment*Position_Bis+ Time +(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position_Bis +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position_Bis +(1|Block/Section)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position_Bis +(1|Block/Section)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position_Bis+Time +(1|Block/Section)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

#Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
#That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "K_CM"
  K_CM <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}





###____________ K_FM    _____________ ###
fit1 <- stan_glmer(
  log(K_FM) ~ Treatment*Time*Position_Bis+(1|Block/Section),
  #family = binomial(link=logit),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  log(K_FM) ~Treatment*Time + Position_Bis+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  log(K_FM) ~Treatment + Time + Position_Bis+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  log(K_FM) ~Treatment*Position_Bis+ Time +(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position_Bis +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position_Bis +(1|Block/Section)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position_Bis +(1|Block/Section)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position_Bis+Time +(1|Block/Section)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

#Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
#That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "K_FM"
  K_FM <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}







# ###____________ Cyanobacteria    _____________ ###
# fit1 <- stan_glmer(
#   Cyanobacteria ~ Treatment*Time*Position_Bis+(1|Block/Section),
#   #family = binomial(link=logit),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit2 <- stan_glmer(
#   Cyanobacteria ~Treatment*Time + Position_Bis+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# fit3 <- stan_glmer(
#   Cyanobacteria ~Treatment + Time + Position_Bis+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit4 <- stan_glmer(
#   Cyanobacteria ~Treatment*Position_Bis+ Time +(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position_Bis +(1|Block/Section)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position_Bis +(1|Block/Section)
# loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position_Bis +(1|Block/Section)
# loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position_Bis+Time +(1|Block/Section)
# loo <- compare(loo1, loo2, loo3, loo4) 
# # The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# # If the difference is be positive then the second model is preferred.
# 
# #Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
# #That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Cyanobacteria"
#   Cyanobacteria <- best
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
# }





# ###____________ Green_algae    _____________ ###
# fit1 <- stan_glmer(
#   log(Green_algae +1) ~ Treatment*Time*Position_Bis+(1|Block/Section),
#   #family = binomial(link=logit),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit2 <- stan_glmer(
#   log(Green_algae +1) ~Treatment*Time + Position_Bis+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# fit3 <- stan_glmer(
#   log(Green_algae +1) ~Treatment + Time + Position_Bis+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit4 <- stan_glmer(
#   log(Green_algae +1) ~Treatment*Position_Bis+ Time +(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position_Bis +(1|Block/Section)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position_Bis +(1|Block/Section)
# loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position_Bis +(1|Block/Section)
# loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position_Bis+Time +(1|Block/Section)
# loo <- compare(loo1, loo2, loo3, loo4) 
# # The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# # If the difference is be positive then the second model is preferred.
# 
# #Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
# #That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# 
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Green_algae"
#   Green_algae <- best
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
# }






# ###____________ Diatoms    _____________ ###
# fit1 <- stan_glmer(
#   log(Diatoms + 1) ~ Treatment*Time*Position_Bis+(1|Block/Section),
#   #family = binomial(link=logit),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit2 <- stan_glmer(
#   log(Diatoms +1) ~Treatment*Time + Position_Bis+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# fit3 <- stan_glmer(
#   log(Diatoms +1) ~Treatment + Time + Position_Bis+(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit4 <- stan_glmer(
#   log(Diatoms +1) ~Treatment*Position_Bis+ Time +(1|Block/Section), 
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = sp,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position_Bis +(1|Block/Section)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position_Bis +(1|Block/Section)
# loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position_Bis +(1|Block/Section)
# loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position_Bis+Time +(1|Block/Section)
# loo <- compare(loo1, loo2, loo3, loo4) 
# # The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# # If the difference is be positive then the second model is preferred.
# 
# #Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
# #That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Diatoms"
#   Diatoms <- best
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
# }







###____________ All_Chloro    _____________ ###
fit1 <- stan_glmer(
  All_Chloro ~ Treatment*Time*Position_Bis+(1|Block/Section),
  #family = binomial(link=logit),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  All_Chloro ~Treatment*Time + Position_Bis+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  All_Chloro ~Treatment + Time + Position_Bis+(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  All_Chloro ~Treatment*Position_Bis+ Time +(1|Block/Section), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position_Bis +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position_Bis +(1|Block/Section)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position_Bis +(1|Block/Section)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position_Bis+Time +(1|Block/Section)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

#Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
#That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "All_Chloro"
  All_Chloro <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}




###____________ Nb_Chironomidae    _____________ ###
fit1 <- stan_glmer(
  Nb_Chironomidae ~ Treatment*Time*Position_Bis+(1|Block/Section),
  family = poisson(link=log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Nb_Chironomidae ~Treatment*Time + Position_Bis+(1|Block/Section), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  Nb_Chironomidae ~Treatment + Time + Position_Bis+(1|Block/Section), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  Nb_Chironomidae ~Treatment*Position_Bis+ Time +(1|Block/Section), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = sp,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position_Bis +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position_Bis +(1|Block/Section)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position_Bis +(1|Block/Section)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position_Bis+Time +(1|Block/Section)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

#Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
#That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Nb_Chironomidae"
  Nb_Chironomidae <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
}



### SAVE ALL
#save(Sum_SRP, Sum_NH4, NH4, PO4, FluxNH4, FluxPO4, FluxNH4, FluxPO4, Fila, K_CM, K_FM ,Cyanobacteria,Green_algae,Diatoms,All_Chloro, Nb_Chironomidae, file=paste0("results/IMS2016.Rdata"))






###______________________________________________###
###______________OUTPUTS_______ __________________###
###______________________________________________###

model.name <- c("Sum_SRP"
                , "Sum_NH4"
                #, "NH4", "PO4"
                , "FluxNH4"
                , "FluxPO4"
                , "Fila"
                , "K_CM"
                , "K_FM"
                #,"Cyanobacteria"
                #,"Green_algae","Diatoms"
                ,"All_Chloro"
                , "Nb_Chironomidae"
                )



## SUMMARY
sink(paste0("results/IMS2016_summary.txt"))
for (myvar in model.name) {

## LOAD MODEL
best = loo = mcmc = NULL
load(paste0("results/IMS2016_",myvar,".Rdata"))
#best <- get(myvar)
mcmc <- as.array(best$stanfit)

## POSTERIOR CHECK
y <- best$y
yrep <- posterior_predict(best)
group <- best$data$Treatment

cat("\n__________",myvar,"__________\n")

if(myvar != "Sum_SRP" & myvar != "Sum_NH4") {  
cat("\n Model selection :\n")
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.
# Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
# That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.

cat( "MODEL 1: Treatment*Position*Time + (1|Block/Section) \n")
cat( "MODEL 2: Treatment*Time+Position + (1|Block/Section) \n")
cat( "MODEL 3: Treatment+Time+Position + (1|Block/Section) \n")
cat( "MODEL 4: Treatment*Position+Time + (1|Block/Section) \n")
cat("\n Loo :\n")
print(loo)
cat( "The difference in ELPD (model1 - model2) will be negative if the expected out-of-sample predictive accuracy of the first model is higher (as AIC, lower is better). If the difference is be positive then the second model is preferred.\n")
# Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
# That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.
}

#launch_shinystan(model)
cat("\nBest model:\n")
print(best$formula)
print(best$family)


#cat("\n____________________\n")
cat("\n Posterior distributions for the effect of Treatments:\n")
mcmc <- as.array(best$stanfit)

if(myvar == "Sum_SRP" || myvar == "Sum_NH4") {  
  SHAM <- mcmc[,1:CORES,"(Intercept)"]
  TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
  # Difference between treatments
  GHvsGHLD <- TreatmentGH - TreatmentGHLD;

  # posterior distributions
  deltas <- cbind(GH = as.vector(TreatmentGH)
                ,GHLD = as.vector(TreatmentGHLD)
                ,GHvsGHLD = as.vector(GHvsGHLD)
  )
}
if(myvar != "Sum_SRP" & myvar != "Sum_NH4") {  
  SHAM <- mcmc[,1:CORES,"(Intercept)"]
  TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
  TreatmentNF <- (mcmc[,1:CORES,"TreatmentNF"])
  # Difference between treatments
  GHvsGHLD <- TreatmentGH - TreatmentGHLD;
  GHvsNF <- TreatmentGH - TreatmentNF; 
  GHLDvsNF <- TreatmentGHLD - TreatmentNF;
  
  if( paste0(best$formula)[3] == "Treatment * Time + Position_Bis + (1 | Block/Section)"){
  TreatmentGH2 <- mcmc[,1:CORES,"TreatmentGH:TimeT2"]
  TreatmentGHLD2 <- mcmc[,1:CORES,"TreatmentGH_LD:TimeT2"]
  TreatmentNF2 <- (mcmc[,1:CORES,"TreatmentNF:TimeT2"])
  # Difference between treatments
  GHvsGHLD2 <- TreatmentGH2 - TreatmentGHLD2;
  GHvsNF2 <- TreatmentGH2 - TreatmentNF2; 
  GHLDvsNF2 <- TreatmentGHLD2 - TreatmentNF2;
  # posterior distributions
  deltas <- cbind(
    GH = as.vector(TreatmentGH)
    ,GHLD = as.vector(TreatmentGHLD)
    ,NF = as.vector(TreatmentNF)
    ,GHvsGHLD = as.vector(GHvsGHLD)
    ,GHvsNF = as.vector(GHvsNF)
    ,GHLDvsNF = as.vector(GHLDvsNF)
    ,GH2 = as.vector(TreatmentGH) + as.vector(TreatmentGH2)
    ,GHLD2 = as.vector(TreatmentGHLD) + as.vector(TreatmentGHLD2)
    ,NF2 = as.vector(TreatmentNF) + as.vector(TreatmentNF2)
    ,GHvsGHLD2 = as.vector(GHvsGHLD2)
    ,GHvsNF2 = as.vector(GHvsNF2)
    ,GHLDvsNF2 = as.vector(GHLDvsNF2)
  )
  } else {
  deltas <- cbind(
                    GH = as.vector(TreatmentGH)
                  ,GHLD = as.vector(TreatmentGHLD)
                  ,NF = as.vector(TreatmentNF)
                  ,GHvsGHLD = as.vector(GHvsGHLD)
                  ,GHvsNF = as.vector(GHvsNF)
                  ,GHLDvsNF = as.vector(GHLDvsNF)
  )
  }
}

res <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),3) # quantiles
p <- round(apply(deltas, 2, function(x) gf(x)),3) # P: confidence that the parameter is positive (P>0)
table<- rbind(res,"P"=p)

print(table)
cat("P>0: Proportion of positive posterior values, i.e. confidence that the effect is positive \n")
# TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
# cat("Effect of GH: ", median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]\n")
# cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentGH)*100,1),"%\n")
# 
# TreatmentGH_LD <- as.vector(mcmc[,1:CORES,"TreatmentGH_LD"])
# cat("Effect of GH_LD: ", median(TreatmentGH_LD),"[",quantile(TreatmentGH_LD, probs=0.025),";",quantile(TreatmentGH_LD, probs=0.975),"]\n")
# cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentGH_LD)*100,1),"%\n")

 # if(myvar != "Sum_SRP" & myvar != "Sum_NH4") { 
 # TreatmentNF <- as.vector(mcmc[,1:CORES,"TreatmentNF"])
 # cat("Effect of NF: ", median(TreatmentNF),"[",quantile(TreatmentNF, probs=0.025),";",quantile(TreatmentNF, probs=0.975),"]\n")
 # cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentNF)*100,1),"%\n")
 # }

# cat( "\n Posterior check: comparison between observed values (y) and predicted values (yrep) \n")
# cat("Ratio (y/yrep) should include the value 1\n")
# # ratio y vs yrep
# r=NULL;for (i in 1:nrow(yrep)){r[i]<-mean(y/yrep[i,], na.rm=TRUE)}
# cat("Ratio (y/yrep): ", median(r),"[",quantile(r, probs=0.025),";",quantile(r, probs=0.975),"]\n")
} #End loop 
sink()


# ## FIGURES POSTERIOR CHECK
# for (myvar in model.name) {
#   #myvar= "Sum_SRP"
# pdf(paste0("results/IMS2016_",myvar,"_check.pdf"),width=8, height = 8)
#   
#   ## LOAD MODEL
#   best = loo = mcmc = NULL
#   load(paste0("results/IMS2016_",myvar,".Rdata"))
#   #best <- get(myvar)
#   #mcmc <- as.array(best$stanfit)
#   
#   ## POSTERIOR CHECK
#   y <- best$y
#   yrep <- posterior_predict(best)
#   group <- best$data$Treatment
#   
# 
# #par(mfrow=c(2,2))
# #par(mar=c(3.5, 4.5, 4.5, 1.5))
# # # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
# plot1 <- ppc_stat(y, yrep)
# plot2 <- ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
# #color_scheme_set("teal")
# #ppc_stat_grouped(y, yrep, group)
# plot3 <- pp_check(y, yrep, ppc_dens_overlay)
# # pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
# # # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
# plot4 <- ppc_scatter_avg(y, yrep)
# # ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)
# 
# grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
# dev.off()
# } 







## FIGURE RESULTS




# pdf(paste0("results/IMS2016_excretion.pdf"),width=8, height = 8)
# #layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
# par(mfrow=c(1,2))
# #par(mar=c(3.5, 4.5, 4.5, 1.5))
# model.name <- c("Sum_SRP", "Sum_NH4")
# for (myvar in model.name) {
#   
#   ## LOAD MODEL
#   best = loo = mcmc = NULL
#   load(paste0("results/IMS2016_",myvar,".Rdata"))
#   #best <- get(myvar)
#   mcmc <- as.array(best$stanfit)
#   
#   # ## POSTERIOR CHECK
#   # y <- best$y
#   # yrep <- posterior_predict(best)
#   # group <- best$data$Treatment
#   
# #if(myvar == "Sum_SRP" & myvar == "Sum_NH4") {
#   SHAM <- mcmc[,1:CORES,"(Intercept)"]
#   TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
#   TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
#   # Difference between treatments
#   GHvsGHLD <- TreatmentGHLD - TreatmentGH;
#   
#   # posterior distributions
#   deltas <- cbind(GH = as.vector(TreatmentGH)
#                   ,GHLD = as.vector(TreatmentGHLD)
#                   ,GHvsGHLD = as.vector(GHvsGHLD)
#   )
#   
#   q <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),3) # quantiles
#   p <- round(apply(deltas, 2, function(x) gf(x)),3)
#   res <- t(rbind(q,p))
#   
#   plot(NULL, xlim=c(0.5,nrow(res)+.5),ylim=c(min(res), max(res)+.2),xaxt="n",ylab=expression(Beta ~ ": treatment effects"), xlab='',main= myvar,bty="n", xaxt="n")
#   axis(1, labels=c("GH", "GHLD", "GH-GHLD"), at = 1:nrow(res),cex=.7)
#   abline(h=0,lty=2)
#   points(1:nrow(res),res[1:nrow(res),"50%"], pch=16,cex=1.5)
#   segments(1:nrow(res),res[1:nrow(res),"2.5%"],1:nrow(res),res[1:nrow(res),"97.5%"])
#   segments(1:nrow(res),res[1:nrow(res),"25%"],1:nrow(res),res[1:nrow(res),"75%"],lwd=3)
#   for (i in 1:nrow(res)){
#     y=max(res[,i])+.15
#     text(i,y-.05,  paste0(res[i,"50%"]," [",res[i,"2.5%"],";",res[i,"97.5%"],"]; P = ",res[i,"p"]),cex=0.5)
#     #text(i,y-.1,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
#     #text(i,y-.15,"(i.e., do not overlap with 0)",cex=.4)
#   } # end loop i
# #} # end if
# 
# } # end loop myvar
# dev.off()




#pdf(paste0("results/IMS2016_nutrient.pdf"),width=8, height = 8)
#layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
par(mfrow=c(1,2))
#par(mar=c(3.5, 4.5, 4.5, 1.5))
model.name <- c("NH4", "PO4")#, "FluxNH4", "FluxPO4")
for (myvar in model.name) {
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2016_",myvar,".Rdata"))
  #best <- get(myvar)
  mcmc <- as.array(best$stanfit)
  
  # ## POSTERIOR CHECK
  # y <- best$y
  # yrep <- posterior_predict(best)
  # group <- best$data$Treatment
  
  
  #if(myvar != "Sum_SRP" & myvar != "Sum_NH4") {  
  SHAM <- mcmc[,1:CORES,"(Intercept)"]
  TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
  TreatmentNF <- (mcmc[,1:CORES,"TreatmentNF"])
  # Difference between treatments
  GHvsGHLD <- TreatmentGHLD - TreatmentGH;
  #GHvsNF <- TreatmentGH - TreatmentNF; 
  #GHLDvsNF <- TreatmentGHLD - TreatmentNF;
  
  # posterior distributions
  deltas <- cbind(GH = as.vector(TreatmentGH)
                  ,GHLD = as.vector(TreatmentGHLD)
                  ,NF = as.vector(TreatmentNF)
                  ,GHvsGHLD = as.vector(GHvsGHLD)
                  #,GHvsNF = as.vector(GHvsNF)
                  #,GHLDvsNF = as.vector(GHLDvsNF)
  )
  
  q <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),3) # quantiles
  p <- round(apply(deltas, 2, function(x) gf(x)),3)
  res <- t(rbind(q,p))
  
  # PLOT
  plot(NULL, xlim=c(0.5,nrow(res)+.5),ylim=c(min(res[,1:5]), max(res[,1:5])+.2),xaxt="n",ylab=expression(Beta ~ ": treatment effects"), xlab='',main= myvar,bty="n", xaxt="n")
  axis(1, labels=c("GH", "GHLD", "NF", "GH-GHLD"), at = 1:nrow(res),cex=.7)
  abline(h=0,lty=2)
  points(1:nrow(res),res[1:nrow(res),"50%"], pch=16,cex=1.5)
  segments(1:nrow(res),res[1:nrow(res),"2.5%"],1:nrow(res),res[1:nrow(res),"97.5%"])
  segments(1:nrow(res),res[1:nrow(res),"25%"],1:nrow(res),res[1:nrow(res),"75%"],lwd=3)
  for (i in 1:nrow(res)){
    y=max(res[i,1:5])+.05
    text(i,y,  paste0(res[i,"50%"]," [",res[i,"2.5%"],";",res[i,"97.5%"],"]; P = ",res[i,"p"]),cex=0.75)
    #text(i,y-.1,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
    #text(i,y-.15,"(i.e., do not overlap with 0)",cex=.4)
  } # end loop i






# Precictive values
nd <- data.frame(Time = c(rep("T0",8), rep("T1",8), rep("T2", 8)), Treatment = c(rep("SHAM",2),rep("GH",2),rep("GH_LD",2),rep("NF",2)), Position = rep(rev(levels(best$data$Position)),12))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
#boxplot(y_rep,col=rep(1:2,5))
y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
y.pred.tmp <- list()
y.pred.tmp[[1]] <- y.pred[,1:8] # predicted values at T0
y.pred.tmp[[2]] <- y.pred[,9:16] # predicted values at T1
y.pred.tmp[[3]] <- y.pred[,17:24] # predicted values at T2

point= TRUE # plot the data?
for (t in 1:3){
  gap <- c(-.3, -.1, .1, .3)
  mycol <- c(1, "#CD2626", "#CD262650", "lightgrey")
  plot(NULL,xlim=c(.5,2.5),ylim=c(min(y.pred),max(y.pred)),ylab=myvar,xlab="",xaxt="n",main=paste0("T",t-1),bty="n")
  axis(1, labels=rev(levels(best$data$Position)), at = 1:2)
  
  
  # SHAM
  points((1:2)+gap[1],y.pred.tmp[[t]]["50%",1:2],pch=16,col=mycol[1], cex=1.5)
  segments((1:2)+gap[1],y.pred.tmp[[t]]["25%",1:2],(1:2)+gap[1],y.pred.tmp[[t]]["75%",1:2],lwd=3,col=mycol[1])
  segments((1:2)+gap[1],y.pred.tmp[[t]]["2.5%",1:2],(1:2)+gap[1],y.pred.tmp[[t]]["97.5%",1:2],lwd=1,col=mycol[1])
  # if (point == TRUE){
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="SHAM" & Position=="Up")
  #   points(rep(1,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="SHAM" & Position=="Down")
  #   points(rep(3,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
  # }
  
  # GH
  points((1:2)+gap[2],y.pred.tmp[[t]]["50%",3:4],pch=16,col=mycol[2], cex=1.5)
  segments((1:2)+gap[2],y.pred.tmp[[t]]["25%",3:4],(1:2)+gap[2],y.pred.tmp[[t]]["75%",3:4],lwd=3,col=mycol[2])
  segments((1:2)+gap[2],y.pred.tmp[[t]]["2.5%",3:4],(1:2)+gap[2],y.pred.tmp[[t]]["97.5%",3:4],lwd=1,col=mycol[2])
  # if (point == TRUE){
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="GH" & Position=="Up")
  #   points(rep(1,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="GH" & Position=="Down")
  #   points(rep(3,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
  # }
  
  # GH_LD
  points((1:2)+gap[3],y.pred.tmp[[t]]["50%",5:6],pch=16,col=mycol[3], cex=1.5)
  segments((1:2)+gap[3],y.pred.tmp[[t]]["25%",5:6],(1:2)+gap[3],y.pred.tmp[[t]]["75%",5:6],lwd=3,col=mycol[3])
  segments((1:2)+gap[3],y.pred.tmp[[t]]["2.5%",5:6],(1:2)+gap[3],y.pred.tmp[[t]]["97.5%",5:6],lwd=1,col=mycol[3])
  # if (point == TRUE){
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="GH_LD" & Position=="Up")
  #   points(rep(1,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="GH_LD" & Position=="Down")
  #   points(rep(3,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
  # }
  
  # NF
  points((1:2)+gap[4],y.pred.tmp[[t]]["50%",7:8],pch=16,col=mycol[4], cex=1.5)
  segments((1:2)+gap[4],y.pred.tmp[[t]]["25%",7:8],(1:2)+gap[4],y.pred.tmp[[t]]["75%",7:8],lwd=3,col=mycol[4])
  segments((1:2)+gap[4],y.pred.tmp[[t]]["2.5%",7:8],(1:2)+gap[4],y.pred.tmp[[t]]["97.5%",7:8],lwd=1,col=mycol[4])
  # if (point == TRUE){
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="NF" & Position=="Up")
  #   points(rep(1,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
  #   tmp <- subset(best$data, Time==paste0("T",t-1) & Treatment =="NF" & Position=="Down")
  #   points(rep(3,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
  # }
  
  legend("topright",legend=c("SHAM", "GH", "GH_LD", "NF"),col=mycol,pch=rep(16,4),bty="n")
} #end loop t 

} # end loop myvar
dev.off()









mycol <- c("#6B6B6B", "#009ACD", "#5CACEE", "#6B6B6B50")

pdf(paste0("results/IMS2016_Fluxnutrient.pdf"),width=8, height = 8)
#layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
par(mfrow=c(3,3))
par(mar=c(3.5, 4.5, 4.5, 1.5))
model.name <- c("FluxNH4", "FluxPO4", "Flux_NNO3", "Flux_N_total")
for (myvar in model.name) {
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2016_",myvar,".Rdata"))
  #best <- get(myvar)
  mcmc <- as.array(best$stanfit)
  
  # ## POSTERIOR CHECK
  # # y <- best$y
  # # yrep <- posterior_predict(best)
  # # group <- best$data$Treatment
  # 
  # 
  # #if(myvar != "Sum_SRP" & myvar != "Sum_NH4") {  
  # SHAM <- mcmc[,1:CORES,"(Intercept)"]
  # TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  # TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
  # TreatmentNF <- (mcmc[,1:CORES,"TreatmentNF"])
  # # Difference between treatments
  # GHvsGHLD <- TreatmentGH - TreatmentGHLD;
  # GHvsNF <- TreatmentGH - TreatmentNF; 
  # GHLDvsNF <- TreatmentGHLD - TreatmentNF;
  # 
  # # posterior distributions
  # deltas <- cbind(GH = as.vector(TreatmentGH)
  #                 ,GHLD = as.vector(TreatmentGHLD)
  #                 ,NF = as.vector(TreatmentNF)
  #                 ,GHvsGHLD = as.vector(GHvsGHLD)
  #                 ,GHvsNF = as.vector(GHvsNF)
  #                 ,GHLDvsNF = as.vector(GHLDvsNF)
  # )
  # 
  # res <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),3) # quantiles
  # p <- round(apply(deltas, 2, function(x) gf(x)),3) # P: confidence that the parameter is positive (P>0)
  # 
  # # PLOT
  # plot(NULL, xlim=c(0.5,ncol(deltas)+.5),ylim=c(min(res), max(res)+.2),xaxt="n",ylab=expression(Beta ~ ": treatment effects"), xlab='',main= myvar,bty="n", xaxt="n")
  # axis(1, labels=c("GH", "GHLD", "NF", "GH-GHLD", "GH-NF","GHLD-NF"), at = 1:ncol(deltas),cex=.7)
  # abline(h=0,lty=2)
  # points(1:ncol(deltas),res["50%",1:ncol(deltas)], pch=16,cex=1.5)
  # segments(1:ncol(deltas),res["2.5%",1:ncol(deltas)],1:ncol(deltas),res["97.5%",1:ncol(deltas)])
  # segments(1:ncol(deltas),res["25%",1:ncol(deltas)],1:ncol(deltas),res["75%",1:ncol(deltas)],lwd=3)
  # for (i in 1:ncol(deltas)){
  #   y=max(res[,i])+.15
  #   text(i,y-.05,  paste0(res["50%",i]," [",res["2.5%",i],";",res["97.5%",i],"]; P = ",p[i]),cex=0.5)
  #   #text(i,y-.1,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
  #   #text(i,y-.15,"(i.e., do not overlap with 0)",cex=.4)
  # } # end loop i
  # # } # end if  
  
  
  
  
  
  
  # Precictive values
  if (myvar == "FluxNH4" || myvar == "FluxPO4") {  
    time=1:3
    nd <- data.frame(Time = c(rep("T0",4), rep("T1",4), rep("T2", 4)), Treatment = c(rep("SHAM",1),rep("GH",1),rep("GH_LD",1),rep("NF",1)))
    y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
    #boxplot(y_rep,col=rep(1:2,5))
    y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
    res <- cbind(nd, t(y.pred))
    # y.pred.tmp <- list()
    # y.pred.tmp[[1]] <- y.pred[,1:4] # predicted values at T0
    # y.pred.tmp[[2]] <- y.pred[,5:8] # predicted values at T1
    # y.pred.tmp[[3]] <- y.pred[,9:12] # predicted values at T2
    }
  
  
  if (myvar == "Flux_NNO3" || myvar == "Flux_N_total") {
  time=3
  nd <- data.frame(Time = c(rep("T2",4)), Treatment = c(rep("SHAM",1),rep("GH",1),rep("GH_LD",1),rep("NF",1)))
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  res <- cbind(nd, t(y.pred))
  # y.pred.tmp <- list()
  # y.pred.tmp[[1]] <- y.pred[,1:4] # predicted values at T0
  # y.pred.tmp[[2]] <- y.pred[,1:4] # predicted values at T1
  # y.pred.tmp[[3]] <- y.pred[,1:4] # predicted values at T2
  }
 
# 
#   plot(NULL,xlim=c(.5,nrow(res)+.5),ylim=c(min(y.pred),max(y.pred)),ylab=myvar,xlab="Flux (Down/up)",xaxt="n", bty="n",main=paste0("T",t-1))
#   axis(1, labels=res$Time, at = 1:nrow(res))
#   points(1:nrow(res),res[,"50%"],pch=16,col=mycol, cex=1.5)
#   segments(1:nrow(res),res[,"25%"],1:nrow(res),res[,"75%"],lwd=3,col=mycol)
#   segments(1:nrow(res),res[,"2.5%"],1:nrow(res),res[,"97.5%"],lwd=1,col=mycol)

  for (t in time){
    plot(NULL,xlim=c(.5,4+.5),ylim=c(min(y.pred),max(y.pred)),ylab=myvar,xlab="Flux (Down/up)",xaxt="n", bty="n",main=paste0("T",t-1))
    #axis(1, labels=res$Time, at = 1:nrow(res))
    points(1:4,res[1:4,"50%"],pch=16,col=mycol, cex=1.5)
    segments(1:4,res[,"25%"],1:4,res[1:4,"75%"],lwd=3,col=mycol)
    segments(1:4,res[,"2.5%"],1:4,res[1:4,"97.5%"],lwd=1,col=mycol)
    
    ## add data points
    j=0
    for (treatment in c("SHAM" , "GH" ,   "GH_LD", "NF")){
      j=j+1
    tmp <- best$y[which(best$data$Time==paste0("T",t-1) & best$data$Treatment ==treatment)]
    points(rep(j,5), tmp, col=mycol[j], pch=3)
    }
  } # end loop t
    #legend("topright",legend=c("SHAM", "GH", "GH_LD", "NF"),col=mycol,pch=rep(16,4),bty="n")

} # end loop myvar
dev.off()








pdf(paste0("results/IMS2016_Filamentous.pdf"),width=8, height = 8)
#layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow = TRUE))
par(mfrow=c(2,1))
#par(mar=c(3.5, 4.5, 4.5, 1.5))
  
  myvar <- "Fila"
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2016_",myvar,".Rdata"))
  #best <- get(myvar)
  mcmc <- as.array(best$stanfit)
  
  ## POSTERIOR CHECK
  y <- best$y
  yrep <- posterior_predict(best)
  group <- best$data$Treatment
  
  
  #if(myvar != "Sum_SRP" & myvar != "Sum_NH4") {  
  SHAM <- mcmc[,1:CORES,"(Intercept)"]
  TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
  TreatmentNF <- (mcmc[,1:CORES,"TreatmentNF"])
  # Difference between treatments
  GHvsGHLD <- TreatmentGHLD - TreatmentGH;
  # GHvsNF <- TreatmentGH - TreatmentNF; 
  # GHLDvsNF <- TreatmentGHLD - TreatmentNF;
  
  # posterior distributions
  deltas <- cbind(GH = as.vector(TreatmentGH)
                  ,GHLD = as.vector(TreatmentGHLD)
                  ,NF = as.vector(TreatmentNF)
                  ,GHvsGHLD = as.vector(GHvsGHLD)
                  #,GHvsNF = as.vector(GHvsNF)
                  #,GHLDvsNF = as.vector(GHLDvsNF)
  )
  
  q <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),3) # quantiles
  p <- round(apply(deltas, 2, function(x) gf(x)),3)
  res <- t(rbind(q,p))
  
  # PLOT
  plot(NULL, xlim=c(0.5,nrow(res)+.5),ylim=c(min(res[,1:5]), max(res[,1:5])+.2),xaxt="n",ylab=expression(Beta ~ ": treatment effects"), xlab='',main= myvar,bty="n", xaxt="n")
  axis(1, labels=c("GH", "GHLD", "NF", "GH-GHLD"), at = 1:nrow(res),cex=.7)
  abline(h=0,lty=2)
  points(1:nrow(res),res[1:nrow(res),"50%"], pch=16,cex=1.5)
  segments(1:nrow(res),res[1:nrow(res),"2.5%"],1:nrow(res),res[1:nrow(res),"97.5%"])
  segments(1:nrow(res),res[1:nrow(res),"25%"],1:nrow(res),res[1:nrow(res),"75%"],lwd=3)
  for (i in 1:nrow(res)){
    y=max(res[i,1:5])+.05
    text(i,y,  paste0(res[i,"50%"]," [",res[i,"2.5%"],";",res[i,"97.5%"],"]; P = ",res[i,"p"]),cex=0.75)
    #text(i,y-.1,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
    #text(i,y-.15,"(i.e., do not overlap with 0)",cex=.4)
  } # end loop i
  
  
  
  
  
  
  # Precictive values
  nd <- data.frame(Treatment = c(rep("SHAM",6),rep("GH",6),rep("GH_LD",6),rep("NF",6)), Position = rep(levels(best$data$Position),1))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  #boxplot(y_rep,col=rep(1:2,5))
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  y.pred.tmp <- list()
  y.pred.tmp[[1]] <- y.pred[,1:6] # predicted values for SHAM
  y.pred.tmp[[2]] <- y.pred[,7:12] # predicted values at T1
  y.pred.tmp[[3]] <- y.pred[,13:18] # predicted values at T2
  y.pred.tmp[[4]] <- y.pred[,19:24] # predicted values at T2
  
  point= TRUE # plot the data?
  
    gap <- c(-.3, -.1, .1, .3)
    mycol <- c(1, "#CD2626", "#CD262650", "lightgrey")
    plot(NULL,xlim=c(.5,6.5),ylim=c(min(y.pred),max(y.pred)),ylab="Filamentous (logit scale)",xlab="",xaxt="n",bty="n")
    axis(1, labels=(levels(best$data$Position)), at = 1:6)
    
    for (t in 1:4){
    # SHAM
    points((1:6)+gap[t],y.pred.tmp[[t]]["50%",1:6],pch=16,col=mycol[t], cex=1.5)
    segments((1:6)+gap[t],y.pred.tmp[[t]]["25%",1:6],(1:6)+gap[t],y.pred.tmp[[t]]["75%",1:6],lwd=3,col=mycol[t])
    segments((1:6)+gap[t],y.pred.tmp[[t]]["2.5%",1:6],(1:6)+gap[t],y.pred.tmp[[t]]["97.5%",1:6],lwd=1,col=mycol[t])
  #  points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], best$y, col= mycol[as.numeric(as.factor(best$data$Treatment))], pch=3)
  } #end loop t 
    #legend("topright",legend=c("SHAM", "GH", "GH_LD", "NF"),col=mycol,pch=rep(16,4),bty="n")
dev.off()






# pdf(paste0("results/IMS2016_Others.pdf"),width=8, height = 8)
# #layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
# par(mfrow=c(1,1))
# #par(mar=c(3.5, 4.5, 4.5, 1.5))
# point= FALSE # plot the data?
# model.name <- c("K_CM", "K_FM","Cyanobacteria","Green_algae","Diatoms","All_Chloro", "Nb_Chironomidae")
# for (myvar in model.name) {
#   
#   ## LOAD MODEL
#   best = loo = mcmc = NULL
#   load(paste0("results/IMS2016_",myvar,".Rdata"))
#   #best <- get(myvar)
#   mcmc <- as.array(best$stanfit)
#   
#   # ## POSTERIOR CHECK
#   # y <- best$y
#   # yrep <- posterior_predict(best)
#   # group <- best$data$Treatment
#   
#   
#   SHAM <- mcmc[,1:CORES,"(Intercept)"]
#   TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
#   TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
#   TreatmentNF <- (mcmc[,1:CORES,"TreatmentNF"])
#   # Difference between treatments
#   #GHvsGHLD <- TreatmentGHLD - TreatmentGH;
#   #GHvsNF <- TreatmentGH - TreatmentNF; 
#   #GHLDvsNF <- TreatmentGHLD - TreatmentNF;
#   
#   deltas <- cbind(
#     GH = as.vector(TreatmentGH)
#     ,GHLD = as.vector(TreatmentGHLD)
#     ,NF = as.vector(TreatmentNF)
#     #,GHvsGHLD = as.vector(GHvsGHLD)
#     #,GHvsNF = as.vector(GHvsNF)
#     #,GHLDvsNF = as.vector(GHLDvsNF)
#   )
#   
#   
#   q <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),2) # quantiles
#   p <- round(apply(deltas, 2, function(x) gf(x)),2)
#   res <- t(rbind(q,p))
#   
#   
#   if( paste0(best$formula)[3] == "Treatment * Time + Position_Bis + (1 | Block/Section)"){
#     TreatmentGH2 <- TreatmentGH + mcmc[,1:CORES,"TreatmentGH:TimeT2"]
#     TreatmentGHLD2 <- TreatmentGHLD + mcmc[,1:CORES,"TreatmentGH_LD:TimeT2"]
#     TreatmentNF2 <- TreatmentNF + (mcmc[,1:CORES,"TreatmentNF:TimeT2"])
#     # Difference between treatments
#     #GHvsGHLD2 <- TreatmentGHLD2 - TreatmentGH2;
#     #GHvsNF2 <- TreatmentGH2 - TreatmentNF2; 
#     #GHLDvsNF2 <- TreatmentGHLD2 - TreatmentNF2;
#     # posterior distributions
#     # deltas <- cbind(
#     #   GH = as.vector(TreatmentGH)
#     #   ,GHLD = as.vector(TreatmentGHLD)
#     #   ,NF = as.vector(TreatmentNF)
#       #,GHvsGHLD = as.vector(GHvsGHLD)
#       #,GHvsNF = as.vector(GHvsNF)
#       #,GHLDvsNF = as.vector(GHLDvsNF)
#       # ,GH2 = as.vector(TreatmentGH2)
#       # ,GHLD2 = as.vector(TreatmentGHLD2)
#       # ,NF2 = as.vector(TreatmentNF2)
#       # ,GHvsGHLD2 = as.vector(GHvsGHLD2)
#       # ,GHvsNF2 = as.vector(GHvsNF2)
#       # ,GHLDvsNF2 = as.vector(GHLDvsNF2)
#     #)
#     # q <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),2) # quantiles
#     # p <- round(apply(deltas, 2, function(x) gf(x)),2)
#     # res <- t(rbind(q,p))
#     
#     
#     deltas2 <- cbind(
#       GH2 = as.vector(TreatmentGH2)
#       ,GHLD2 = as.vector(TreatmentGHLD2)
#       ,NF2 = as.vector(TreatmentNF2)
#       #,GHvsGHLD2 = as.vector(GHvsGHLD2)
#       #,GHvsNF2 = as.vector(GHvsNF2)
#       #,GHLDvsNF2 = as.vector(GHLDvsNF2)
#     )
#     q2 <- round(apply(deltas2, 2, quantile , probs=c(.025, .25, .5, .75, .975)),2) # quantiles
#     p2 <- round(apply(deltas2, 2, function(x) gf(x)),2)
#     res2 <- t(rbind(q2,p2))
#     
#     
# 
#     # PLOT
#     gap <- c(-.05, .05)
#     ylim = c(min(res[,1:5], res2[,1:5]), max(res[,1:5], res2[,1:5])+.05)
#     plot(NULL, xlim=c(0.5,3.5),ylim=ylim,xaxt="n",ylab=expression(Beta ~ ": treatment effects"), xlab='',main= myvar,bty="n", xaxt="n")
#     axis(1, labels=c("GH", "GHLD", "NF"), at = 1:3,cex=.7)
#     abline(h=0,lty=2)
#     points((1:nrow(res))+gap[1],res[1:nrow(res),"50%"], pch=16,cex=1.5)
#     segments((1:nrow(res))+gap[1],res[1:nrow(res),"2.5%"],(1:nrow(res))+gap[1],res[1:nrow(res),"97.5%"])
#     segments((1:nrow(res))+gap[1],res[1:nrow(res),"25%"],(1:nrow(res))+gap[1],res[1:nrow(res),"75%"],lwd=3)
#     for (i in 1:nrow(res)){
#       y=max(res[i,1:5])+.05
#       text(i,y-.05,  paste0(res[i,"50%"]," [",res[i,"2.5%"],";",res[i,"97.5%"],"]; P = ",p[i]),cex=0.75)
#       #text(i,y-.1,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
#       #text(i,y-.15,"(i.e., do not overlap with 0)",cex=.4)
#     } # end loop i
# 
#     points((1:nrow(res2))+gap[2],res2[1:nrow(res2),"50%"], pch=16,cex=1.5)
#     segments((1:nrow(res2))+gap[2],res2[1:nrow(res2),"2.5%"],(1:nrow(res2))+gap[2],res2[1:nrow(res2),"97.5%"])
#     segments((1:nrow(res2))+gap[2],res2[1:nrow(res2),"25%"],(1:nrow(res2))+gap[2],res2[1:nrow(res2),"75%"],lwd=3)
#     for (i in 1:nrow(res)){
#       y=max(res2[i,1:5])+.05
#       text(i,y-.05,  paste0(res2[i,"50%"]," [",res2[i,"2.5%"],";",res2[i,"97.5%"],"]; P = ",p2[i]),cex=0.75)
#       #text(i,y-.1,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
#       #text(i,y-.15,"(i.e., do not overlap with 0)",cex=.4)
#     } # end loop i
#     
#     
#     
#     #legend("topright",legend=c("T1", "T2"),col=c(1,"#CD2626"),pch=rep(16,2),bty="n",cex=.75)
#     
#   } else {
#     
#     
#     
#     # # PLOT
#     # plot(NULL, xlim=c(0.5,ncol(deltas)+.5),ylim=c(min(res), max(res)+.2),xaxt="n",ylab=expression(Beta ~ ": treatment effects"), xlab='',main= myvar,bty="n", xaxt="n")
#     plot(NULL, xlim=c(0.5,3.5),ylim=c(min(res[,1:5]), max(res[,1:5])+0.05),xaxt="n",ylab=expression(Beta ~ ": treatment effects"), xlab='',main= myvar,bty="n", xaxt="n")
#     axis(1, labels=c("GH", "GHLD", "NF"), at = 1:3,cex=.7)
#     abline(h=0,lty=2)
#     points((1:nrow(res))+gap[1],res[1:nrow(res),"50%"], pch=16,cex=1.5)
#     segments((1:nrow(res))+gap[1],res[1:nrow(res),"2.5%"],(1:nrow(res))+gap[1],res[1:nrow(res),"97.5%"])
#     segments((1:nrow(res))+gap[1],res[1:nrow(res),"25%"],(1:nrow(res))+gap[1],res[1:nrow(res),"75%"],lwd=3)
#     for (i in 1:nrow(res)){
#       y=max(res[i,1:5])+.05
#       text(i,y-.05,  paste0(res[i,"50%"]," [",res[i,"2.5%"],";",res[i,"97.5%"],"]; P = ",p[i]),cex=0.75)
#       #text(i,y-.1,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
#       #text(i,y-.15,"(i.e., do not overlap with 0)",cex=.4)
#     } # end loop i
#   
#   } # end else
#   
# 
# 
# 
#   # # Precictive values
#   # nd <- data.frame(Time = c(rep("T1",12), rep("T2", 12)), Treatment = c(rep("SHAM",3),rep("GH",3),rep("GH_LD",3),rep("NF",3)), Position_Bis = rep(rev(levels(best$data$Position_Bis)),8))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
#   # y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
#   # #boxplot(y_rep,col=rep(1:2,5))
#   # y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
#   # y.pred.tmp <- list()
#   # y.pred.tmp[[1]] <- y.pred[,1:12] # predicted values at T1
#   # y.pred.tmp[[2]] <- y.pred[,13:24] # predicted values at T2
#   # 
#   # 
#   # for (t in 1:2){
#   # gap <- c(-.3, -.1, .1, .3)
#   # mycol <- c("#6B6B6B", "#009ACD", "#5CACEE", "lightgrey")
#   # plot(NULL,xlim=c(.5,3.5),ylim=c(min(y.pred),max(y.pred)),ylab=myvar,xlab="",xaxt="n",main=paste0("T",t),bty="n")
#   # axis(1, labels=rev(levels(best$data$Position_Bis)), at = 1:3)
#   # 
#   # 
#   # # SHAM
#   # points((1:3)+gap[1],y.pred.tmp[[t]]["50%",1:3],pch=16,col=mycol[1], cex=1.5)
#   # segments((1:3)+gap[1],y.pred.tmp[[t]]["25%",1:3],(1:3)+gap[1],y.pred.tmp[[t]]["75%",1:3],lwd=3,col=mycol[1])
#   # segments((1:3)+gap[1],y.pred.tmp[[t]]["2.5%",1:3],(1:3)+gap[1],y.pred.tmp[[t]]["97.5%",1:3],lwd=1,col=mycol[1])
#   # if (point == TRUE){
#   # tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="SHAM" & Position_Bis=="Upstream")
#   # points(rep(1,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
#   # tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="SHAM" & Position_Bis=="Middle")
#   # points(rep(2,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
#   # tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="SHAM" & Position_Bis=="Downstream")
#   # points(rep(3,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
#   # }
#   # 
#   # # GH
#   # points((1:3)+gap[2],y.pred.tmp[[t]]["50%",4:6],pch=16,col=mycol[2], cex=1.5)
#   # segments((1:3)+gap[2],y.pred.tmp[[t]]["25%",4:6],(1:3)+gap[2],y.pred.tmp[[t]]["75%",4:6],lwd=3,col=mycol[2])
#   # segments((1:3)+gap[2],y.pred.tmp[[t]]["2.5%",4:6],(1:3)+gap[2],y.pred.tmp[[t]]["97.5%",4:6],lwd=1,col=mycol[2])
#   # if (point == TRUE){
#   # tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH" & Position_Bis=="Upstream")
#   # points(rep(1,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
#   # tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH" & Position_Bis=="Middle")
#   # points(rep(2,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
#   # tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH" & Position_Bis=="Downstream")
#   # points(rep(3,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
#   # }
#   # 
#   # # GH_LD
#   # points((1:3)+gap[3],y.pred.tmp[[t]]["50%",7:9],pch=16,col=mycol[3], cex=1.5)
#   # segments((1:3)+gap[3],y.pred.tmp[[t]]["25%",7:9],(1:3)+gap[3],y.pred.tmp[[t]]["75%",7:9],lwd=3,col=mycol[3])
#   # segments((1:3)+gap[3],y.pred.tmp[[t]]["2.5%",7:9],(1:3)+gap[3],y.pred.tmp[[t]]["97.5%",7:9],lwd=1,col=mycol[3])
#   # if (point == TRUE){
#   #   tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH_LD" & Position_Bis=="Upstream")
#   #   points(rep(1,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
#   #   tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH_LD" & Position_Bis=="Middle")
#   #   points(rep(2,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
#   #   tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH_LD" & Position_Bis=="Downstream")
#   #   points(rep(3,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
#   # }
#   # 
#   # # NF
#   # points((1:3)+gap[4],y.pred.tmp[[t]]["50%",10:12],pch=16,col=mycol[4], cex=1.5)
#   # segments((1:3)+gap[4],y.pred.tmp[[t]]["25%",10:12],(1:3)+gap[4],y.pred.tmp[[t]]["75%",10:12],lwd=3,col=mycol[4])
#   # segments((1:3)+gap[4],y.pred.tmp[[t]]["2.5%",10:12],(1:3)+gap[4],y.pred.tmp[[t]]["97.5%",10:12],lwd=1,col=mycol[4])
#   # if (point == TRUE){
#   #   tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="NF" & Position_Bis=="Upstream")
#   #   points(rep(1,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
#   #   tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="NF" & Position_Bis=="Middle")
#   #   points(rep(2,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
#   #   tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="NF" & Position_Bis=="Downstream")
#   #   points(rep(3,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
#   # }
#   # 
#   # legend("topright",legend=c("SHAM", "GH", "GH_LD", "NF"),col=mycol,pch=rep(16,4),bty="n", cex=.75)
#   # } #end loop t 
#   
# 
# 
# } # end loop
# dev.off()



