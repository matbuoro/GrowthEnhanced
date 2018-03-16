rm(list=ls())   # Clear memory

SEED <- 1234
set.seed(SEED)

#### PACKAGES ####
library(readxl)
library(readr)
# see https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library("rstanarm")
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(bayesplot)

#### FUNCTIONS ####
invlogit<-function(x) {1/(1+exp(-(x)))} # inverse logit
logit<-function(x) {log(x/(1-x))} # logisitic transformation
gf <- function(x){if(median(x)>=0){mean(x<=0)}else{mean(x>0)}} # function to calculate the proportion of posterior with different sign as mean

#### DIRECTORY ####
setwd("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/")
#setwd("~/Documents/mbuoro/Ims/")

#### ANALYSIS ####
CHAINS  = 3 # number of chains
CORES = CHAINS # number of core
ITER = 10000 # number of iterations
WARM = 2000 # warmup
ndelta = 0.99








############ Statistics of the phenotypic effects of GH-treatment #############

#### Body Mass ####

#### Indoor
data_pheno_indoor <- read_excel("data/Growth_Mass_Length_Ims_2015_all.xlsx", 
                                sheet = "Indoor")
data_pheno_indoor <- as.data.frame(data_pheno_indoor)
data_pheno_indoor$Treatment <- as.factor(data_pheno_indoor$Treatment)
levels(data_pheno_indoor$Treatment)[levels(data_pheno_indoor$Treatment)=="sham"] <- "SHAM"
data_pheno_indoor <- within(data_pheno_indoor, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference


## ANALYSIS
fit_mass_indoor <- stan_glmer(
  Final_Mass ~ Treatment  + (1| Tank_ID),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = data_pheno_indoor,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

## RESULTS
if (any(summary(fit_mass_indoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_Mass_Indoor"
  ### SAVE
  best <- fit_mass_indoor
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}

mcmc <- as.array(fit_mass_indoor$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)


#### OUTDOOR
data_pheno_outdoor <- read_excel("data/Growth_Mass_Length_Ims_2015_all.xlsx", 
                                 sheet = "River park")
data_pheno_outdoor <- as.data.frame(data_pheno_outdoor)
data_pheno_outdoor$Treatment <- as.factor(data_pheno_outdoor$Treatment)
data_pheno_outdoor$Fil_Mass <- as.numeric(data_pheno_outdoor$Fil_Mass)
data_pheno_outdoor$Growth_Rate_SGRW_Tagging_Release <- as.numeric(data_pheno_outdoor$Growth_Rate_SGRW_Tagging_Release)
data_pheno_outdoor$Growth_Rate_SGRW_Release_Recapture <- as.numeric(data_pheno_outdoor$Growth_Rate_SGRW_Release_Recapture)
levels(data_pheno_outdoor$Treatment)[levels(data_pheno_outdoor$Treatment)=="sham"] <- "SHAM"
data_pheno_outdoor <- na.omit(data_pheno_outdoor)
data_pheno_outdoor <- subset(data_pheno_outdoor, data_pheno_outdoor$Treatment !="no tag")
data_pheno_outdoor <- droplevels(data_pheno_outdoor)
data_pheno_outdoor <- within(data_pheno_outdoor, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference

## Remove all movers
data_pheno_outdoor <- subset(data_pheno_outdoor, Status == "Stay")

## ANALYSIS
fit_mass_outdoor <- stan_glmer(
  Fil_Mass ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = data_pheno_outdoor,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## RESULTS
if (any(summary(fit_mass_outdoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_Mass_Outdoor"
  ### SAVE
  best <- fit_mass_outdoor
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}

mcmc <- as.array(fit_mass_outdoor$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)


#### GROWTH ####

## INDOOR
fit_growth_indoor <- stan_glmer(
  Growth_Rate_SGRW ~ Treatment  + (1| Tank_ID  ),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = data_pheno_indoor,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

## RESULTS
if (any(summary(fit_growth_indoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Growth_Indoor"
  ### SAVE
  best <- fit_growth_indoor
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}

mcmc <- as.array(fit_growth_indoor$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)


## OUTDOOR

## 1. Growth measured between tagging and release in river park
Growth_Outdoor_TaggingRelease <- stan_glmer(
  Growth_Rate_SGRW_Tagging_Release ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = data_pheno_outdoor,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## RESULTS
if (any(summary(Growth_Outdoor_TaggingRelease)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Growth_Outdoor_TaggingRelease"
  ### SAVE
  best <- Growth_Outdoor_TaggingRelease
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}

mcmc <- as.array(Growth_Outdoor_TaggingRelease$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)

# 2. Growth measured between release in river park and recapture
Growth_Outdoor_ReleaseRecapture <- stan_glmer(
  Growth_Rate_SGRW_Release_Recapture ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = data_pheno_outdoor,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## RESULTS
if (any(summary(Growth_Outdoor_ReleaseRecapture)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Growth_Outdoor_ReleaseRecapture"
  ### SAVE
  best <- Growth_Outdoor_ReleaseRecapture
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}

mcmc <- as.array(Growth_Outdoor_ReleaseRecapture$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)


#### MORPHOLOGY ####
Data_Morpho_2015_Experimental_Streams <- read_csv("data/Data_Morpho_2015_Experimental_Streams.csv")
Data_Morpho_2015_Experimental_Streams <- as.data.frame(Data_Morpho_2015_Experimental_Streams)
Data_Morpho_2015_Experimental_Streams$Treatment <- as.factor(Data_Morpho_2015_Experimental_Streams$Treatment)
Data_Morpho_2015_Experimental_Streams <- within(Data_Morpho_2015_Experimental_Streams, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference


# 1. WARP1_outdoor
fit1 <- stan_glmer(
  WARP1 ~ Treatment * log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = Data_Morpho_2015_Experimental_Streams,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

fit2 <- stan_glmer(
  WARP1 ~ Treatment + log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = Data_Morpho_2015_Experimental_Streams,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
loo <- compare(loo1, loo2) 

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
WARP1_outdoor <- best
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "WARP1_outdoor"
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
}

mcmc <- as.array(WARP1_outdoor$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)


# 2. WARP2_outdoor
fit1 <- stan_glmer(
  WARP2 ~ Treatment * log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = Data_Morpho_2015_Experimental_Streams,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

fit2 <- stan_glmer(
  WARP2 ~ Treatment + log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = Data_Morpho_2015_Experimental_Streams,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
loo <- compare(loo1, loo2) 

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
WARP2_outdoor <- best
if (any(summary(WARP2_outdoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "WARP2_outdoor"
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
}

mcmc <- as.array(WARP2_outdoor$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)



# ####### ACTIVITY ######
# ####laboratory scoring of activity
# salmoGH <- read.table("data/lab_activity_2.csv", header=TRUE, sep=";", na.strings="NA", dec=",", strip.white=TRUE)
# salmoGH$group <- as.factor(salmoGH$group)
# salmoGH$ID <- as.factor(salmoGH$ID)
# salmoGH$activity_round <- as.factor(salmoGH$activity_round)
# salmoGH$activity_tank <- as.factor(salmoGH$activity_tank)
# salmoGH <- within(salmoGH, hormone <- relevel(hormone, ref = 2)) # SHAM as reference
# 
# #####11TH ROUND IS SIGNIFICANTLY DIFFERENT FROM THE REST -  THIS WAS ADDITIONAL ROUND AND CONTAINED INDIVIDUAS TREATED IN NON-STANDARD WAY -> REMOVED FROM THE DATASET
# ##AFTER REMOVING 11TH SCORING ROUND EFFECT OF THE ROUND BECOME NON-SIGNIFICANT
# ###effect of GH on laboratory activity in the open field test
# salmoGH <-subset(salmoGH, Distance_moved !="NA")
# salmoGH <-subset(salmoGH, activity_round !="11")
# 
# 
# 
# ## MCMC (modified by Mat)
# #hist(salmoGH$Distance_moved)
# fit1 <- stan_glmer(
#   Distance_moved ~ log(initial_weight)+hormone+scoring + (1|ID),
#   #family = poisson(link = log),
#   data = salmoGH,
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit2 <- stan_glmer(
#   Distance_moved ~ log(initial_weight)*hormone+scoring + (1|ID),
#   #family = poisson(link = log),
#   data = salmoGH,
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit3 <- stan_glmer(
#   Distance_moved ~ log(initial_weight)+hormone*scoring + (1|ID),
#   #family = poisson(link = log),
#   data = salmoGH,
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit4 <- stan_glmer(
#   Distance_moved ~ log(initial_weight)*hormone*scoring + (1|ID),
#   #family = poisson(link = log),
#   data = salmoGH,
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # log(initial_weight)+hormone+scoring
# loo2 <- loo(fit2,k_threshold = 0.7) # log(initial_weight)*hormone+scoring
# loo3 <- loo(fit3,k_threshold = 0.7) # log(initial_weight)+hormone*scoring
# loo4 <- loo(fit4,k_threshold = 0.7) # log(initial_weight)*hormone*scoring
# loo <- compare(loo1, loo2, loo3, loo4) 
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# ## here, model with (fit2) and without (fit1) interaction provide same values of loo. We chose the model without interaction (fit1) for simplicity.
# Activity <- fit1
# if (any(summary(Activity)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Activity"
#   ### SAVE
#   save(Activity,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
# }
# 
# 
# ## RESULTS
# cat( "\n Posterior distribution: \n")
# mcmc <- as.array(Activity$stanfit)
# TreatmentGH <- as.vector(mcmc[,1:CORES,"hormoneGHA"])
# cat(
#   round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
#   ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
#   , ", P="
#   , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
# )



#### EXCRETION #####

## INDOOR
excretion_indoor <- read.table("data/EXCRETION_Indoor_Mat.txt",h=TRUE) # INDOOR
excretion_indoor <- within(excretion_indoor, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference

## SRP
SRP_indoor <- stan_glmer(
  log(Per_Capita_SRP) ~ Treatment + log(Mass_t2) + (1|Tank), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = excretion_indoor,
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## RESULTS
if (any(summary(SRP_indoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "SRP_indoor"
  ### SAVE
  best <- SRP_indoor
  save(best, file=paste0("results/IMS2016_",model.name,".Rdata"))
}

cat( "\n Posterior distribution: \n")
mcmc <- as.array(SRP_indoor$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)



## NH4
NH4_indoor <- stan_glmer(
  #log(Per_Gram_NH4) ~ Treatment + (1|Tank),
  log(Per_Capita_NH4) ~ Treatment + log(Mass_t2) + (1|Tank),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = excretion_indoor,
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)
## RESULTS
if (any(summary(NH4_indoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "NH4_indoor"
  ### SAVE
  best <- NH4_indoor
  save(best, file=paste0("results/IMS2016_",model.name,".Rdata"))
}

cat( "\n Posterior distribution: \n")
mcmc <- as.array(best$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)




## OUTDOOR
excretion_outdoor <- read.table("data/EXCRETION_Mat.txt",h=TRUE) # OUTDOOR
excretion_outdoor <- as.data.frame(excretion_outdoor)

excretion_outdoor <- within(excretion_outdoor, Treatment <- relevel(Treatment, ref = 3)) # SHAM as reference
excretion_outdoor <- within(excretion_outdoor, Treatment_bis <- relevel(Treatment_bis, ref = 2)) # SHAM as reference
excretion_outdoor$Treatment <- excretion_outdoor$Treatment_bis


SRP_outdoor <- stan_glmer(
  log(Per_Capita_SRP) ~ Treatment + log(Mass_t2) + (1|Block/Section), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = excretion_outdoor,
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## RESULTS
if (any(summary(SRP_outdoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "SRP_outdoor"
  ### SAVE
  best <- SRP_outdoor
  save(best, file=paste0("results/IMS2016_",model.name,".Rdata"))
}

cat( "\n Posterior distribution: \n")
mcmc <- as.array(SRP_outdoor$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),3),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)




## NH4
NH4_outdoor <- stan_glmer(
  log(Per_Capita_NH4) ~ Treatment + log(Mass_t2) + (1|Block/Section), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = excretion_outdoor,
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)
## RESULTS
if (any(summary(NH4_outdoor)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "NH4_outdoor"
  ### SAVE
  best <- NH4_outdoor
  save(best, file=paste0("results/IMS2016_",model.name,".Rdata"))
}

cat( "\n Posterior distribution: \n")
mcmc <- as.array(NH4_outdoor$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)






