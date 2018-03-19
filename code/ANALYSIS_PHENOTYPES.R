############ Statistics of the phenotypic effects of GH-treatment #############

rm(list=ls())   # Clear memory

SEED <- 1234
set.seed(SEED)

### PACKAGES ###
library(readxl)
library(readr)
# see https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library("rstanarm")
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
#library(bayesplot)

### FUNCTIONS ###
invlogit<-function(x) {1/(1+exp(-(x)))} # inverse logit
logit<-function(x) {log(x/(1-x))} # logisitic transformation
gf <- function(x){if(median(x)>=0){mean(x<=0)}else{mean(x>0)}} # function to calculate the proportion of posterior with different sign as mean

### DIRECTORY ###
setwd("~/Documents/RESEARCH/PROJECTS/SALMOINVADE/GrowthEnhancement")
#setwd("~/Documents/mbuoro/Ims/")

### ANALYSIS ###
CHAINS  = 3 # number of chains
CORES = CHAINS # number of core
ITER = 10000 # number of iterations
WARM = 2000 # warmup
ndelta = 0.99




#### Exp_1_Mass_Growth_Hatchery_Conditions ####
# Exp_1_Mass_Growth_Hatchery_Conditions => Fig 2a and 2b + Extended data Table 1
# Response variables: Body_Mass and Growth_Rate_SGRW

## Load dataset
df <- read_csv("data/Exp_1_Mass_Growth_Hatchery_Conditions.csv")
df <- as.data.frame(df)
df$Treatment <- as.factor(df$Treatment)
#levels(df$Treatment)[levels(df$Treatment)=="sham"] <- "SHAM"
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference



# 1. Body Mass

## Analysis
fit_mass_indoor <- stan_glmer(
  Final_Mass ~ Treatment  + (1| Tank_ID),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

## Results
best <- fit_mass_indoor
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE; TRUE if convergence failure
  model.name = "Exp_1_Mass_Hatchery_Conditions"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
}

## Posterior for GH treatment
mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)


# 2.  Growth

## Analysis
fit_growth_indoor <- stan_glmer(
  Growth_Rate_SGRW ~ Treatment  + (1| Tank_ID  ),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

## RESULTS
best <- fit_growth_indoor
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE; TRUE if convergence failure
  model.name = "Exp_1_Growth_Hatchery_Conditions"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
}


## Posterior for GH treatment
mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)





#### Exp_1_Mass_Growth_Experimental_Streams ####
# Exp_1_Mass_Growth_Experimental_Streams => Fig 2a and 2b + Extended data Table 1
# Response variables: Body_Mass, Growth_Rate_SGRW_Before_Release and Growth_Rate_SGRW_After_Release

# Load dataset
df <- read_csv("data/Exp_1_Mass_Growth_Experimental_Streams.csv")
df <- as.data.frame(df)
df$Treatment <- as.factor(df$Treatment)
df <- subset(df, Treatment != "no tag") # Remove unidentified fish
#levels(df$Treatment)[levels(df$Treatment)=="Sham"] <- "SHAM"
df <- droplevels(df)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference
df <- subset(df, Status == "Recaptured") # keep recaptured fish only



# 1. Body Mass

## Analysis
fit_mass_outdoor <- stan_glmer(
  Body_Mass ~ Treatment  + (1| Channel/Section),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## Results
best <- fit_mass_outdoor
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE; TRUE if convergence failure
  model.name = "Exp_1_Mass_Experimental_Streams"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
}

## Posterior for GH treatment
mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)



# 2. Growth Rate BEFORE Release

## Analysis
Growth_Rate_SGRW_Before_Release <- stan_glmer(
  Growth_Rate_SGRW_Before_Release ~ Treatment  + (1| Channel/Section),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## Results
best <- Growth_Rate_SGRW_Before_Release
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE; TRUE if convergence failure
  model.name = "Exp_1_Growth_Experimental_Streams_BEFORE_Release"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
}

## Posterior for GH treatment
mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)



# 3. Growth Rate AFTER Release

## Analysis
Growth_Rate_SGRW_After_Release <- stan_glmer(
  Growth_Rate_SGRW_After_Release ~ Treatment  + (1| Channel/Section),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## Results
best <- Growth_Rate_SGRW_After_Release
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE; TRUE if convergence failure
  model.name = "Exp_1_Growth_Experimental_Streams_AFTER_Release"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
}

## Posterior for GH treatment
mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)







#### MORPHOLOGY ####
# Exp_1_Morphology_Experimental_Streams => Fig 2c + Extended data Table 1 + Extended data Figure 2
# Response variables: WARP1 and WARP2

# Load dataset
df <- read_csv("data/Exp_1_Morphology_Experimental_Streams.csv")
df <- as.data.frame(df)
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference


# 1. WARP1

## Analysis
# WITH Interaction :
fit1 <- stan_glmer(
  WARP1 ~ Treatment * log(Body_Mass) + (1| Channel/Section), 
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

# WITHOUT Interaction :
fit2 <- stan_glmer(
  WARP1 ~ Treatment + log(Body_Mass) + (1| Channel/Section), 
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## Compare models using Leave-one-out cross-validation (LOO)
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment * log(Body_Mass)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment + log(Body_Mass)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## Results
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_1_Morphology_Experimental_Streams_WARP1"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}

## Posterior for GH treatment
mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)


# 2. WARP2

## Analysis
# WITH Interaction :
fit1 <- stan_glmer(
  WARP2 ~ Treatment * log(Body_Mass) + (1| Channel/Section),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

# WITHOUT Interaction :
fit2 <- stan_glmer(
  WARP2 ~ Treatment + log(Body_Mass) + (1| Channel/Section),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)


## Compare models using Leave-one-out cross-validation (LOO)
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment * log(Body_Mass)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment + log(Body_Mass)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## Results
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_1_Morphology_Experimental_Streams_WARP2"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}

## Posterior for GH treatment
mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)










######## ACTIVITY ######
# Exp_2_Activity_Open_Field_Tests => Fig 2d + Extended data Table 1 
# Response variables: Activity

# Load dataset
df <- read_csv("data/Exp_2_Activity_Open_Field_Tests.csv")
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference



## Analysis
fit1 <- stan_glmer(
  Activity ~ log(Body_Mass)+Treatment+Scoring_Session + (1|Fish_ID),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Activity ~ log(Body_Mass)*Treatment+Scoring_Session + (1|Fish_ID),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit3 <- stan_glmer(
  Activity ~ log(Body_Mass)+Treatment*Scoring_Session + (1|Fish_ID),
  data = df,
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  Activity ~ log(Body_Mass)*Treatment*Scoring_Session + (1|Fish_ID),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models using Leave-one-out cross-validation (LOO)
loo1 <- loo(fit1,k_threshold = 0.7) # log(Body_Mass)+Treatment+scoring
loo2 <- loo(fit2,k_threshold = 0.7) # log(Body_Mass)*Treatment+scoring
loo3 <- loo(fit3,k_threshold = 0.7) # log(Body_Mass)+Treatment*scoring
loo4 <- loo(fit4,k_threshold = 0.7) # log(Body_Mass)*Treatment*scoring
loo <- compare(loo1, loo2, loo3, loo4)

# Select the best model
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
## here, model with (fit2) and without (fit1) interaction provide same values of loo. We chose the model without interaction (fit1) for simplicity.
best <- fit1
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_2_Activity_Open_Field_Tests"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}


## Posterior for GH treatment
cat( "\n Posterior distribution: \n")
mcmc <- as.array(best$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)










#### Exp_2_Movement_Habitat_Use_Stream_Mesocosms ####
# Exp_2_Movement_Habitat_Use_Stream_Mesocosms => Fig 2e + Extended data Table 1
# Response variables: Movement and Habitat_Use

# Load dataset
df <- read_csv("data/Exp_2_Movement_Habitat_Use_Stream_Mesocosms.csv")
df <- as.data.frame(df)
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference

# 1. Movement

## Analysis
fit1 <- stan_glmer(
  Movement ~ Treatment+Time_of_the_day+log(Body_mass) + (1|Tracking_session:ID) + (1|Stream_mesocosm),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Movement ~ Treatment*log(Body_mass) + Time_of_the_day + (1|Tracking_session:ID) + (1|Stream_mesocosm),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit3 <- stan_glmer(
  Movement ~ Treatment* Time_of_the_day + log(Body_mass) + (1|Tracking_session:ID) + (1|Stream_mesocosm),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  Movement ~ Treatment * Time_of_the_day * log(Body_mass) + (1|Tracking_session:ID) + (1|Stream_mesocosm),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # log(Body_mass)+treatment+time
loo2 <- loo(fit2,k_threshold = 0.7) # log(Body_mass)*treatment+time
loo3 <- loo(fit3,k_threshold = 0.7) # log(Body_mass)+treatment*time
loo4 <- loo(fit4,k_threshold = 0.7) # log(Body_mass)*treatment*time
loo <- compare(loo1, loo2, loo3, loo4) 

# Select the best model
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_2_Movement_Stream_Mesocosms"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}


cat( "\n Posterior distribution: \n")
mcmc <- as.array(best$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat(median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]","\n")

cat( "\n The proportion of the posterior with the same sign as the mean: ")
f <- gf(TreatmentGH)
cat(round(f*100,1),"% \n")





# 2. Habitat Use

fit1 <- stan_glmer(
  Habitat_Use ~ Treatment+Time_of_the_day+log(Body_mass) + (1|Tracking_session:ID) + (1|Stream_mesocosm),
  family = binomial(),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Habitat_Use ~ Treatment*Time_of_the_day+log(Body_mass) + (1|Tracking_session:ID) + (1|Stream_mesocosm),
  family = binomial(),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # log(Body_mass)+Treatment+time
loo2 <- loo(fit2,k_threshold = 0.7) # log(Body_mass)*Treatment+time
loo <- compare(loo1, loo2)

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_2_HabitatUse_Stream_Mesocosms"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
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












#### Exp_3_Excretion_Hatchery_Conditions #####
# Exp_3_Excretion_Hatchery_Conditions => Extended data Table 1
# Response variables: N_Excretion and P_Excretion 

# Load dataset
df <- read_csv("data/Exp_3_Excretion_Hatchery_Conditions.csv")
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference



# 1. P_Excretion
P_Excretion <- stan_glmer(
  log(P_Excretion) ~ Treatment + log(Body_Mass) + (1|Tank), # with mass as covariate
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## RESULTS
best <- P_Excretion
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Excretion_Hatchery_Conditions_P_Excretion"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
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



# 2. N_Excretion
N_Excretion <- stan_glmer(
  log(N_Excretion) ~ Treatment + log(Body_Mass) + (1|Tank), # with mass as covariate
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## RESULTS
best <- N_Excretion
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Excretion_Hatchery_Conditions_N_Excretion"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
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







#### Exp_3_Excretion_Stream_Mesocosms ####
# Exp_3_Excretion_Stream_Mesocosms => Fig 2f + Extended data Table 1
# Response variables: N_Excretion and P_Excretion 

# Load dataset
df <- read_csv("data/Exp_3_Excretion_Stream_Mesocosms.csv")
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference

# excretion_outdoor <- within(excretion_outdoor, Treatment <- relevel(Treatment, ref = 3)) # SHAM as reference
# excretion_outdoor <- within(excretion_outdoor, Treatment_bis <- relevel(Treatment_bis, ref = 2)) # SHAM as reference
# excretion_outdoor$Treatment <- excretion_outdoor$Treatment_bis


# 1. P_Excretion

P_Excretion <- stan_glmer(
  log(P_Excretion) ~ Treatment + log(Body_Mass) + (1|Block/Stream_Mesocosm), # with mass as covariate
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## RESULTS
best <- P_Excretion
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Excretion_Stream_Mesocosms_P_excretion"
  ### SAVE
  save(best, file=paste0("results/",model.name,".Rdata"))
}

cat( "\n Posterior distribution: \n")
mcmc <- as.array(best$stanfit)
TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
cat(
  round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
  ," [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),3),"]"
  , ", P="
  , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
)




# 2. N_Excretion
N_Excretion <- stan_glmer(
  log(N_Excretion) ~ Treatment + log(Body_Mass) + (1|Block/Stream_Mesocosm), # with mass as covariate
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## RESULTS
best <- N_Excretion
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Excretion_Stream_Mesocosms_N_excretion"
  ### SAVE
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






