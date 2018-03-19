############ GH-TREATMENT EFFECTS ON COMMUNITY AND ECOSYSTEMS #############

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
library(bayesplot)

### FUNCTIONS ###
invlogit<-function(x) {1/(1+exp(-(x)))} # inverse logit
logit<-function(x) {log(x/(1-x))} # logisitic transformation
gf <- function(x){if(median(x)>=0){mean(x<=0)}else{mean(x>0)}} # function to calculate the proportion of posterior with different sign as mean

### DIRECTORY ###
#setwd("~/path/GrowthEnhancement")

### ANALYSIS ###
CHAINS  = 3 # number of chains
CORES = CHAINS # number of core
ITER = 10000 # number of iterations
WARM = 2000 # warmup
ndelta = 0.99

### DATA ###
## First, Choose experiment's dataset
data.name <- "With_Fish" # Exp_1_Community_Ecosystem_With_Fish_Experimental_Streams => Fig 3, Extended data Figure 4 + Extended data Table 2
#data.name <- "Surber_With_Fish" # Exp_1_Surber_With_Fish_Experimental_Streams => Extended data Figure 3+ Extended data Table 3
#data.name <- "Without_Fish" # Exp_1_Community_Ecosystem_Without_Fish_Experimental_Streams => Extended data Figure 5 + Extended data Table 2
#data.name <- "Surber_Without_Fish" #Exp_1_Surber_Without_Fish_Experimental_Streams => Extended data Table 3

## Then, load dataset
df <- read_csv(paste0("data/Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams.csv"))
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference


# 1. Rhyacophilidae
fit1<- stan_glmer(
  Rhyacophilidae ~ Treatment * Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  Rhyacophilidae ~ Treatment + Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Compare models using Leave-one-out cross-validation (LOO)
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE; TRUE if convergence failure
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Rhyacophilidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
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





# 2. Polycentropodidae
fit1<- stan_glmer(
  Polycentropodidae ~ Treatment * Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  Polycentropodidae ~ Treatment + Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Compare models using Leave-one-out cross-validation (LOO)
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE; TRUE if convergence failure
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Polycentropodidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
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



# 3. Total_Prim_Cons
fit1 <- stan_glmer(
  Total_Prim_Cons ~ Treatment*Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Total_Prim_Cons ~ Treatment+Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Total_Prim_Cons")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}



# 4. Chironomidae
fit1 <- stan_glmer(
  Chironomidae ~ Treatment*Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Chironomidae ~ Treatment + Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Chironomidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}


# 5. Hydropsychidae
fit1 <- stan_glmer(
  Hydropsychidae ~ Treatment*Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Hydropsychidae ~ Treatment + Position + (1 | Channel), 
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Hydropsychidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}


# 6. Baetidae
fit1 <- stan_glmer(
  Baetidae ~ Treatment*Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Baetidae ~ Treatment + Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Baetidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}


fit1 <- stan_glmer(
  Simuliidae ~ Treatment*Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Simuliidae ~ Treatment + Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2)#, loo3) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)#, loo3$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Baetidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}


# 7. Simuliidae
fit1 <- stan_glmer(
  Simuliidae ~ Treatment*Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Simuliidae ~ Treatment + Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2)

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Simuliidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}



# 8. Planorbidae
fit1 <- stan_glmer(
  Planorbidae ~ Treatment*Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Planorbidae ~ Treatment + Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2)

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Planorbidae")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}



# 9. Primary_production
fit1 <- stan_glmer(
  log(Primary_production) ~ Treatment*Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  log(Primary_production) ~ Treatment+Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2)

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Primary_production")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}



# 10. Total_Decomposition
fit1 <- stan_glmer(
  Total_Decomposition ~ Treatment * Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  Total_Decomposition ~ Treatment + Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2)

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Total_Decomposition")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}


# 11. Microbial_Decomposition 
fit1 <- stan_glmer(
  Microbial_Decomposition ~ Treatment * Position + (1 | Channel), 
  family = poisson(link = log),
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  Microbial_Decomposition ~ Treatment + Position + (1 | Channel), 
  data = df,
  na.action = "na.omit",
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2)

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = paste0("Exp_1_Community_Ecosystem_",data.name,"_Experimental_Streams_Microbial_Decomposition")
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}






#### Exp_3_Community_Ecosystem_Stream_Mesocosms ####
# Exp_3_Community_Ecosystem_Stream_Mesocosms => Table 1, Extended data Figure 6 and Extended data Table 4
# Response variables: Chironomidae, Primary_Production, Total_Decomposition and Microbial_Decomposition 

## Load dataset
df <- read_csv("data/Exp_3_Community_Ecosystem_Stream_Mesocosms.csv")
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 4)) # SHAM as reference


# 1. Chironomidae
fit1 <- stan_glmer(
  Chironomidae ~ Treatment*Time*Position+(1|Block/Stream_Mesocosm),
  family = poisson(link=log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Chironomidae ~Treatment*Time + Position+(1|Block/Stream_Mesocosm), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  Chironomidae ~Treatment + Time + Position+(1|Block/Stream_Mesocosm), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  Chironomidae ~Treatment*Position+ Time +(1|Block/Stream_Mesocosm), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Stream_Mesocosm)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Stream_Mesocosm)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position +(1|Block/Stream_Mesocosm)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position+Time +(1|Block/Stream_Mesocosm)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

# Select the best model
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Community_Ecosystem_Stream_Mesocosms_Chironomidae"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}



# 2. Primary_Production
fit1 <- stan_glmer(
  Primary_Production ~ Treatment*Time*Position+(1|Block/Stream_Mesocosm),
  family = binomial(link=logit),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Primary_Production ~Treatment*Time + Position+(1|Block/Stream_Mesocosm), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  Primary_Production ~Treatment + Time + Position+(1|Block/Stream_Mesocosm), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  Primary_Production ~Treatment*Position+ Time +(1|Block/Stream_Mesocosm), 
  family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Stream_Mesocosm)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Stream_Mesocosm)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position +(1|Block/Stream_Mesocosm)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position+Time +(1|Block/Stream_Mesocosm)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

# Select the best model
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Community_Ecosystem_Stream_Mesocosms_Primary_Production"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}




# 3. Total_Decomposition
fit1 <- stan_glmer(
  log(Total_Decomposition) ~ Treatment*Time*Position+(1|Block/Stream_Mesocosm),
  #family = binomial(link=logit),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  log(Total_Decomposition) ~Treatment*Time + Position+(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  log(Total_Decomposition) ~Treatment + Time + Position+(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  log(Total_Decomposition) ~Treatment*Position+ Time +(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Stream_Mesocosm)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Stream_Mesocosm)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position +(1|Block/Stream_Mesocosm)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position+Time +(1|Block/Stream_Mesocosm)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

# Select the best model
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Community_Ecosystem_Stream_Mesocosms_Total_Decomposition"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}



# 5. Microbial_Decomposition
fit1 <- stan_glmer(
  log(Microbial_Decomposition) ~ Treatment*Time*Position+(1|Block/Stream_Mesocosm),
  #family = binomial(link=logit),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  log(Microbial_Decomposition) ~Treatment*Time + Position+(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit3 <- stan_glmer(
  log(Microbial_Decomposition) ~Treatment + Time + Position+(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit4 <- stan_glmer(
  log(Microbial_Decomposition) ~Treatment*Position+ Time +(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Stream_Mesocosm)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Stream_Mesocosm)
loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Time+Position +(1|Block/Stream_Mesocosm)
loo4 <- loo(fit4,k_threshold = 0.7) # Treatment*Position+Time +(1|Block/Stream_Mesocosm)
loo <- compare(loo1, loo2, loo3, loo4) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

# Select the best model
looic <- c(loo1$looic, loo2$looic, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Community_Ecosystem_Stream_Mesocosms_Microbial_Decomposition"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}




#### Exp_3_Excretion_Population_Stream_Mesocosms ####
# Exp_3_Excretion_Population_Stream_Mesocosms => Table 1 and Extended data Table 4
# Response variables: N_Excretion and P_Excretion

## Load dataset
df <- read_csv("data/Exp_3_Excretion_Population_Stream_Mesocosms.csv")
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 3)) # SHAM as reference

# 1. Excretion_N
fit1 <- stan_glmer(
  log(Excretion_N) ~ Treatment +(1|Block/Stream_Mesocosm), # with mass as covariate
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## RESULTS
best <- fit1
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Excretion_Population_Stream_Mesocosms_Excretion_N"
  ### SAVE
  save(best,file=paste0("results/",model.name,".Rdata"))
}


# 2. Excretion_P
fit1 <- stan_glmer(
  log(Excretion_P) ~ Treatment +(1|Block/Stream_Mesocosm), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Excretion_Population_Stream_Mesocosms_Excretion_P"
  ### SAVE
  save(best,file=paste0("results/",model.name,".Rdata"))
}





#### Exp_3_Filamentous_Algae_Stream_Mesocosms ####
# Exp_3_Filamentous_Algae_Stream_Mesocosms => Table 1 and Extended data Table 4
#Response variables: Filamentous_Algae

## Load dataset
df <- read_csv("data/Exp_3_Filamentous_Algae_Stream_Mesocosms.csv")
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 4)) # SHAM as reference
# data transformed to logit scale
df$Filamentous_Algae[df$Filamentous_Algae==1] <- .99 # must be <1
df$logit.Fila <- logit(df$Filamentous_Algae)


fit1 <- stan_glmer(
  logit.Fila ~ Treatment*Position+(1|Block/Stream_Mesocosm), 
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


fit2 <- stan_glmer(
  logit.Fila ~Treatment + Position+(1|Block/Stream_Mesocosm), 
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position +(1|Block/Stream_Mesocosm)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position +(1|Block/Stream_Mesocosm)
loo <- compare(loo1, loo2) 
# The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
# If the difference is be positive then the second model is preferred.

# Select the best model
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)


## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Filamentous_Algae_Stream_Mesocosms"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}





#### Exp_3_Nutrient_Fluxes_Stream_Mesocosms ####
# Exp_3_Nutrient_Fluxes_Stream_Mesocosms=> Table 1 and Extended data Table 4
# Response variables: NH4_concentration and PO4_concentration

# Load dataset
df <- read_csv("data/Exp_3_Nutrient_Fluxes_Stream_Mesocosms.csv")
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 4)) # SHAM as reference

# 1. NH4_concentration
fit1 <- stan_glmer(
  log(NH4_Concentration) ~Treatment+Time+(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  log(NH4_Concentration) ~Treatment*Time +(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Stream_Mesocosm)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Stream_Mesocosm)
loo <- compare(loo1, loo2)#, loo3, loo4) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)#, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Nutrient_Fluxes_Stream_Mesocosms_NH4_concentration"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}



# 2. PO4_Concentration
fit1 <- stan_glmer(
  log(PO4_Concentration) ~Treatment+Time+(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  log(PO4_Concentration) ~Treatment*Time +(1|Block/Stream_Mesocosm), 
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Stream_Mesocosm)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Stream_Mesocosm)
loo <- compare(loo1, loo2)#, loo3, loo4) 

# Select the best model
looic <- c(loo1$looic, loo2$looic)#, loo3$looic, loo4$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Exp_3_Nutrient_Fluxes_Stream_Mesocosms_PO4_concentration"
  ### SAVE
  save(best,loo, file=paste0("results/",model.name,".Rdata"))
}
