rm(list=ls())   # Clear memory

#--------------- PACKAGES ---------------#
library(readxl)
# see https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library("rstanarm")
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(betareg)
library(bayesplot)
library(gridExtra)

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
gf <- function(x){if(mean(x)>=0){mean(x>=0)}else{mean(x<0)}} # function to calculate the proportion of posterior with same sign as mean


#--------------- DIRECTORY --------------#
setwd("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/")
#setwd("~/Documents/mbuoro/Ims/")

#-----------------ANALYSIS---------------#
CHAINS  = 3 # number of chains
CORES = CHAINS # number of core
ITER = 10000 # number of iterations
WARM = 2000 # warmup
ndelta = 0.99









library(readxl)
data <- read_excel("data/Data_for_Path_2015_Last.xlsx")
View(data)


# Cange format of data:

data$Treatment <- as.factor(data$Treatment)
data$Position <- as.factor(data$Position)
data$Perc_Herb_Chir <- as.numeric(data$Perc_Herb_Chir)
data$Perc_Pred_Chir <- as.numeric(data$Perc_Pred_Chir)
data$Herb_Chir <- as.numeric(data$Herb_Chir)
data$Pred_Chir <- as.numeric(data$Pred_Chir)
data <- as.data.frame(data)
#data$Position <- as.numeric(as.factor(data$Position))



# Data type to use:
data.name = "fish" # "fish" vs "nofish" # Run analysis for fish or no fish
# Datasets
if (data.name == "fish"){
  data.new <- subset(data, Treatment == "GH" | Treatment == "SHAM")
  df <- droplevels(data.new) # cleaning levels
  df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference
} else {
  data.new <- subset(data, Treatment == "GH_NOFISH" | Treatment == "SHAM_NOFISH")
  df <- droplevels(data.new) # cleaning levels
  df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference
  
  # Rename by name: change "GH_NOFISH" to "GH" and "SHAM_NOFISH" to "SHAM"
  levels(df$Treatment)[levels(df$Treatment)=="GH_NOFISH"] <- "GH"
  levels(df$Treatment)[levels(df$Treatment)=="SHAM_NOFISH"] <- "SHAM"
}



df$Perc_Pred_Chir.logit <- logit(df$Perc_Pred_Chir/100)


#########################
###### PREDATORS ########
#########################

## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Perc_Pred_Chir.logit ~ Treatment * Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Perc_Pred_Chir.logit ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  data = df,
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
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Perc_Pred_Chir"
  K_FM <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}



mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
res <- round(quantile(TreatmentGH, probs=c(.025, .25, .5, .75, .975)),3) # quantiles
p <- round(mean(TreatmentGH>0),3) # P: confidence that the parameter is positive (P>0)











## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Pred_Chir ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Pred_Chir ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
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
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Pred_Chir"
  Polycentropodidae <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}


mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
res <- round(quantile(TreatmentGH, probs=c(.025, .25, .5, .75, .975)),3) # quantiles
p <- round(mean(TreatmentGH>0),3) # P: confidence that the parameter is positive (P>0)



#########################
###### HERBIVORES #######
#########################

## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Perc_Herb_Chir ~ Treatment * Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Perc_Herb_Chir ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  data = df,
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
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Perc_Herb_Chir"
  K_FM <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}


mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
res <- round(quantile(TreatmentGH, probs=c(.025, .25, .5, .75, .975)),3) # quantiles
p <- round(mean(TreatmentGH>0),3) # P: confidence that the parameter is positive (P>0)






## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Herb_Chir ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Herb_Chir ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
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
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Herb_Chir"
  Polycentropodidae <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}


mcmc <- as.array(best$stanfit)
TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
res <- round(quantile(TreatmentGH, probs=c(.025, .25, .5, .75, .975)),3) # quantiles
p <- round(mean(TreatmentGH>0),3) # P: confidence that the parameter is positive (P>0)
