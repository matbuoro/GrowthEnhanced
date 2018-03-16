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
#if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Nb_Chironomidae"
  Nb_Chironomidae <- best
  ### SAVE
  save(best,loo, file=paste0("results/IMS2016_",model.name,".Rdata"))
#}

