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
# Load dataset:
library(readxl)
data <- read_excel("data/Data_for_Path_2015.xlsx")

# Cange format of data:
#data <- as.data.frame(data)
mode(data$Chlo_Torch_All) <-"numeric" # transform vector of characters to numeric
#data <- data[order(data$Treatment),]
data$Treatment <- as.factor(data$Treatment)
data$Position <- as.factor(data$Position)
data <- as.data.frame(data)
#data$Position <- as.numeric(as.factor(data$Position))



# Data type to use:
data.name = "nofish" # "fish" vs "nofish" # Run analysis for fish or no fish
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





###______________________________________________###
###_________________ANALYSIS_________________________###
###______________________________________________###


############ 1. Total_Prim_Cons (somme de tous les consommateurs )
fit1 <- stan_glmer(
  Total_Prim_Cons ~ Treatment*Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
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
model.name = "Total_Prim_Cons"
Total_Prim_Cons <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}




########### 2. Total_Sec_Cons (somme de tous les prédateurs (Rhyaco + Polycentrodidae)
# fit1<- stan_glmer(
#   Total_Sec_Cons ~ Treatment*Position + (1 | Channel),
#   family = poisson(link = log),
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# fit2<- stan_glmer(
#   Total_Sec_Cons ~ Treatment + Position + (1 | Channel),
#   family = poisson(link = log),
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
# loo <- compare(loo1, loo2) 
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
# model.name = "Total_Sec_Cons"
# Total_Sec_Cons <- best
# ### SAVE
# save(Total_Sec_Cons,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
# }




######### 3. Pred_Prey_Ratio 
# # ratio entre Total Sec Cons et Total Prim Cons
# #df$Pred_Prey_Ratio <- log(df$Pred_Prey_Ratio) # log-transformed
# 
# ## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
# fit1<- stan_glmer(
#   log(Pred_Prey_Ratio) ~ Treatment * Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# fit2<- stan_glmer(
#   log(Pred_Prey_Ratio) ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
# loo <- compare(loo1, loo2) 
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
# model.name = "Pred_Prey_Ratio"
# Pred_Prey_Ratio <- best
# ### SAVE
# save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
# }


###______________________________________________###
###_____________ PRIMARY CONSUMERS______________ ###
###______________________________________________###


########### 3. Rhyacophilidae
fit1<- stan_glmer(
  Rhyacophilidae ~ Treatment * Position + (1 | Channel),
  family = poisson(link = log),
  data = df,
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
model.name = "Rhyacophilidae"
Rhyacophilidae <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}




###______________________________________________###
###____________ SECONDARY CONSUMERS_____________ ###
###______________________________________________###


## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Chironomidae ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Chironomidae ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
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
model.name = "Chironomidae"
Chironomidae <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}





## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Polycentropodidae ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Polycentropodidae ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
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
model.name = "Polycentropodidae"
Polycentropodidae <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}




## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Hydropsychidae ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Hydropsychidae ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
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
model.name = "Hydropsychidae"
Hydropsychidae <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}



## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Baetidae ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Baetidae ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

# fit3 <- stan_glmer(
#   Baetidae ~ Treatment:Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
#   family = poisson(link = log),
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
# )


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
model.name = "Baetidae"
Baetidae <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}




## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Simuliidae ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
)

fit2 <- stan_glmer(
  Simuliidae ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

# fit3 <- stan_glmer(
#   Simuliidae ~ Treatment:Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
#   family = poisson(link = log),
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
# )


## Compare models
loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
#loo3 <- loo(fit3,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
loo <- compare(loo1, loo2)#, loo3) 

# Select the best model
# elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# best <- get(paste0("fit",which.max(elpd_loo)))
looic <- c(loo1$looic, loo2$looic)#, loo3$looic)
best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
model.name = "Simuliidae"
Simuliidae <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}


## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Planorbidae ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  family = poisson(link = log),
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  Planorbidae ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
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
model.name = "Planorbidae"
Planorbidae <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}




###______________________________________________###
###____________ Perc_Contr_Inv    _____________ ###
###______________________________________________###
# contribution des Inv au processus de decompo

# ## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
# fit1<- stan_glmer(
#   Perc_Contr_Inv ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# fit2 <- stan_glmer(
#   Perc_Contr_Inv ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta),
#   QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Position + (1 | Channel)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment+Position + (1 | Channel)
# loo <- compare(loo1, loo2) 
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
# model.name = "Perc_Contr_Inv"
# Perc_Contr_Inv <- best
# ### SAVE
# save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
# }


###______________________________________________###
###______________ PRIMARY PRODUCTION ____________###
###______________________________________________###

df$Chlo_Torch_All <- log(df$Chlo_Torch_All) # log-transformed


## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  Chlo_Torch_All ~ Treatment*Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit2 <- stan_glmer(
  Chlo_Torch_All ~ Treatment+Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
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
model.name = "Chlo_Torch_All"
Chlo_Torch_All <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}


###______________________________________________###
###__________GLOBAL DECOMPOSITION________________###
###______________________________________________###

#df$K_CM <- sqrt(df$K_CM) # log-transformed

## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  K_CM ~ Treatment * Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit2 <- stan_glmer(
  K_CM ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
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
model.name = "Global_decompo"
K_CM <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}




###______________________________________________###
###__________MICROBIAL DECOMPOSITION_____________###
###______________________________________________###

#df$K_FM <- sqrt(df$K_FM) # log-transformed

## Fit the model with treatment (SHAM vs GH) predictor AND position interaction
fit1 <- stan_glmer(
  K_FM ~ Treatment * Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta),
  QR = TRUE
)

fit2 <- stan_glmer(
  K_FM ~ Treatment + Position + (1 | Channel), # y = A*B : Multiplicative fully factorial analysis of variance of y against A and B
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
model.name = "Micro_decompo"
K_FM <- best
### SAVE
save(best,loo, file=paste0("results/IMS2015_",model.name,"_",data.name,".Rdata"))
}


### SAVE ALL
#save(Total_Prim_Cons, Total_Sec_Cons,Pred_Prey_Ratio,Rhyacophilidae,Chironomidae,Polycentropodidae,Hydropsychidae,Baetidae,Simuliidae,Planorbidae,Chlo_Torch_All, Perc_Contr_Inv, K_CM, K_FM, file=paste0("results/IMS2015_",data.name,".Rdata"))






###______________________________________________###
###______________OUTPUTS_______ __________________###
###______________________________________________###

# model.name <- c("Total_Prim_Cons", "Total_Sec_Cons","Pred_Prey_Ratio","Rhyacophilidae","Chironomidae","Polycentropodidae","Hydropsychidae","Baetidae"
#                 ,"Simuliidae","Planorbidae","Chlo_Torch_All", "Perc_Contr_Inv", "K_CM", "K_FM")
# 
# #--------------- FUNCTIONS ---------------#
# gf <- function(x){if(mean(x)>=0){mean(x>=0)}else{mean(x<0)}} # function to calculate the proportion of posterior with same sign as mean
# 
# plot.effect <- function(model, main, formula){
#   fit <- as.array(model$stanfit)
#   TreatmentSHAM <- as.vector(fit[,1:3,"(Intercept)"])
#   TreatmentGHA <- as.vector(fit[,1:3,"TreatmentGH"])
#   
#   if( paste0(best$formula)[3] == "Treatment * Position + (1 | Channel)"){
#   TreatmentGHB <- TreatmentGHA + as.vector(fit[,1:3,"TreatmentGH:PositionB"])
#   TreatmentGHC <- TreatmentGHA + as.vector(fit[,1:3,"TreatmentGH:PositionC"])
#   TreatmentGHD <- TreatmentGHA + as.vector(fit[,1:3,"TreatmentGH:PositionD"])
#   TreatmentGHE <- TreatmentGHA + as.vector(fit[,1:3,"TreatmentGH:PositionE"])
#   
# deltas <- cbind(
#   TreatmentGHA=TreatmentGHA
#   , TreatmentGHB=TreatmentGHB
#   , TreatmentGHC=TreatmentGHC
#   , TreatmentGHD=TreatmentGHD
#   , TreatmentGHE=TreatmentGHE
# )
# 
# q <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),2) # quantiles
# p <- round(apply(deltas, 2, function(x) gf(x)),2)
# res <- t(rbind(q,p))
#   }
#   
#   
#   res <- quantile(deltas, probs=c(.025, .25, .5, .75, .975))
#   
#   plot(NULL, xlim=c(0,2),ylim=c(min(res), max(res)),xaxt="n",ylab=paste("Effect of GH (",expression(beta),")",sep=""),main= main,bty="n") #range(mcmc)
#   abline(h=0,lty=2)
#   points(1,res["50%"], pch=16,cex=1.5)
#   segments(1,res["2.5%"],1,res["97.5%"])
#   segments(1,res["25%"],1,res["75%"],lwd=3)
#   
#   text(1,1.5, paste0(formula),cex=0.75)
#   text(1,1.4,  (paste0("beta = ",round(median(mcmc),2)," [",round(quantile(mcmc, probs=0.025),2),";",round(quantile(mcmc, probs=0.975),2),"]; P>0 = ",round(mean(mcmc>0),3))),cex=0.75)
#   text(1,1.3,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
#   text(1,1.25,"(i.e., do not overlap with 0)",cex=.4)
#   
# }
# 
# plot.pred <- function(model, ylab, xlab, main, fun){
#   
#   # Precictive values
#   nd <- data.frame(Treatment = c(rep("SHAM",5),rep("GH",5)), Position = rep(c("A", "B", "C", "D", "E"),2))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
#   y_rep <- posterior_predict(model,newdata =nd, fun = fun , re.form=NA)
#   #boxplot(y_rep,col=rep(1:2,5))
#   y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
#   
#   y.pred <- ifelse(y.pred =="-Inf",0, y.pred) 
#   
#   #  formula <- "Treatment * Position + (1 | Channel)"
#   #  if (paste0(best$formula)[3] == formula) {
#   gap <- .05
#   plot(NULL,xlim=c(1,5),ylim=c(min(y.pred),max(y.pred)),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
#   axis(1, labels=LETTERS[1:5], at = 1:5)
#   points((1:5)-gap,y.pred["50%",1:5],pch=16,col=1, cex=1.5)
#   points((1:5)+gap,y.pred["50%",6:10],pch=16,col=2, cex = 1.5)
#   segments((1:5)-gap,y.pred["25%",1:5],(1:5)-gap,y.pred["75%",1:5],lwd=3,col=1)
#   segments((1:5)+gap,y.pred["25%",6:10],(1:5)+gap,y.pred["75%",6:10],lwd=3,col=2)
#   segments((1:5)-gap,y.pred["2.5%",1:5],(1:5)-gap,y.pred["97.5%",1:5],lwd=1,col=1)
#   segments((1:5)+gap,y.pred["2.5%",6:10],(1:5)+gap,y.pred["97.5%",6:10],lwd=1,col=2)
#   legend("topright",legend=c("SHAM", "GH"),col=1:2,pch=c(16,16),bty="n")
#   #points(as.numeric(as.factor(model$data$Position)), log(model$data$model+.1), col=as.numeric(as.factor(model$data$Treatment)), pch=as.numeric(as.factor(model$data$Channel)))
#   # }
  
  # formula <- "Treatment + Position + (1 | Channel)"
  # if (paste0(best$formula)[3] == formula) {
  #   gap <- .05
  #   plot(NULL,xlim=c(0.5,1.5),ylim=c(min(y.pred),max(y.pred)),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
  #   axis(1, labels=LETTERS[1], at = 1)
  #   points(1-gap,y.pred["50%",1],pch=16,col=1, cex=2)
  #   points(1+gap,y.pred["50%",6],pch=16,col=2, cex = 2)
  #   segments(1-gap,y.pred["25%",1],1-gap,y.pred["75%",1],lwd=3,col=1)
  #   segments(1+gap,y.pred["25%",6],1+gap,y.pred["75%",6],lwd=3,col=2)
  #   segments(1-gap,y.pred["2.5%",1],1-gap,y.pred["97.5%",1],lwd=1,col=1)
  #   segments(1+gap,y.pred["2.5%",6],1+gap,y.pred["97.5%",6],lwd=1,col=2)
  #   legend("topright",legend=c("SHAM", "GH"),col=1:2,pch=c(16,16),bty="n")
  #   #points(as.numeric(as.factor(model$data$Position)), log(model$data$model+.1), col=as.numeric(as.factor(model$data$Treatment)), pch=as.numeric(as.factor(model$data$Channel)))
  # }
}
# 
# plot.res <- function(model, ylab, xlab, main){
#   gap <- .05
#   plot(NULL,xlim=c(0,5),ylim=range(model$stan_summary[3:10,c("2.5%","97.5%")]),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
#   axis(1, labels=LETTERS[2:5], at = 1:4)
#   abline(h=0,lty=2)
#   points((1:4)-gap,model$stan_summary[3:6,"50%"],pch=16,col=rep(1,4), cex=1.5)
#   points((1:4)+gap,model$stan_summary[7:10,"50%"],pch=16,col=rep(2,4), cex=1.5)
#   segments((1:4)-gap,model$stan_summary[3:6,"2.5%"],(1:4)-gap,model$stan_summary[3:6,"97.5%"],lwd=1,col=rep(1,4))
#   segments((1:4)-gap,model$stan_summary[3:6,"25%"],(1:4)-gap,model$stan_summary[3:6,"75%"],lwd=3,col=rep(1,4))
#   segments((1:4)+gap,model$stan_summary[7:10,"2.5%"],(1:4)+gap,model$stan_summary[7:10,"97.5%"],lwd=1,col=rep(2,4))
#   segments((1:4)+gap,model$stan_summary[7:10,"25%"],(1:4)+gap,model$stan_summary[7:10,"75%"],lwd=3,col=rep(2,4))
# } # plot of position effects


# Table
#table <- array(,dim=c(14,4));colnames(table)<-c("levels","best model","Effect of GH", "proportion of the posterior with the same sign as the mean")
## Table
# mcmc <- as.array(best$stanfit)
# TreatmentGH <- as.vector(mcmc[,1:3,"TreatmentGH"])
# table[14,1]<- "K_FM"
# table[14,2]<- paste0(best$formula)[3]
# table[14,3]<-  paste0(round(median(TreatmentGH),3),"[",round(quantile(TreatmentGH, probs=0.025),3),";",round(quantile(TreatmentGH, probs=0.975),3),"]")
# table[14,4]<-  paste0(round(gf(TreatmentGH)*100,1),"%")
# write.csv(table,file=paste0("results/Synthesis-",data.name,".csv"),row.names = FALSE)








## SUMMARY
# sink(paste0("results/IMS2015_summary_",data.name,".txt"))
# for (myvar in model.name) {
#   
#   ## LOAD MODEL
#   best = loo = mcmc = NULL
#   load(paste0("results/IMS2015_",myvar,"_",data.name,".Rdata"))
#   #best <- get(myvar)
#   mcmc <- as.array(best$stanfit)
#   
# 
#   
#   ## POSTERIOR CHECK
#   y <- best$y
#   yrep <- posterior_predict(best)
#   group <- best$data$Treatment
#   
#   cat("\n__________",myvar,"__________\n")
#   
# 
#   
#     cat("\n Model selection :\n")
#     cat( "MODEL 1:  Treatment*Position + (1 | Channel) \n")
#     cat( "MODEL 2:  Treatment+Position + (1 | Channel) \n")
#     cat("\n Loo :\n")
#     print(loo)
#     cat( "The difference in ELPD (model1 - model2) will be negative if the expected out-of-sample predictive accuracy of the first model is higher (as AIC, lower is better). If the difference is be positive then the second model is preferred.\n")
#     # Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
#     # That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.
#     
#   #launch_shinystan(model)
#   cat("\nBest model:\n")
#   print(best$formula)
#   print(best$family)
#   
#   ## Print summary
#   print(round(best$stan_summary[,c("2.5%", "50%","97.5%","Rhat")],2))
#   
#   #cat("\n____________________\n")
#   cat("\n Posterior distributions for the effect of Treatments (median [95% CI]):\n")
#   #mcmc <- as.array(best$stanfit)
# 
#   TreatmentSHAM <- as.vector(mcmc[,1:3,"(Intercept)"])
#   TreatmentGHA <- as.vector(mcmc[,1:3,"TreatmentGH"])
#   
#   if( paste0(best$formula)[3] == "Treatment * Position + (1 | Channel)"){
#     TreatmentGHB <- TreatmentGHA + as.vector(mcmc[,1:3,"TreatmentGH:PositionB"])
#     TreatmentGHC <- TreatmentGHA + as.vector(mcmc[,1:3,"TreatmentGH:PositionC"])
#     TreatmentGHD <- TreatmentGHA + as.vector(mcmc[,1:3,"TreatmentGH:PositionD"])
#     TreatmentGHE <- TreatmentGHA + as.vector(mcmc[,1:3,"TreatmentGH:PositionE"])
#     
#     deltas <- cbind(
#       TreatmentGHA=TreatmentGHA
#       , TreatmentGHB=TreatmentGHB
#       , TreatmentGHC=TreatmentGHC
#       , TreatmentGHD=TreatmentGHD
#       , TreatmentGHE=TreatmentGHE
#     )
#     
#     q <- round(apply(deltas, 2, quantile , probs=c(.025, .25, .5, .75, .975)),2) # quantiles
#     p <- round(apply(deltas, 2, function(x) gf(x)),2)
#     res <- t(rbind(q,p))
#   } else {
#     
#     quantile(TreatmentGHA, probs=c(.025, .25, .5, .75, .975),2) # quantiles
#     gf(TreatmentGHA)
#   }
#   
# 
#     
#     # posterior distributions
#     #deltas <- cbind(GH = as.vector(TreatmentGH))
# 
#   # res <- quantile(TreatmentGH , probs=c(.025, .25, .5, .75, .975), na.rm=TRUE) # quantiles
#   # print(paste0(round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),2)
#   #              , " [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),2),"]"
#   #              , ", P[>0] = "
#   #              , round(mean(TreatmentGH>0),3) # P: confidence that the parameter is positive (P>0)))
#   # ))
#   #p <- round(mean(TreatmentGH>0),3) # P: confidence that the parameter is positive (P>0)
#   #mytable<- c(res,"P>0"=p)
#   #print(mytable)
#   #cat("P>0: Proportion of positive posterior values, i.e. confidence that the effect is positive \n")
#   # TreatmentGH <- as.vector(mcmc[,1:CORES,"TreatmentGH"])
#   # cat("Effect of GH: ", median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]\n")
#   # cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentGH)*100,1),"%\n")
# 
#   
#   cat( "\n Posterior check: comparison between observed values (y) and predicted values (yrep) \n")
#   cat("Ratio (y/yrep) should include the value 1\n")
#   # ratio y vs yrep
#   r=NULL;for (i in 1:nrow(yrep)){r[i]<-mean(y/yrep[i,],na.rm=TRUE)}
#   cat("Ratio (y/yrep): ", median(r),"[",quantile(r, probs=0.025),";",quantile(r, probs=0.975),"]\n")
# } #End loop 
# sink()









### CREATE TABLE
model.name <- c("Rhyacophilidae"
                ,"Polycentropodidae"
                ,"Total_Prim_Cons"
                ,"Chironomidae"
                ,"Hydropsychidae"
                ,"Baetidae"
                ,"Simuliidae"
                ,"Planorbidae"
                ,"Chlo_Torch_All"
                ,"Global_decompo"
                , "Micro_decompo"
                )



sink(paste0("results/IMS2015_table2_",data.name,".csv"))
for (myvar in model.name) {
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2015_",myvar,"_",data.name,".Rdata"))
  
  
  if( paste0(best$formula)[3] == "Treatment + Position + (1 | Channel)"){
  mcmc <- as.array(best$stanfit)
  TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  
tmp <-  paste0(
          myvar
          , ","
          , "-"
          , ","
         , best$formula
         , ","
         , round(quantile(TreatmentGH , probs=c(.5), na.rm=TRUE),3)
         , " [", round(quantile(TreatmentGH , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGH , probs=c(.975), na.rm=TRUE),3),"]"
         , ", "
          , round(gf(TreatmentGH),3) # P: confidence that the parameter is positive (P>0)
  )

print(tmp[3])
} # end if


if( paste0(best$formula)[3] == "Treatment * Position + (1 | Channel)"){
  
  mcmc <- as.array(best$stanfit)
  TreatmentGHA <- mcmc[,1:CORES,"TreatmentGH"]
  TreatmentGHB <- TreatmentGHA + mcmc[,1:CORES,"TreatmentGH:PositionB"]
  TreatmentGHC <- TreatmentGHA + mcmc[,1:CORES,"TreatmentGH:PositionC"]
  TreatmentGHD <- TreatmentGHA + mcmc[,1:CORES,"TreatmentGH:PositionD"]
  TreatmentGHE <- TreatmentGHA + mcmc[,1:CORES,"TreatmentGH:PositionE"]
  
  tmp <-  paste0(
    myvar
    , ","
    , "A"
    , ","
    , best$formula
    , ","
    , round(quantile(TreatmentGHA , probs=c(.5), na.rm=TRUE),3)
    , " [", round(quantile(TreatmentGHA , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGHA , probs=c(.975), na.rm=TRUE),3),"]"
    , ", "
    , round(gf(TreatmentGHA),4) # P: confidence that the parameter is positive (P>0)
  )
  print(tmp[3])
  
  tmp <-  paste0(
    myvar
    , ","
    , "B"
    , ","
    , best$formula
    , ","
    , round(quantile(TreatmentGHB , probs=c(.5), na.rm=TRUE),3)
    , " [", round(quantile(TreatmentGHB , probs=c(.025), na.rm=TRUE),3), "; ", round(quantile(TreatmentGHB , probs=c(.975), na.rm=TRUE),3),"]"
    , ", "
    , round(gf(TreatmentGHB),4) # P: confidence that the parameter is positive (P>0)
  )
  print(tmp[3])
  
  tmp <-  paste0(
    myvar
    , ","
    , "C"
    , ","
    , best$formula
    , ","
    , round(quantile(TreatmentGHC , probs=c(.5), na.rm=TRUE),2)
    , " [", round(quantile(TreatmentGHC , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGHC , probs=c(.975), na.rm=TRUE),3),"]"
    , ", "
    , round(gf(TreatmentGHC),4) # P: confidence that the parameter is positive (P>0)
  )
  print(tmp[3])
  
  tmp <-  paste0(
    myvar
    , ","
    , "D"
    , ","
    , best$formula
    , ","
    , round(quantile(TreatmentGHD , probs=c(.5), na.rm=TRUE),2)
    , " [", round(quantile(TreatmentGHD , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGHD , probs=c(.975), na.rm=TRUE),3),"]"
    , ", "
    , round(gf(TreatmentGHD),4) # P: confidence that the parameter is positive (P>0)
  )
  print(tmp[3])
  
  tmp <-  paste0(
    myvar
    , ","
    , "E"
    , ","
    , best$formula
    , ","
    , round(quantile(TreatmentGHE , probs=c(.5), na.rm=TRUE),2)
    , " [", round(quantile(TreatmentGHE , probs=c(.025), na.rm=TRUE),2), "; ", round(quantile(TreatmentGHE , probs=c(.975), na.rm=TRUE),3),"]"
    , ", "
    , round(gf(TreatmentGHE),4) # P: confidence that the parameter is positive (P>0)
  )
  print(tmp[3])
} # end if
  
  
} # end loop
sink()














# ## FIGURES EFFETCS
# pdf(paste("results/IMS2015_",data.name,".pdf",sep=""),width=8, height=10,)
# par(mfrow=c(3,3))
# par(mar=c(2.5, 4.5, 2.5, 1.5))
# for (myvar in model.name) {
# 
#   ## LOAD MODEL
#   best = loo = mcmc = NULL
#   load(paste0("results/IMS2015_",myvar,"_",data.name,".Rdata"))
#   #best <- get(myvar)
#   mcmc <- as.array(best$stanfit)
#   
#   ## POSTERIOR CHECK
#   y <- best$y
#   yrep <- posterior_predict(best)
#   group <- best$data$Treatment
#   
# 
# main = myvar
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab=myvar,xlab="Position",main=main,fun=NULL)
# if(myvar != "Chlo_Torch_All"){
#   points(as.numeric(as.factor(best$data$Position)), best$y, col=as.numeric(as.factor(best$data$Treatment)), pch=as.numeric(as.factor(best$data$Channel)))
# }
# 
# } # end loop myvar
# dev.off()






# ## FIGURES POSTERIOR CHECK
# for (myvar in model.name) {
#   #myvar= "Sum_SRP"
#   pdf(paste0("results/IMS2015_",myvar,"_check_",data.name,".pdf"),width=8, height = 8)
#   
#   ## LOAD MODEL
#   best = loo = mcmc = NULL
#   load(paste0("results/IMS2015_",myvar,"_",data.name,".Rdata"))
#   #best <- get(myvar)
#   mcmc <- as.array(best$stanfit)
#   
#   ## POSTERIOR CHECK
#   y <- best$y
#   yrep <- posterior_predict(best)
#   group <- best$data$Treatment
#   
#   
#   #par(mfrow=c(2,2))
#   #par(mar=c(3.5, 4.5, 4.5, 1.5))
#   # # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
#   plot1 <- ppc_stat(y, yrep)
#   plot2 <- ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
#   #color_scheme_set("teal")
#   #ppc_stat_grouped(y, yrep, group)
#   plot3 <- pp_check(y, yrep, ppc_dens_overlay)
#   # pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
#   # # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
#   plot4 <- ppc_scatter_avg(y, yrep)
#   # ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)
#   
#   grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
#   dev.off()
# } 





# ## Total_Prim_Cons
# load(paste0("results/Total_Prim_Cons-",data.name,".Rdata"))
# main = "Total_Prim_Cons" #"TOTAL PRIM CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Total_Prim_Cons (log)",xlab="Position",main=main,fun=log)
# points(as.numeric(as.factor(df$Position)), log(df$Total_Prim_Cons+.01), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Total_Sec_Cons
# load(paste0("results/Total_Sec_Cons-",data.name,".Rdata"))
# main =  "Total_Sec_Cons" #"TOTAL SEC CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Total_Sec_Cons (log)",xlab="Position",main=main,fun=log)
# points(as.numeric(as.factor(df$Position)), log(df$Total_Sec_Cons+.01), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# 
# ## Pred_Prey_Ratio
# load(paste0("results/Pred_Prey_Ratio-",data.name,".Rdata"))
# main = "Pred_Prey_Ratio"   #"RATIO SEC / PRIM CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Pred_Prey_Ratio (log)",xlab="Position",main=main, fun = NULL)
# points(as.numeric(as.factor(df$Position)), log(df$Pred_Prey_Ratio), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# 
## Rhyacophilidae
# load(paste0("results/Rhyacophilidae-",data.name,".Rdata"))
# main = "Rhyacophilidae" #"PRIMARY CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Rhyacophilidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Rhyacophilidae), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))

# 
# ## Chironomidae
# load(paste0("results/Chironomidae-",data.name,".Rdata"))
# main = "Chironomidae" #"SECONDARY CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Chironomidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Chironomidae), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Polycentropodidae
# load(paste0("results/Polycentropodidae-",data.name,".Rdata"))
# main = "Polycentropodidae" #"SECONDARY CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Polycentropodidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Polycentropodidae), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Hydropsychidae
# load(paste0("results/Hydropsychidae-",data.name,".Rdata"))
# main = "Hydropsychidae" #"SECONDARY CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Hydropsychidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Hydropsychidae), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Baetidae
# load(paste0("results/Baetidae-",data.name,".Rdata"))
# main =  "Baetidae" #"SECONDARY CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Baetidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Baetidae), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Simuliidae
# load(paste0("results/Simuliidae-",data.name,".Rdata"))
# main = "Simuliidae"   #"SECONDARY CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Simuliidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Simuliidae), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Planorbidae
# load(paste0("results/Planorbidae-",data.name,".Rdata"))
# main = "Planorbidae" #"SECONDARY CONSUMERS"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Planorbidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Planorbidae), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Perc_Contr_Inv
# load(paste0("results/Perc_Contr_Inv-",data.name,".Rdata"))
# main = "Perc_Contr_Inv" #"INVERTEBRATES DECOMPO"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Perc_Contr_Inv",xlab="Position",main=main, fun = NULL)
# points(as.numeric(as.factor(df$Position)), (df$Perc_Contr_Inv), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# ## Chlo_Torch_All
# load(paste0("results/Chlo_Torch_All-",data.name,".Rdata"))
# main = "Chlo_Torch_All" #"PRIMARY PRODUCTION"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="Chlo_Torch_All",xlab="Position",main=main, fun = exp)
# points(as.numeric(as.factor(df$Position)), (df$Chlo_Torch_All), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# 
# ## K_CM
# load(paste0("results/K_CM-",data.name,".Rdata"))
# main = "K_CM" #"GLOBAL DECOMPOSITION"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="K_CM",xlab="Position",main=main, fun = NULL)
# points(as.numeric(as.factor(df$Position)), (df$K_CM), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# 
# ## K_FM
# load(paste0("results/K_FM-",data.name,".Rdata"))
# main = "K_FM" #"MICRO. DECOMPOSITION"
# plot.effect(best, main=main, paste0(best$formula)[3])
# plot.res(best, ylab="Position effect",xlab="",main=main)
# plot.pred(best,ylab="K_FM",xlab="Position",main=main, fun = NULL)
# points(as.numeric(as.factor(df$Position)), (df$K_FM), col=as.numeric(as.factor(df$Treatment)), pch=as.numeric(as.factor(df$Channel)))
# 
# 
# dev.off()

