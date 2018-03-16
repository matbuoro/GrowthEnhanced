rm(list=ls())   # Clear memory

#--------------- PACKAGES ---------------#
#install.packages("devtools")
library(readxl)
# see https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library("rstanarm")
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(betareg)
library(bayesplot)
library(gridExtra)
#devtools::install_github("rasmusab/bayesian_first_aid")
library(BayesianFirstAid)

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

### Calcul du % de valeurs positives #######
prop<-function(x){
  # Calcul des % des paramËtres
  x.pos<-x[x>0]    # valeurs > ‡ 0
  x.neg<-x[x<0]    # valeurs < ‡ 0
  # % de valeurs > ‡ 0
  prop.x.neg<-length(x.neg)/length(x)  
  prop.x.pos<-length(x.pos)/length(x) 
  return(prop.x.neg,prop.x.pos)
}

gf <- function(x){if(mean(x)>=0){mean(x>=0)}else{mean(x<0)}} # function to calculate the proportion of posterior with same sign as mean

K=function(W,L) {(1e5 * W)/(L^3)}

#--------------- DIRECTORY --------------#
setwd("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/")
#setwd("~/Documents/mbuoro/Ims/")

#-----------------ANALYSIS---------------#
CHAINS  = 3 # number of chains
CORES = CHAINS # number of core
ITER = 30000 # number of iterations
WARM = 5000 # warmup
ndelta = 0.99








###______________________________________________###
###_________________INDOOR_______________________###
###______________________________________________###

Growth_Mass_Length_Ims_2015_all <- read_excel("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/data/Growth_Mass_Length_Ims_2015_all.xlsx", sheet = "Indoor")
df <- as.data.frame(Growth_Mass_Length_Ims_2015_all)
df$Treatment <- as.factor(df$Treatment)
levels(df$Treatment)[levels(df$Treatment)=="sham"] <- "SHAM"
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference


fit <- lm(log(df$Initial_Mass)~log(df$Initial_FL))
df$Ki <- residuals.lm(fit)#K(df$Initial_Mass , df$Initial_FL)
fit <- lm(log(df$Final_Mass)~log(df$Final_Length))
df$Kf <- residuals.lm(fit)#K(df$Final_Mass , df$Final_Length)

## INITIAL MASS & LENGTH
par(mfrow=c(1,3), byrow = TRUE)

boxplot(df$Initial_Mass~df$Treatment, main="Mass at tagging")
t.test(df$Initial_Mass[df$Treatment=="SHAM"],df$Initial_Mass[df$Treatment=="GH"], paired = TRUE)
#fit <- bayes.t.test(df$Initial_Mass[df$Treatment=="SHAM"],df$Initial_Mass[df$Treatment=="GH"], paired = TRUE)
#plot(fit)

boxplot(df$Initial_FL~df$Treatment, main = "Length at tagging")
t.test(df$Initial_FL[df$Treatment=="SHAM"],df$Initial_FL[df$Treatment=="GH"], paired = TRUE)
#fit <- bayes.t.test(df$Initial_FL[df$Treatment=="SHAM"],df$Initial_FL[df$Treatment=="GH"], paired = TRUE)
#plot(fit)


boxplot(df$Ki~df$Treatment, main = "Condition at tagging")
t.test(df$Ki[df$Treatment=="SHAM"],df$Ki[df$Treatment=="GH"], paired = TRUE)
#fit <- bayes.t.test(df$Ki[df$Treatment=="SHAM"],df$Ki[df$Treatment=="GH"], paired = TRUE)
#plot(fit)

############ MASS / GROWTH #############
# Mass_Growth_Indoor_2015 <- read_excel("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/data/Mass_Growth_Indoor_2015.xlsx", sheet = "Indoor 2015")
# df <- as.data.frame(Mass_Growth_Indoor_2015)
# df$Treatment <- as.factor(df$Treatment)
# levels(df$Treatment)[levels(df$Treatment)=="sham"] <- "SHAM"
# df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference



# fit1=best=NULL
# fit1 <- stan_glmer(
#   Initial_Mass ~ Treatment  + (1| Tank_ID  ),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# best <- fit1
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Ini_Mass_Indoor"
#   ### SAVE
#   save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
# }



fit1=best=NULL
fit1 <- stan_glmer(
  Final_Mass ~ Treatment  + (1| Tank_ID  ),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_Mass_Indoor"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}



fit1=best=NULL
fit1 <- stan_glmer(
  Final_Length ~ Treatment  + (1| Tank_ID  ),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_Length_Indoor"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}


fit1=best=NULL
fit1 <- stan_glmer(
  Kf ~ Treatment  + (1| Tank_ID  ),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_K_Indoor"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}



fit1=best=NULL
fit1 <- stan_glmer(
  Growth_Rate_SGRW ~ Treatment  + (1| Tank_ID  ),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Growth_Indoor"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}







# fit1=best=NULL
# fit1 <- stan_glmer(
#   Ki ~ Treatment  + (1| Tank_ID  ),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# best <- fit1
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Ini_K_Indoor"
#   ### SAVE
#   save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
# }





############ MORPHO #############
# 
# library(readr)
# Data_Morpho_2015_Indoor <- read_csv("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/data/Data_Morpho_2015_Indoor.csv")
# df <- as.data.frame(Data_Morpho_2015_Indoor)
# df$Treatment <- as.factor(df$Treatment)
# df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference
# 
# fit1=best=NULL
# fit1 <- stan_glmer(
#   WARP1 ~ Treatment * log(Mass) + (1|Tank),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# fit2 <- stan_glmer(
#   WARP1 ~ Treatment + log(Mass) + (1|Tank),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
# loo <- compare(loo1, loo2) 
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "WARP1_indoor"
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
# }
# 
# 
# 
# 
# fit1=best=NULL
# fit1 <- stan_glmer(
#   WARP2 ~ Treatment * log(Mass) + (1|Tank),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# fit2 <- stan_glmer(
#   WARP2 ~ Treatment + log(Mass) + (1|Tank),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
# loo <- compare(loo1, loo2) 
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "WARP2_indoor"
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
# }
# 
# 
# 
# 
# 
# 
# fit1=best=NULL
# fit1 <- stan_glmer(
#   WARP3 ~ Treatment * log(Mass) + (1|Tank),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# fit2 <- stan_glmer(
#   WARP3 ~ Treatment + log(Mass) + (1|Tank),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# 
# ## Compare models
# loo1 <- loo(fit1,k_threshold = 0.7) # Treatment*Time*Position +(1|Block/Section)
# loo2 <- loo(fit2,k_threshold = 0.7) # Treatment*Time+Position +(1|Block/Section)
# loo <- compare(loo1, loo2) 
# 
# # Select the best model
# # elpd_loo <- c(loo1$elpd_loo, loo2$elpd_loo, loo3$elpd_loo, loo4$elpd_loo)
# # best <- get(paste0("fit",which.max(elpd_loo)))
# looic <- c(loo1$looic, loo2$looic)
# best <- get(paste0("fit",which.min(looic))) # lower is better (as AIC)
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "WARP3_indoor"
#   ### SAVE
#   save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
# }







###______________________________________________###
###_________________OUTDOOR______________________###
###______________________________________________###




############ MASS / GROWTH #############
#Growth_Mass_Length_Ims_2015_all <- read_excel("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/data/Growth_Mass_Length_Ims_2015_all.xlsx", sheet = "River park")
Growth_Mass_Length_Ims_2015_all <- read_delim("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/data/Growth_Mass_Length_Ims_2015_all.csv", 
                                              ";", escape_double = FALSE, trim_ws = TRUE)
#Mass_Growth_River_Park_2015 <- read_excel("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/data/Mass_Growth_River_Park_2015.xlsx")
df <- as.data.frame(Growth_Mass_Length_Ims_2015_all)
#df$Initial_Mass <- as.numeric(df$Initial_Mass)
#df$Mass_at_release <- as.numeric(df$Mass_at_release)
#df$Fil_Mass <- as.numeric(df$Fil_Mass)
df$Treatment <- as.factor(df$Treatment)
#levels(df$Treatment)[levels(df$Treatment)=="SHAM"] <- "SHAM"
df <- na.omit(df)
df <- subset(df, df$Treatment !="no tag")
df <- within(df, Treatment <- relevel(Treatment, ref = 3)) # SHAM as reference
df <- droplevels(df)




#df$Ki <- K(df$Mass_at_release , df$Initial_FL)
#df$Kf <- K(df$Fil_Mass , df$Fil_FL)
fit <- lm(log(df$Initial_Mass)~log(df$Initial_FL))
df$Ki <- residuals.lm(fit)#K(df$Initial_Mass , df$Initial_FL)
fit <- lm(log(df$Fil_Mass)~log(df$Fil_FL))
df$Kf <- residuals.lm(fit)#K(df$Final_Mass , df$Final_Length)



## INITIAL MASS & LENGTH
par(mfrow=c(2,2), byrow = TRUE)

boxplot(df$Initial_Mass~df$Treatment, main="Mass at tagging",ylim=c(1,8))
t.test(df$Initial_Mass[df$Treatment=="SHAM"],df$Initial_Mass[df$Treatment=="GH"], data= df)
#fit <- bayes.t.test(df$Initial_Mass[df$Treatment=="SHAM"],df$Initial_Mass[df$Treatment=="GH"])
#plot(fit)

boxplot(df$Mass_at_release~df$Treatment, main="Mass at release",ylim=c(1,8))
t.test(df$Mass_at_release[df$Treatment=="SHAM"],df$Mass_at_release[df$Treatment=="GH"], data= df)
# #fit <- bayes.t.test(df$Mass_at_release[df$Treatment=="SHAM"],df$Mass_at_release[df$Treatment=="GH"])
# #plot(fit)

boxplot(df$Initial_FL~df$Treatment, main="Length at tagging")
t.test(df$Initial_FL[df$Treatment=="SHAM"],df$Initial_FL[df$Treatment=="GH"])
#fit <- bayes.t.test(df$Initial_FL[df$Treatment=="SHAM"],df$Initial_FL[df$Treatment=="GH"])
#plot(fit)


boxplot(df$Ki~df$Treatment, main="Condition at tagging")
t.test(df$Ki[df$Treatment=="SHAM"],df$Ki[df$Treatment=="GH"])
#fit <- bayes.t.test(df$Ki[df$Treatment=="SHAM"],df$Ki[df$Treatment=="GH"])
#plot(fit)


# fit1=best=NULL
# fit1 <- stan_glmer(
#   Mass_at_release ~ Treatment  + (1| Channel/Section),
#   #family = poisson(link = log),
#   na.action = "na.omit",
#   data = df,
#   prior = student_t(df = 7), 
#   prior_intercept = student_t(df = 7),
#   chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
#   control=list(adapt_delta=ndelta)
#   #QR = TRUE
# )
# 
# best <- fit1
# 
# ## RESULTS
# #any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
# if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
#   model.name = "Ini_Mass_Outdoor"
#   ### SAVE
#   save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
# }


## Remove all movers
df <- subset(df, Status == "Stay")

fit1=best=NULL
fit1 <- stan_glmer(
  Fil_Mass ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_Mass_Outdoor"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}


fit1=best=NULL
fit1 <- stan_glmer(
  Fil_FL ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_Length_Outdoor"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}



fit1=best=NULL
fit1 <- stan_glmer(
  Kf ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Final_K_Outdoor"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}




fit1=best=NULL
fit1 <- stan_glmer(
  Growth_Rate_SGRW_Global ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Growth_Outdoor_Global"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}


fit1=best=NULL
fit1 <- stan_glmer(
  Growth_Rate_SGRW_Tagging_Release ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Growth_Outdoor_TaggingRelease"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}


fit1=best=NULL
fit1 <- stan_glmer(
  Growth_Rate_SGRW_Release_Recapture ~ Treatment  + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

best <- fit1

## RESULTS
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "Growth_Outdoor_ReleaseRecapture"
  ### SAVE
  save(best, file=paste0("results/IMS2015_",model.name,".Rdata"))
}











############ MORPHO OUTDOOR #############

Data_Morpho_2015_Experimental_Streams <- read_csv("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/data/Data_Morpho_2015_Experimental_Streams.csv")
df <- as.data.frame(Data_Morpho_2015_Experimental_Streams)
df$Treatment <- as.factor(df$Treatment)
df <- within(df, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference


# WARPX = Treatment * log(Mass) + (1| Channel/Section    )
fit1=best=NULL
fit1 <- stan_glmer(
  WARP1 ~ Treatment * log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
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
  data = df,
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
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "WARP1_outdoor"
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
}



fit1=best=NULL
fit1 <- stan_glmer(
  WARP2 ~ Treatment * log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
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
  data = df,
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
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "WARP2_outdoor"
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
}






fit1=best=NULL
fit1 <- stan_glmer(
  WARP3 ~ Treatment * log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  prior = student_t(df = 7), 
  prior_intercept = student_t(df = 7),
  chains = CHAINS, cores = CORES,   iter = ITER, warmup = WARM, seed = SEED,
  control=list(adapt_delta=ndelta)
  #QR = TRUE
)

fit2 <- stan_glmer(
  WARP3 ~ Treatment + log(Mass) + (1| Channel/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
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
#any(summary(best)[, "Rhat"] > 1.1) # should be FALSE
if (any(summary(best)[, "Rhat"] > 1.1)==FALSE){ # should be FALSE
  model.name = "WARP3_outdoor"
  ### SAVE
  save(best,loo, file=paste0("results/IMS2015_",model.name,".Rdata"))
}














model.name <- c("Final_Mass_Indoor","Final_Length_Indoor","Final_K_Indoor","Growth_Indoor"
                ,"Final_Mass_Outdoor","Final_Length_Outdoor","Final_K_Outdoor","Growth_Outdoor_Global","Growth_Outdoor_TaggingRelease","Growth_Outdoor_ReleaseRecapture"
                ,"WARP1_indoor","WARP2_indoor","WARP3_indoor","WARP1_outdoor","WARP2_outdoor","WARP3_outdoor"
                )

## SUMMARY
sink(paste0("results/IMS2015_PHENOTYPES.txt"))
for (myvar in model.name) {
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2015_",myvar,".Rdata"))
  #best <- get(myvar)
  mcmc <- as.array(best$stanfit)
  
  ## POSTERIOR CHECK
  y <- best$y
  yrep <- posterior_predict(best)
  group <- best$data$Treatment
  
  cat("\n__________",myvar,"__________\n")

  if (!is.null(loo)) {
    
    cat("\n Model selection :\n")
    # The difference in ELPD will be negative if the expected out-of-sample predictive accuracy of the first model is higher. 
    # If the difference is be positive then the second model is preferred.
    # Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
    # That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.
    
    cat( "MODEL 1: Treatment * log(Mass) + (1|) \n")
    cat( "MODEL 2: Treatment + log(Mass) + (1|) \n")
    cat("\n Loo :\n")
    print(loo)
    cat( "The difference in ELPD (model1 - model2) will be negative if the expected out-of-sample predictive accuracy of the first model is higher (as AIC, lower is better). If the difference is be positive then the second model is preferred.\n")
    # Evaluating the expected log predictive distribution using loo reveals that the second of the two models is slightly preferred. 
    # That said, in this case the standard error of the difference in elpd is large enough (relative to the difference itself) that we can’t say there is a definitive preference for the second model.
  
  #launch_shinystan(model)
  }
  
  cat("\nBest model:\n")
  print(best$formula)
  print(best$family)
  
  
  #cat("\n____________________\n")
  cat("\n Posterior distributions for the effect of Treatments:\n")
  mcmc <- as.array(best$stanfit)
  
  TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  
  res <- round(quantile(TreatmentGH, probs=c(.025, .25, .5, .75, .975)),3) # quantiles
  p <- round(mean(TreatmentGH>0),3) # P: confidence that the parameter is positive (P>0)
  #table<- rbind(res,"P>0"=p)
  print(res)
  print(paste0("P>0 = ", p))
  
  cat("P>0: Proportion of positive posterior values, i.e. confidence that the effect is positive \n")
  
  cat( "\n Posterior check: comparison between observed values (y) and predicted values (yrep) \n")
  cat("Ratio (y/yrep) should include the value 1\n")
  # ratio y vs yrep
  r=NULL;for (i in 1:nrow(yrep)){r[i]<-mean(y/yrep[i,], na.rm=TRUE)}
  cat("Ratio (y/yrep): ", median(r),"[",quantile(r, probs=0.025),";",quantile(r, probs=0.975),"]\n")
} #End loop 
sink()





## FIGURES POSTERIOR CHECK
for (myvar in model.name) {
  #myvar= "Sum_SRP"
  pdf(paste0("results/IMS2015_",myvar,"_check.pdf"),width=8, height = 8)
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2015_",myvar,".Rdata"))
  #best <- get(myvar)
  #mcmc <- as.array(best$stanfit)
  
  ## POSTERIOR CHECK
  y <- best$y
  yrep <- posterior_predict(best)
  group <- best$data$Treatment
  
  
  #par(mfrow=c(2,2))
  #par(mar=c(3.5, 4.5, 4.5, 1.5))
  # # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
  plot1 <- ppc_stat_grouped(y, yrep, group) #ppc_stat(y, yrep)
  plot2 <- ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
  #color_scheme_set("teal")
  plot3 <- pp_check(y, yrep, ppc_dens_overlay)
  # pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
  # # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
  plot4 <- ppc_scatter_avg(y, yrep)
  # ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)
  
  grid.arrange(plot1, plot2, plot3, plot4, nrow=2, ncol=2)
  dev.off()
} 






### PLOT RESULTS

mycol <- c("#6B6B6B", "#009ACD", "#009ACD50", "#6B6B6B50")

model.name <- c(
                #"Final_Mass_Indoor","Final_Length_Indoor","Final_K_Indoor","Growth_Indoor"
                "WARP1_indoor","WARP2_indoor","WARP3_indoor"
                #"Final_Mass_Outdoor","Final_Length_Outdoor","Final_K_Outdoor"
                #,"Growth_Outdoor_Global","Growth_Outdoor_TaggingRelease","Growth_Outdoor_ReleaseRecapture"
                ,"WARP1_outdoor","WARP2_outdoor","WARP3_outdoor"
)


pdf(paste0("results/IMS2015_Morpho.pdf"),width=8, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,2), byrow = TRUE)
#par(mar=c(3.5, 4.5, 4.5, 1.5))

for (myvar in model.name) {
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2015_",myvar,".Rdata"))
  #best <- get(myvar)
  #mcmc <- as.array(best$stanfit)
  
  ## POSTERIOR CHECK
  #y <- best$y
  # yrep <- posterior_predict(best)
  # group <- best$data$Treatment
  
  # Precictive values
  nd <- data.frame(Treatment = c("SHAM","GH"))
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  #boxplot(y_rep,col=rep(1:2,5))
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  
  
gap <- c(-.3, -.1, .1, .3)

ylim <- c(0,4)
#if (myvar=="Ini_Mass_Indoor" | myvar=="Ini_Mass_Outdoor") {ylim <- c(0,6)}
if (myvar=="Final_Mass_Indoor" | myvar=="Final_Mass_Outdoor") {ylim <- c(0,25)}
if (myvar=="Final_Length_Indoor" | myvar=="Final_Length_Outdoor") {ylim <- c(50,120)}
if (myvar=="Growth_Indoor" | myvar=="Growth_Outdoor") {ylim <- c(0,4)}
plot(NULL,xlim=c(.5,3.5),ylim=ylim,ylab=myvar,xlab="",xaxt="n",main="",bty="n")
#axis(1, labels=rev(levels(best$data$Position_Bis)), at = 1:3)


points((1:2)+gap[1],y.pred["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
segments((1:2)+gap[1],y.pred["25%",1:2],(1:2)+gap[1],y.pred["75%",1:2],lwd=2,col=mycol[1:2])
segments((1:2)+gap[1],y.pred["2.5%",1:2],(1:2)+gap[1],y.pred["97.5%",1:2],lwd=1,col=mycol[1:2])

}

dev.off()











pdf(paste0("results/IMS2015_Morpho.pdf"),width=8, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,3), byrow = TRUE)
#par(mar=c(3.5, 4.5, 4.5, 1.5))


model.name <- c(
  #"WARP1_indoor","WARP2_indoor","WARP3_indoor",
  "WARP1_outdoor","WARP2_outdoor","WARP3_outdoor"
  )


for (myvar in model.name) {
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2015_",myvar,".Rdata"))
  #best <- get(myvar)
  #mcmc <- as.array(best$stanfit)
  
  ## POSTERIOR CHECK
  #y <- best$y
  # yrep <- posterior_predict(best)
  # group <- best$data$Treatment
  
  # Precictive values
  weight <- seq(2,20,1) # Final weight predicted
  nd <- data.frame(Treatment = c(rep("SHAM",length(weight)),rep("GH",length(weight))),Mass= rep(log(weight),2)) # new data
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  


plot(NULL,xlim=range(log(weight)),ylim=c(-0.05,0.05),ylab=myvar,xlab="Final weight (log)",bty="n")
points(1:length((weight)),y.pred["50%",1:length((weight))],pch=16,col=mycol[1], lwd=2,type='l')
points(1:length((weight)),y.pred["2.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
points(1:length((weight)),y.pred["50%",(1:length(weight))+10],pch=16,col=mycol[2], lwd = 2,type='l')
points(1:length((weight)),y.pred["2.5%",(1:length(weight))+10],pch=16,col=mycol[2], lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",(1:length(weight))+10],pch=16,col=mycol[2], lty=2,type='l')
#points(1:length((weight)),y.pred["97.5%",(1:length(weight))+10],pch=16,col=2, lwd = 2,type='l')
# segments(1:length((weight)),y.pred["25%",1:length((weight))],1:length((weight)),y.pred["75%",1:length((weight))],lwd=3,col=1)
# segments(1:length((weight)),y.pred["25%",(1:length(weight))+10],1:length((weight)),y.pred["75%",(1:length(weight))+10],lwd=3,col=2)
# segments(1:length((weight)),y.pred["2.5%",1:length((weight))],1:length((weight)),y.pred["97.5%",1:length((weight))],lwd=1,col=1)
# segments(1:length((weight)),y.pred["2.5%",(1:length(weight))+10],1:length((weight)),y.pred["97.5%",(1:length(weight))+10],lwd=1,col=2)
#legend("topright",legend=c("GH", "SHAM"),col=1:2,pch=c(16,16),bty="n")

}
dev.off()
