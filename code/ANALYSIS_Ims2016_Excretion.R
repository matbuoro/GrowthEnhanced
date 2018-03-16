rm(list=ls())   # Clear memory

#--------------- DIRECTORY --------------#
setwd("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/")

#--------------- PACKAGES ---------------#
# see https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library("rstanarm")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
set.seed(1)

# # Install from Github (development version)
# if (!require(devtools)) {
#   install.packages("devtools")
#   library(devtools)
# }
# install_github("stan-dev/rstanarm", args = "--preclean", build_vignettes = FALSE)


#--------------- FUNCTIONS ---------------#
gf <- function(x){if(mean(x)>=0){mean(x>=0)}else{mean(x<0)}} # function to calculate the proportion of posterior with same sign as mean

plot.effect <- function(model, main, formula){
  fit <- as.array(model$stanfit)
  mcmc <- as.vector(fit[,1:nchains,"TreatmentGH"])
  res <- quantile(mcmc, probs=c(.025, .25, .5, .75, .975))
  ylim=c(-1, 1)
  plot(NULL, xlim=c(0,2),ylim=ylim,xaxt="n",ylab=paste("Effect of GH (",expression(Beta),")",sep=""),main= main,bty="n") #range(mcmc)
  abline(h=0,lty=2)
  points(1,res["50%"], pch=16,cex=1.5)
  segments(1,res["2.5%"],1,res["97.5%"])
  segments(1,res["25%"],1,res["75%"],lwd=3)
  
  text(1,max(ylim), paste0(formula),cex=0.75)
  text(1,max(ylim)-.1,  (paste0(expression(beta)," = ",round(median(mcmc),2)," [",round(quantile(mcmc, probs=0.025),2),";",round(quantile(mcmc, probs=0.975),2),"]; P>0 = ",round(mean(mcmc>0),3))),cex=0.75)
  text(1,max(ylim)-.2,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
  text(1,max(ylim)-.25,"(i.e., do not overlap with 0)",cex=.4)
  
}



###______________________________________________###
###_________________DATA_________________________###
###______________________________________________###
# Load dataset:
data <- read.table("data/EXCRETION_Mat.txt",h=TRUE) # OUTDOOR
data <- as.data.frame(data)

df <- within(data, Treatment <- relevel(Treatment, ref = 3)) # SHAM as reference
df <- within(df, Treatment_bis <- relevel(Treatment_bis, ref = 2)) # SHAM as reference
df$Treatment <- df$Treatment_bis
View(df)

#-----------------ANALYSIS-------------------#
nchains  = 3 # number of chains
ncores = nchains # number of core
niter = 10000 # number of iterations
nwarm = 2000 # warmup
ndelta = 0.99



###__________________________________###
###____________ SRP    _____________ ###
###__________________________________###

SRP <- stan_glmer(
 # log(Per_Gram_SRP) ~ Treatment + (1|Block/Section),
  log(Per_Capita_SRP) ~ Treatment + log(Mass_t2) + (1|Block/Section), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  cores = ncores,
  iter = niter, warmup = nwarm, chains = nchains,
  control=list(adapt_delta=ndelta)
#  QR = TRUE
)


## RESULTS
any(summary(SRP)[, "Rhat"] > 1.1) # should be FALSE
#launch_shinystan(SRP)
summary(SRP)
fit <- as.array(SRP$stanfit)
TreatmentGH <- as.vector(fit[,1:ncores,"TreatmentGH"])
cat("Effect of GH: ", median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]")
cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentGH)*100,1),"%")


## POSTERIOR CHECK
library(bayesplot)
#toremove <- which(is.na(df$Per_Gram_SRP))
#y <- log(as.vector(na.omit(df$Per_Gram_SRP)))
toremove <- which(is.na(df$Per_Capita_SRP))
y <- log(as.vector(na.omit(df$Per_Capita_SRP)))
yrep <- posterior_predict(SRP)
group <- df$Treatment_bis[-toremove]

# # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
ppc_stat(y, yrep)
ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
#color_scheme_set("teal")
#ppc_stat_grouped(y, yrep, group)

# pp_check(y, yrep, ppc_dens_overlay)
# pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
# 
# # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
ppc_scatter_avg(y, yrep)
# ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)







###__________________________________###
###____________ NH4    _____________ ###
###__________________________________###

NH4 <- stan_glmer(
  #log(Per_Gram_NH4) ~ Treatment + (1|Block/Section),
  log(Per_Capita_NH4) ~ Treatment + log(Mass_t2) + (1|Block/Section),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  cores = ncores,
  iter = niter, warmup = nwarm, chains = nchains,
  control=list(adapt_delta=ndelta)
  #  QR = TRUE
)


## RESULTS
any(summary(NH4)[, "Rhat"] > 1.1)
#launch_shinystan(NH4)
summary(NH4)
fit <- as.array(NH4$stanfit)
TreatmentGH <- as.vector(fit[,1:3,"TreatmentGH"])
cat("Effect of GH: ", median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]")
cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentGH)*100,1),"%")


## POSTERIOR CHECK
library(bayesplot)
#toremove <- which(is.na(df$Per_Gram_NH4))
#y <- log(as.vector(na.omit(df$Per_Gram_NH4)))
toremove <- which(is.na(df$Per_Capita_NH4))
y <- log(as.vector(na.omit(df$Per_Capita_NH4)))
yrep <- posterior_predict(NH4)
group <- df$Treatment_bis[-toremove]

# # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
ppc_stat(y, yrep)
#ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
#color_scheme_set("teal")
#ppc_stat_grouped(y, yrep, group)

#pp_check(y, yrep, ppc_dens_overlay)
# pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
# 
# # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
ppc_scatter_avg(y, yrep)
# ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)








###__________________________________###
###____________ FIGURES    _____________ ###
###__________________________________###

pdf(paste0("results/Figure_Excretion_outdoor.pdf"),width=8, height=10)
par(mfrow=c(2,2))
par(mar=c(4.5, 4.5, 2.5, 1.5))

plot.effect(SRP, main="Per_Capita_SRP", paste0(SRP$formula)[3])


#weight <- seq(2,9,1) # initial weight predicted
weight <- seq(.6,2,.2) # log scale
nd <- data.frame(Treatment = c(rep("SHAM",length(weight)),rep("GH",length(weight))),Mass_t2= rep(weight,2)) # new data
yrep <- posterior_predict(SRP,newdata =nd, fun = NULL , re.form=NA)
y.pred <- apply(yrep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))

plot(NULL,xlim=c(1,length(weight)),ylim=range(y.pred),ylab="Per_Capita_SRP (log)",xlab="Weight (log)",bty="n",xaxt='n')
axis(side = 1, at = 1:length(weight), labels = weight)
# lines
points(1:length((weight)),y.pred["50%",1:length((weight))],pch=16,col=1, lwd=2,type='l')
points(1:length((weight)),y.pred["2.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["50%",(1:length(weight))+8],pch=16,col=2, lwd = 2,type='l')
points(1:length((weight)),y.pred["2.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')

# # Or points
# gap=c(-.1,.1)
# points((1:length(weight))+gap[1],y.pred["50%",(1:length(weight))],pch=16,col=1)
# points((1:length(weight))+gap[2],y.pred["50%",(1:length(weight))+8],pch=16,col=2)
# segments((1:length(weight))+gap[1],y.pred["25%",1:length((weight))],(1:length(weight))+gap[1],y.pred["75%",1:length((weight))],lwd=3,col=1)
# segments((1:length(weight))+gap[2],y.pred["25%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["75%",(1:length(weight))+8],lwd=3,col=2)
# segments((1:length(weight))+gap[1],y.pred["2.5%",1:length((weight))],(1:length(weight))+gap[1],y.pred["97.5%",1:length((weight))],lwd=1,col=1)
# segments((1:length(weight))+gap[2],y.pred["2.5%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["97.5%",(1:length(weight))+8],lwd=1,col=2)
legend("topright",legend=c("SHAM", "GH"),col=1:2,pch=c(16,16),bty="n")



plot.effect(NH4, main="Per_Capita_NH4", paste0(SRP$formula)[3])

nd <- data.frame(Treatment = c(rep("SHAM",length(weight)),rep("GH",length(weight))),Mass_t2= rep(weight,2)) # new data
yrep <- posterior_predict(NH4,newdata =nd, fun = NULL , re.form=NA)
y.pred <- apply(yrep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))

plot(NULL,xlim=c(1,length(weight)),ylim=range(y.pred),ylab="Per_Capita_NH4 (log)",xlab="Weight (log)",bty="n",xaxt='n')
axis(side = 1, at = 1:length(weight), labels = weight)
# lines
points(1:length((weight)),y.pred["50%",1:length((weight))],pch=16,col=1, lwd=2,type='l')
points(1:length((weight)),y.pred["2.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["50%",(1:length(weight))+8],pch=16,col=2, lwd = 2,type='l')
points(1:length((weight)),y.pred["2.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')

## Or points
# gap=c(-.1,.1)
# points((1:length(weight))+gap[1],y.pred["50%",(1:length(weight))],pch=16,col=1)
# points((1:length(weight))+gap[2],y.pred["50%",(1:length(weight))+8],pch=16,col=2)
# segments((1:length(weight))+gap[1],y.pred["25%",1:length((weight))],(1:length(weight))+gap[1],y.pred["75%",1:length((weight))],lwd=3,col=1)
# segments((1:length(weight))+gap[2],y.pred["25%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["75%",(1:length(weight))+8],lwd=3,col=2)
# segments((1:length(weight))+gap[1],y.pred["2.5%",1:length((weight))],(1:length(weight))+gap[1],y.pred["97.5%",1:length((weight))],lwd=1,col=1)
# segments((1:length(weight))+gap[2],y.pred["2.5%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["97.5%",(1:length(weight))+8],lwd=1,col=2)
legend("topright",legend=c("SHAM", "GH"),col=1:2,pch=c(16,16),bty="n")

dev.off()




#########################################################################################



## Load dataset
data <- read.table("data/EXCRETION_Indoor_Mat.txt",h=TRUE) # INDOOR
df <- within(data, Treatment <- relevel(Treatment, ref = 2)) # SHAM as reference
#View(df)


###__________________________________###
###____________ SRP    _____________ ###
###__________________________________###

SRP <- stan_glmer(
  # log(Per_Gram_SRP) ~ Treatment + (1|Tank),
  log(Per_Capita_SRP) ~ Treatment + log(Mass_t2) + (1|Tank), # with mass as covariate
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  cores = ncores,
  iter = niter, warmup = nwarm, chains = nchains,
  control=list(adapt_delta=ndelta)
  #  QR = TRUE
)


## RESULTS
any(summary(SRP)[, "Rhat"] > 1.1) # should be FALSE
#launch_shinystan(SRP)
summary(SRP)
fit <- as.array(SRP$stanfit)
TreatmentGH <- as.vector(fit[,1:ncores,"TreatmentGH"])
cat("Effect of GH: ", median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]")
cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentGH)*100,1),"%")


## POSTERIOR CHECK
library(bayesplot)
#toremove <- which(is.na(df$Per_Gram_SRP))
#y <- log(as.vector(na.omit(df$Per_Gram_SRP)))
toremove <- which(is.na(df$Per_Capita_SRP))
y <- log(as.vector(na.omit(df$Per_Capita_SRP)))
yrep <- posterior_predict(SRP)
group <- df$Treatment_bis[-toremove]

# # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
ppc_stat(y, yrep)
ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
#color_scheme_set("teal")
#ppc_stat_grouped(y, yrep, group)

# pp_check(y, yrep, ppc_dens_overlay)
# pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
# 
# # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
# ppc_scatter_avg(y, yrep)
# ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)







###__________________________________###
###____________ NH4    _____________ ###
###__________________________________###

NH4 <- stan_glmer(
  #log(Per_Gram_NH4) ~ Treatment + (1|Tank),
  log(Per_Capita_NH4) ~ Treatment + log(Mass_t2) + (1|Tank),
  #family = poisson(link = log),
  na.action = "na.omit",
  data = df,
  cores = ncores,
  iter = niter, warmup = nwarm, chains = nchains,
  control=list(adapt_delta=ndelta)
  #  QR = TRUE
)


## RESULTS
any(summary(NH4)[, "Rhat"] > 1.1)
#launch_shinystan(NH4)
summary(NH4)
fit <- as.array(NH4$stanfit)
TreatmentGH <- as.vector(fit[,1:3,"TreatmentGH"])
cat("Effect of GH: ", median(TreatmentGH),"[",quantile(TreatmentGH, probs=0.025),";",quantile(TreatmentGH, probs=0.975),"]")
cat("Proportion of the posterior with the same sign as the mean: ", round(gf(TreatmentGH)*100,1),"%")


## POSTERIOR CHECK
library(bayesplot)
#toremove <- which(is.na(df$Per_Gram_NH4))
#y <- log(as.vector(na.omit(df$Per_Gram_NH4)))
toremove <- which(is.na(df$Per_Capita_NH4))
y <- log(as.vector(na.omit(df$Per_Capita_NH4)))
yrep <- posterior_predict(NH4)
group <- df$Treatment_bis[-toremove]

# # The distribution of a test statistic T(yrep), or a pair of test statistics, over the simulated datasets in yrep, compared to the observed value T(y)
ppc_stat(y, yrep)
#ppc_stat_2d(y, yrep, stat = c("median", "mean")) + legend_move("bottom")
#color_scheme_set("teal")
#ppc_stat_grouped(y, yrep, group)

#pp_check(y, yrep, ppc_dens_overlay)
# pp_check(y, yrep, fun = "stat_grouped", group = group, stat = "median")
# 
# # Scatterplots of the observed data y vs. simulated/replicated data yrep from the posterior predictive distribution
ppc_scatter_avg(y, yrep)
# ppc_scatter_avg_grouped(y, yrep, group, alpha = 0.7)








###__________________________________###
###____________ FIGURES    _____________ ###
###__________________________________###

pdf(paste0("results/Figure_Excretion_Indoor.pdf"),width=8, height=10)
par(mfrow=c(2,2))
par(mar=c(4.5, 4.5, 2.5, 1.5))

plot.effect(SRP, main="Per_Capita_SRP", paste0(SRP$formula)[3])


#weight <- seq(2,9,1) # initial weight predicted
weight <- seq(.6,2,.2) # log scale
nd <- data.frame(Treatment = c(rep("SHAM",length(weight)),rep("GH",length(weight))),Mass_t2= rep(weight,2)) # new data
yrep <- posterior_predict(SRP,newdata =nd, fun = NULL , re.form=NA)
y.pred <- apply(yrep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))

plot(NULL,xlim=c(1,length(weight)),ylim=range(y.pred),ylab="Per_Capita_SRP (log)",xlab="Weight (gr)",bty="n",xaxt='n')
axis(side = 1, at = 1:length(weight), labels = weight)
# lines
points(1:length((weight)),y.pred["50%",1:length((weight))],pch=16,col=1, lwd=2,type='l')
points(1:length((weight)),y.pred["2.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["50%",(1:length(weight))+8],pch=16,col=2, lwd = 2,type='l')
points(1:length((weight)),y.pred["2.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')

# # Or points
# gap=c(-.1,.1)
# points((1:length(weight))+gap[1],y.pred["50%",(1:length(weight))],pch=16,col=1)
# points((1:length(weight))+gap[2],y.pred["50%",(1:length(weight))+8],pch=16,col=2)
# segments((1:length(weight))+gap[1],y.pred["25%",1:length((weight))],(1:length(weight))+gap[1],y.pred["75%",1:length((weight))],lwd=3,col=1)
# segments((1:length(weight))+gap[2],y.pred["25%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["75%",(1:length(weight))+8],lwd=3,col=2)
# segments((1:length(weight))+gap[1],y.pred["2.5%",1:length((weight))],(1:length(weight))+gap[1],y.pred["97.5%",1:length((weight))],lwd=1,col=1)
# segments((1:length(weight))+gap[2],y.pred["2.5%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["97.5%",(1:length(weight))+8],lwd=1,col=2)
legend("topright",legend=c("SHAM", "GH"),col=1:2,pch=c(16,16),bty="n")



plot.effect(NH4, main="Per_Capita_NH4", paste0(SRP$formula)[3])


nd <- data.frame(Treatment = c(rep("SHAM",length(weight)),rep("GH",length(weight))),Mass_t2= rep(weight,2)) # new data
yrep <- posterior_predict(NH4,newdata =nd, fun = NULL , re.form=NA)
y.pred <- apply(yrep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))

plot(NULL,xlim=c(1,length(weight)),ylim=range(y.pred),ylab="Per_Capita_NH4 (log)",xlab="Weight (gr)",bty="n",xaxt='n')
axis(side = 1, at = 1:length(weight), labels = weight)
# lines
points(1:length((weight)),y.pred["50%",1:length((weight))],pch=16,col=1, lwd=2,type='l')
points(1:length((weight)),y.pred["2.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",1:length((weight))],pch=16,col=1, lty=2,type='l')
points(1:length((weight)),y.pred["50%",(1:length(weight))+8],pch=16,col=2, lwd = 2,type='l')
points(1:length((weight)),y.pred["2.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')
points(1:length((weight)),y.pred["97.5%",(1:length(weight))+8],pch=16,col=2, lty=2,type='l')

## Or points
# gap=c(-.1,.1)
# points((1:length(weight))+gap[1],y.pred["50%",(1:length(weight))],pch=16,col=1)
# points((1:length(weight))+gap[2],y.pred["50%",(1:length(weight))+8],pch=16,col=2)
# segments((1:length(weight))+gap[1],y.pred["25%",1:length((weight))],(1:length(weight))+gap[1],y.pred["75%",1:length((weight))],lwd=3,col=1)
# segments((1:length(weight))+gap[2],y.pred["25%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["75%",(1:length(weight))+8],lwd=3,col=2)
# segments((1:length(weight))+gap[1],y.pred["2.5%",1:length((weight))],(1:length(weight))+gap[1],y.pred["97.5%",1:length((weight))],lwd=1,col=1)
# segments((1:length(weight))+gap[2],y.pred["2.5%",(1:length(weight))+8],(1:length(weight))+gap[2],y.pred["97.5%",(1:length(weight))+8],lwd=1,col=2)
legend("topright",legend=c("SHAM", "GH"),col=1:2,pch=c(16,16),bty="n")

dev.off()


