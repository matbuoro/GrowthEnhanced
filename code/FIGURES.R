rm(list=ls())   # Clear memory

#--------------- DIRECTORY --------------#
setwd("~/Documents/RESEARCH/PROJECTS/Biodiversa/SalmoInvade/Ims/")
#setwd("~/Documents/mbuoro/Ims/")


#--------------- PACKAGES ---------------#
library(readxl)
# see https://cran.r-project.org/web/packages/rstanarm/vignettes/rstanarm.html
library("rstanarm")
#rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(loo)
library(bayesplot)



#--------------- FUNCTIONS ---------------#
invlogit<-function(x) {1/(1+exp(-(x)))} # inverse logit
logit<-function(x) {log(x/(1-x))} # logisitic transformation
gf <- function(x){if(median(x)>=0){mean(x<=0)}else{mean(x>0)}} # function to calculate the proportion of posterior with different sign as mean


mycol <- c("black", "#009ACD", "#5CACEE")




############ Figure 2 ###############
data.points = FALSE

pdf(paste("results/Figure2.pdf",sep=""),width=5, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(3,2))
par(mar=c(4.5, 4.5, 3.5, 1))



model.name <- c("Final_Mass_Indoor","Final_Mass_Outdoor"
                #,"Final_Length_Indoor","Final_Length_Outdoor"
                ,"Growth_Indoor","Growth_Outdoor_ReleaseRecapture" #"Growth_Outdoor_Global"
                #,"Final_K_Indoor" ,"Final_K_Outdoor"
                #,"Growth_Outdoor_TaggingRelease",
                #,"WARP1_indoor","WARP2_indoor","WARP3_indoor","WARP1_outdoor","WARP2_outdoor","WARP3_outdoor"
)

res=data=list()
for (myvar in model.name) {
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2015_",myvar,".Rdata"))

  # Precictive values
  nd <- data.frame(Treatment = c("SHAM","GH"))
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  #boxplot(y_rep,col=rep(1:2,5))
  y.pred=NULL
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  res[[myvar]] <- y.pred
  data[[myvar]] <- best$data
}
  
  

  ylim <- c(0,25)
  plot(NULL,xlim=c(.5,3.5),ylim=ylim,ylab="Body mass (g)",xlab="",xaxt="n",main="",bty="n")
  mtext(c("Hatchery", "Experimental"), side=1, line=1, at= c(.8,2.2), cex=0.7)
  mtext(c("conditions", "stream"), side=1, line=1.75, at= c(.8,2.2), cex=0.7)
  abline(v=1.5,lty=3)
  mtext(paste0("",letters[1],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  ## INDOOR
  gap <- c(-.3, -.1)
  gap2 <- c(-.35, -.05)
  points(1+gap[1:2],res[[1]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(1+gap[1:2],res[[1]]["25%",1:2],1+gap[1:2],res[[1]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(1+gap[1:2],res[[1]]["2.5%",1:2],1+gap[1:2],res[[1]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  if(data.points){
  for (i in 1:length(data[[1]]$Final_Mass)){
    points(1+gap2[as.numeric(data[[1]]$Treatment[i])],data[[1]]$Final_Mass[i], col=mycol[as.numeric(data[[1]]$Treatment[i])],pch=1, cex=.25)
  }}
    
  
  ## OUTDOOR
  gap <- c( .1, .3)
  gap2 <- c( .05, .35)
  points(2+gap[1:2],res[[2]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(2+gap[1:2],res[[2]]["25%",1:2],2+gap[1:2],res[[2]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(2+gap[1:2],res[[2]]["2.5%",1:2],2+gap[1:2],res[[2]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  if(data.points){
  for (i in 1:length(data[[2]]$Fil_Mass)){
    points(2+gap2[as.numeric(data[[2]]$Treatment[i])],data[[2]]$Fil_Mass[i], col=mycol[as.numeric(data[[2]]$Treatment[i])],pch=1, cex=.25)
  }}

  ylim <- c(0,4)
  plot(NULL,xlim=c(.5,3.5),ylim=ylim,ylab=expression("Growth rate (SGR," ~ "%" ~ d^{-1} ~")"), xlab="",xaxt="n",main="",bty="n")
  mtext(c("Hatchery", "Experimental"), side=1, line=1, at= c(.8,2.2), cex=0.7)
  mtext(c("conditions", "stream"), side=1, line=1.75, at= c(.8,2.2), cex=0.7)
  abline(v=1.5,lty=3)
  mtext(paste0("",letters[2],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  ## INDOOR
  gap <- c(-.3, -.1)
  gap2 <- c(-.35, -.05)
  points(1+gap[1:2],res[[3]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(1+gap[1:2],res[[3]]["25%",1:2],1+gap[1:2],res[[3]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(1+gap[1:2],res[[3]]["2.5%",1:2],1+gap[1:2],res[[3]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  if(data.points){
  for (i in 1:length(data[[3]]$Growth_Rate_SGRW)){
    points(1+gap2[as.numeric(data[[3]]$Treatment[i])],data[[3]]$Growth_Rate_SGRW[i], col=mycol[as.numeric(data[[3]]$Treatment[i])],pch=1, cex=.25)
  }}
  
  ## OUTDOOR
  gap <- c( .1, .3)
  gap2 <- c( .05, .35)
  points(2+gap[1:2],res[[4]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(2+gap[1:2],res[[4]]["25%",1:2],2+gap[1:2],res[[4]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(2+gap[1:2],res[[4]]["2.5%",1:2],2+gap[1:2],res[[4]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  if(data.points){
  for (i in 1:length(data[[4]]$Growth_Rate_SGRW_Release_Recapture)){
    points(2+gap2[as.numeric(data[[4]]$Treatment[i])],data[[4]]$Growth_Rate_SGRW_Release_Recapture[i], col=mycol[as.numeric(data[[4]]$Treatment[i])],pch=1, cex=.25)
  }}
#  dev.off()


  
### MORPHO ###
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2015_WARP2_outdoor.Rdata"))
  
  # Precictive values
  weight <- seq(4,22,1) # Final weight predicted
  nd <- data.frame(Treatment = c(rep("SHAM",length(weight)),rep("GH",length(weight))),Mass= rep(log(weight),2)) # new data
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL, re.form=NA)
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  plot(NULL,xlim=range(weight),ylim=c(-0.05,0.05),ylab="Second partial warp",xlab="Body mass (g)",bty="n")#,xaxt='n')
  mtext(paste0("",letters[3],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  # lines
  points(weight,y.pred["50%",1:length((weight))],pch=16,col=mycol[1], lwd=2,type='l')
  points(weight,y.pred["2.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
  points(weight,y.pred["97.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
  
  points(weight,y.pred["50%",(1:length(weight))+length(weight)],pch=16,col=mycol[2], lwd = 2,type='l')
  points(weight,y.pred["2.5%",(1:length(weight))+length(weight)],pch=16,col=mycol[2], lty=2,type='l')
  points(weight,y.pred["97.5%",(1:length(weight))+length(weight)],pch=16,col=mycol[2], lty=2,type='l')
  # data points
  if(data.points){
  for (i in 1:length(best$data$Mass)){
    points(best$data$Mass[i], best$data$WARP2[i], col=mycol[as.numeric(best$data$Treatment[i])],pch=1,cex=.5)
  }}
#axis(side = 1, at = 1:length(weight), labels = weight)



  ### ACTIVITY ###
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2016_Activity_lab.Rdata"))
  
  weight <- seq(3,10,1) # initial weight predicted
  nd <- data.frame(hormone = c(rep("sham",length(weight)),rep("GHA",length(weight))),initial_weight= rep(log(weight),2), scoring = "T0") # new data
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  plot(NULL,xlim=range(weight),ylim=c(1000, 8000),ylab=expression("Activity" ~ ( ~ cm.10min^{-1})),xlab="Body mass (g)",bty="n")#,xaxt="n")
  mtext(paste0("",letters[4],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  # data points
  if(data.points){
  points(best$data$initial_weight, best$data$Distance_moved, col=mycol[as.numeric(best$data$hormone)],pch=1,cex=.5)
  }
  
  points(weight,y.pred["50%",1:length((weight))],pch=16,col=mycol[1], lwd=2,type='l')
  points(weight,y.pred["2.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
  points(weight,y.pred["97.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
  
  points(weight,y.pred["50%",(1:length(weight))+length(weight)],pch=16,col=mycol[2], lwd = 2,type='l')
  points(weight,y.pred["2.5%",(1:length(weight))+length(weight)],pch=16,col=mycol[2], lty=2,type='l')
  points(weight,y.pred["97.5%",(1:length(weight))+length(weight)],pch=16,col=mycol[2], lty=2,type='l')

  #axis(side = 1, line=1, labels=weight, at=1:length(weight))
  

  #### HABITAT USE ###
  plot(NULL,xlim=c(0,1),ylim=c(0,1),ylab="",xlab="",bty="n",xaxt='n',axes=FALSE)
  mtext(paste0("",letters[5],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  

  ### EXCRETION OUTDOOR ###
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load("results/IMS2016_SRP_outdoor.Rdata")
  
  weight <- seq(3,10,1) # initial weight predicted
  #weight <- seq(.5,2,.1) # log scale
  nd <- data.frame(Treatment = c(rep("SHAM",length(weight)),rep("GH",length(weight))),Mass_t2= rep(log(weight),2)) # new data
  yrep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  y.pred <- apply(yrep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  plot(NULL,xlim=range(weight),ylim=c(-3,1),ylab=expression("P excretion" ~ (log ~ mu ~ mol.h^{-1})),xlab="Body mass (g)",bty="n")#,xaxt='n')
  # data points
  if(data.points){
    points(best$data$Mass_t2, log(best$data$Per_Capita_SRP), col=mycol[as.numeric(best$data$Treatment)],pch=1,cex=.5)
  }
  
  mtext(paste0("",letters[6],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  # lines
  points(weight,y.pred["50%",1:length((weight))],pch=16,col=mycol[1], lwd=2,type='l')
  points(weight,y.pred["2.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
  points(weight,y.pred["97.5%",1:length((weight))],pch=16,col=mycol[1], lty=2,type='l')
  points(weight,y.pred["50%",(1:length(weight))+8],pch=16,col=mycol[2], lwd = 2,type='l')
  points(weight,y.pred["2.5%",(1:length(weight))+8],pch=16,col=mycol[2], lty=2,type='l')
  points(weight,y.pred["97.5%",(1:length(weight))+8],pch=16,col=mycol[2], lty=2,type='l')
  
  #axis(side = 1, at = 1:length(weight), labels = weight)
  dev.off() # end figure 2
  
  
  
  
  
  
  
  
  
  

############## Figure 3 ################
  
  plot.pred <- function(model, ylab, xlab, main, fun){
    
    # Precictive values
    nd <- data.frame(Treatment = c(rep("SHAM",5),rep("GH",5)), Position = rep(c("A", "B", "C", "D", "E"),2))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
    y_rep <- posterior_predict(model,newdata =nd, fun = fun , re.form=NA)
    #boxplot(y_rep,col=rep(1:2,5))
    y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
    
    y.pred <- ifelse(y.pred =="-Inf",0, y.pred) 
    
    #  formula <- "Treatment * Position + (1 | Channel)"
    #  if (paste0(best$formula)[3] == formula) {
    gap <- .1
    #col <- c("#6B6B6B", "#009ACD", "#5CACEE")
    plot(NULL,xlim=c(.5,5.5),ylim=c(min(y.pred),max(y.pred)),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
    axis(1, labels=LETTERS[1:5], at = 1:5)
    points((1:5)-gap,y.pred["50%",1:5],pch=16,col=mycol[1], cex=1.5)
    points((1:5)+gap,y.pred["50%",6:10],pch=16,col=mycol[2], cex=1.5)
    segments((1:5)-gap,y.pred["25%",1:5],(1:5)-gap,y.pred["75%",1:5],lwd=3,col=mycol[1])
    segments((1:5)+gap,y.pred["25%",6:10],(1:5)+gap,y.pred["75%",6:10],lwd=3,col=mycol[2])
    segments((1:5)-gap,y.pred["2.5%",1:5],(1:5)-gap,y.pred["97.5%",1:5],lwd=1,col=mycol[1])
    segments((1:5)+gap,y.pred["2.5%",6:10],(1:5)+gap,y.pred["97.5%",6:10],lwd=1,col=mycol[2])
    #legend("topright",legend=c("SHAM", "GH"),col=col[1:2],pch=c(16,16),bty="n")
    #points(as.numeric(as.factor(model$data$Position)), log(model$data$model+.1), col=as.numeric(as.factor(model$data$Treatment)), pch=as.numeric(as.factor(model$data$Channel)))
  }
  
  data.name = "fish" # "fish" vs "nofish" # Run analysis for fish or no fish

pdf(paste("results/Figure3_",data.name,".pdf",sep=""),width=10, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,3))
par(mar=c(5.5, 4.5, 4.5, 1))

gap=c(-.2, .2)
xname="Position"

## Rhyacophilidae
load(paste0("results/IMS2015_Rhyacophilidae_",data.name,".Rdata"))
main = "" 
plot.pred(best,ylab="Rhyacophilidae (abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[1],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
points(as.numeric(as.factor(best$data$Position)) + gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Rhyacophilidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## Polycentropodidae
load(paste0("results/IMS2015_Polycentropodidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Polycentropodidae (abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[2],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Polycentropodidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## Total_Prim_Cons
load(paste0("results/IMS2015_Total_Prim_Cons_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Primary consumers (total abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[3],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Total_Prim_Cons), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## Chlo_Torch_All
load(paste0("results/IMS2015_Chlo_Torch_All_",data.name,".Rdata"))
main = "" 
plot.pred(best,ylab=expression("Primary production" ~ (log ~ mu~g ~ chlo ~ a.cm^{-2})),xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[4],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Chlo_Torch_All), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## K_CM
load(paste0("results/IMS2015_Global_decompo_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab=expression("Global decomposition" ~ (K[Inv])),xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[5],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$K_CM), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## K_FM
load(paste0("results/IMS2015_Micro_decompo_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab=expression("Microbial decomposition" ~ (K[Micro])),xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[6],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$K_FM), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

dev.off()





############## Extended data Figure 3 ################
plot.pred <- function(model, ylab, xlab, main, fun){
  
  # Precictive values
  nd <- data.frame(Treatment = c(rep("SHAM",3),rep("GH",3)), Position = rep(c("A", "B", "D"),2))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
  y_rep <- posterior_predict(model,newdata =nd, fun = fun , re.form=NA)
  #boxplot(y_rep,col=rep(1:2,5))
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  y.pred <- ifelse(y.pred =="-Inf",0, y.pred) 
  
  #  formula <- "Treatment * Position + (1 | Channel)"
  #  if (paste0(best$formula)[3] == formula) {
  gap <- .1
  #col <- c("#6B6B6B", "#009ACD", "#5CACEE")
  plot(NULL,xlim=c(.5,3.5),ylim=c(min(y.pred),max(y.pred)),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
  axis(1, labels=c("A", "B", "D"), at = 1:3)
  points((1:3)-gap,y.pred["50%",1:3],pch=16,col=mycol[1], cex=1.5)
  points((1:3)+gap,y.pred["50%",4:6],pch=16,col=mycol[2], cex=1.5)
  segments((1:3)-gap,y.pred["25%",1:3],(1:3)-gap,y.pred["75%",1:3],lwd=3,col=mycol[1])
  segments((1:3)+gap,y.pred["25%",4:6],(1:3)+gap,y.pred["75%",4:6],lwd=3,col=mycol[2])
  segments((1:3)-gap,y.pred["2.5%",1:3],(1:3)-gap,y.pred["97.5%",1:3],lwd=1,col=mycol[1])
  segments((1:3)+gap,y.pred["2.5%",4:6],(1:3)+gap,y.pred["97.5%",4:6],lwd=1,col=mycol[2])
  #legend("topright",legend=c("SHAM", "GH"),col=col[1:2],pch=c(16,16),bty="n")
  #points(as.numeric(as.factor(model$data$Position)), log(model$data$model+.1), col=as.numeric(as.factor(model$data$Treatment)), pch=as.numeric(as.factor(model$data$Channel)))
}

pdf(paste("results/Extended_data_Figure3.pdf",sep=""),width=10, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,4))
par(mar=c(5.5, 4.5, 4.5, 1))

data.name="fish"
gap=c(-.2, .2)
xname="Position"

## Rhyacophilidae
load(paste0("results/IMS2015_Surber_Rhyacophilidae_",data.name,".Rdata"))
main = "" 
plot.pred(best,ylab="Rhyacophilidae (abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
mtext(paste0("",letters[1],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Rhyacophilidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

## Polycentropodidae
load(paste0("results/IMS2015_Surber_Polycentropodidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Polycentropodidae (abundance)",xlab="",main=main, fun = NULL)
mtext(paste0("",letters[2],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Polycentropodidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

## Total_Prim_Cons
load(paste0("results/IMS2015_Surber_Total_Prim_Cons_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Primary consumers (total abundance)",xlab="",main=main, fun = NULL)
mtext(paste0("",letters[3],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Total_Prim_Cons), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)


## Chironomidae
load(paste0("results/IMS2015_Surber_Chironomidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Chironomidae (abundance)",xlab="",main=main, fun = NULL)
mtext(paste0("",letters[4],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Chironomidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

# ## Polycentropodidae
# load(paste0("results/IMS2015_Polycentropodidae_",data.name,".Rdata"))
# main = "Polycentropodidae"
# plot.pred(best,ylab="Polycentropodidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(best$data$Position)), log(best$data$Polycentropodidae), col=col[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)


## Hydropsychidae
load(paste0("results/IMS2015_Surber_Hydropsychidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Hydropsychidae (abundance)",xlab="",main=main, fun = NULL)
mtext(paste0("",letters[5],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Hydropsychidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

## Baetidae
load(paste0("results/IMS2015_Surber_Baetidae_",data.name,".Rdata"))
main =""
plot.pred(best,ylab="Baetidae (abundance)",xlab="",main=main, fun = NULL)
mtext(paste0("",letters[6],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Baetidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

## Simuliidae
load(paste0("results/IMS2015_Surber_Simuliidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Simuliidae (abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
mtext(paste0("",letters[7],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Simuliidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

## Planorbidae
load(paste0("results/IMS2015_Surber_Planorbidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Planorbidae (abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 2, cex=.7)
mtext(paste0("",letters[8],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Planorbidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

dev.off()












############## Extended data Figure 4 ################
plot.pred <- function(model, ylab, xlab, main, fun){
  
  # Precictive values
  nd <- data.frame(Treatment = c(rep("SHAM",5),rep("GH",5)), Position = rep(c("A", "B", "C", "D", "E"),2))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
  y_rep <- posterior_predict(model,newdata =nd, fun = fun , re.form=NA)
  #boxplot(y_rep,col=rep(1:2,5))
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  y.pred <- ifelse(y.pred =="-Inf",0, y.pred) 
  
  #  formula <- "Treatment * Position + (1 | Channel)"
  #  if (paste0(best$formula)[3] == formula) {
  gap <- .1
  #col <- c("#6B6B6B", "#009ACD", "#5CACEE")
  plot(NULL,xlim=c(.5,5.5),ylim=c(min(y.pred),max(y.pred)),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
  axis(1, labels=LETTERS[1:5], at = 1:5)
  points((1:5)-gap,y.pred["50%",1:5],pch=16,col=mycol[1], cex=1.5)
  points((1:5)+gap,y.pred["50%",6:10],pch=16,col=mycol[2], cex=1.5)
  segments((1:5)-gap,y.pred["25%",1:5],(1:5)-gap,y.pred["75%",1:5],lwd=3,col=mycol[1])
  segments((1:5)+gap,y.pred["25%",6:10],(1:5)+gap,y.pred["75%",6:10],lwd=3,col=mycol[2])
  segments((1:5)-gap,y.pred["2.5%",1:5],(1:5)-gap,y.pred["97.5%",1:5],lwd=1,col=mycol[1])
  segments((1:5)+gap,y.pred["2.5%",6:10],(1:5)+gap,y.pred["97.5%",6:10],lwd=1,col=mycol[2])
  #legend("topright",legend=c("SHAM", "GH"),col=col[1:2],pch=c(16,16),bty="n")
  #points(as.numeric(as.factor(model$data$Position)), log(model$data$model+.1), col=as.numeric(as.factor(model$data$Treatment)), pch=as.numeric(as.factor(model$data$Channel)))
}

pdf(paste("results/Extended_data_Figure4.pdf",sep=""),width=10, height = 8)

#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,3))
par(mar=c(5.5, 4.5, 4.5, 1))

data.name="fish"
gap=c(-.2, .2)
xname="Position"

## Rhyacophilidae
# load(paste0("results/IMS2015_Rhyacophilidae_",data.name,".Rdata"))
# main = "Rhyacophilidae" 
# plot.pred(best,ylab="Rhyacophilidae",xlab="Position",main=main, fun = NULL)
# points(as.numeric(as.factor(best$data$Position)), (best$data$Rhyacophilidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

## Chironomidae
load(paste0("results/IMS2015_Chironomidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Chironomidae (abundance)",xlab="",main=main, fun = NULL)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Chironomidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[1],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)

# ## Polycentropodidae
# load(paste0("results/IMS2015_Polycentropodidae_",data.name,".Rdata"))
# main = "Polycentropodidae"
# plot.pred(best,ylab="Polycentropodidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(best$data$Position)), log(best$data$Polycentropodidae), col=col[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)

## Hydropsychidae
load(paste0("results/IMS2015_Hydropsychidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Hydropsychidae (abundance)",xlab="",main=main, fun = NULL)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Hydropsychidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[2],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)

## Baetidae
load(paste0("results/IMS2015_Baetidae_",data.name,".Rdata"))
main =""
plot.pred(best,ylab="Baetidae (abundance)",xlab="",main=main, fun = NULL)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Baetidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[3],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)

## Simuliidae
load(paste0("results/IMS2015_Simuliidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Simuliidae (abundance)",xlab="",main=main, fun = NULL)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Simuliidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[4],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)

## Planorbidae
load(paste0("results/IMS2015_Planorbidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Planorbidae (abundance)",xlab="",main=main, fun = NULL)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Planorbidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1,cex=0.8)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[5],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)

dev.off()








############## Extended data Figure 5 ################

plot.pred <- function(model, ylab, xlab, main, fun){
  
  # Precictive values
  nd <- data.frame(Treatment = c(rep("SHAM",5),rep("GH",5)), Position = rep(c("A", "B", "C", "D", "E"),2))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
  y_rep <- posterior_predict(model,newdata =nd, fun = fun , re.form=NA)
  #boxplot(y_rep,col=rep(1:2,5))
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  
  y.pred <- ifelse(y.pred =="-Inf",0, y.pred) 
  
  #  formula <- "Treatment * Position + (1 | Channel)"
  #  if (paste0(best$formula)[3] == formula) {
  gap <- .1
  #col <- c("#6B6B6B", "#009ACD", "#5CACEE")
  plot(NULL,xlim=c(.5,5.5),ylim=c(min(y.pred),max(y.pred)),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
  axis(1, labels=LETTERS[1:5], at = 1:5)
  points((1:5)-gap,y.pred["50%",1:5],pch=16,col=mycol[1], cex=1.5)
  points((1:5)+gap,y.pred["50%",6:10],pch=16,col=mycol[2], cex=1.5)
  segments((1:5)-gap,y.pred["25%",1:5],(1:5)-gap,y.pred["75%",1:5],lwd=3,col=mycol[1])
  segments((1:5)+gap,y.pred["25%",6:10],(1:5)+gap,y.pred["75%",6:10],lwd=3,col=mycol[2])
  segments((1:5)-gap,y.pred["2.5%",1:5],(1:5)-gap,y.pred["97.5%",1:5],lwd=1,col=mycol[1])
  segments((1:5)+gap,y.pred["2.5%",6:10],(1:5)+gap,y.pred["97.5%",6:10],lwd=1,col=mycol[2])
  #legend("topright",legend=c("SHAM", "GH"),col=col[1:2],pch=c(16,16),bty="n")
  #points(as.numeric(as.factor(model$data$Position)), log(model$data$model+.1), col=as.numeric(as.factor(model$data$Treatment)), pch=as.numeric(as.factor(model$data$Channel)))
}

data.name = "nofish" # "fish" vs "nofish" # Run analysis for fish or no fish

pdf(paste("results/Extended_data_Figure5_",data.name,".pdf",sep=""),width=10, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,3))
par(mar=c(5.5, 4.5, 4.5, 1))

gap=c(-.2, .2)
xname="Position"

## Rhyacophilidae
load(paste0("results/IMS2015_Rhyacophilidae_",data.name,".Rdata"))
main = "" 
plot.pred(best,ylab="Rhyacophilidae (abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[1],""), side = 3, line = 1, outer = FALSE, at = 0, cex=1.5)
points(as.numeric(as.factor(best$data$Position)) + gap[as.numeric(as.factor(best$data$Treatment))]
       , (best$data$Rhyacophilidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## Polycentropodidae
load(paste0("results/IMS2015_Polycentropodidae_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Polycentropodidae (abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[2],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Polycentropodidae), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## Total_Prim_Cons
load(paste0("results/IMS2015_Total_Prim_Cons_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab="Primary consumers (total abundance)",xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[3],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Total_Prim_Cons), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## Chlo_Torch_All
load(paste0("results/IMS2015_Chlo_Torch_All_",data.name,".Rdata"))
main = "" 
plot.pred(best,ylab=expression("Primary production" ~ (log ~ mu~g ~ chlo ~ a.cm^{-2})),xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[4],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$Chlo_Torch_All), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## K_CM
load(paste0("results/IMS2015_Global_decompo_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab=expression("Global decomposition" ~ (K[Inv])),xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[5],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$K_CM), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

## K_FM
load(paste0("results/IMS2015_Micro_decompo_",data.name,".Rdata"))
main = ""
plot.pred(best,ylab=expression("Microbial decomposition" ~ (K[Micro])),xlab="",main=main, fun = NULL)
mtext(xname, side = 1, line = 2, outer = FALSE, at = 3, cex=.7)
mtext(paste0("",letters[6],""), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(best$data$Position))+ gap[as.numeric(as.factor(best$data$Treatment))], (best$data$K_FM), col=mycol[as.numeric(as.factor(best$data$Treatment))], pch=1)#as.numeric(as.factor(best$data$Channel)))

dev.off()









##################### Extended Data Figure 6 ###########################
#Response variables measured at T1 and T2 in the stream mesocosms. Reported values are predicted from the models

pdf(paste0("results/Extended_data_Figure6.pdf"),width=6, height = 10)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(3,2))
#par(mar=c(3.5, 4.5, 4.5, 1.5))

point= FALSE # plot the data?

model.name <- c("Nb_Chironomidae", "All_Chloro", "K_CM")
name <- c("Chironomidae (abundance)", expression("Primary production" ~ (log ~ mu~g ~ chlo ~ a.cm^{-2})), expression("Global decomposition" ~ (log ~ K[Inv])))
#model.name <- c("K_CM", "K_FM","Cyanobacteria","Green_algae","Diatoms","All_Chloro", "Nb_Chironomidae")
# expression("Global decomposition" ~ (K[Inv]))
CORES = 3

i=j=0
for (myvar in model.name) {
  i=i+1
  
  ## LOAD MODEL
  best = loo = mcmc = NULL
  load(paste0("results/IMS2016_",myvar,".Rdata"))
  #best <- get(myvar)
  mcmc <- as.array(best$stanfit)
  
  # ## POSTERIOR CHECK
  # y <- best$y
  # yrep <- posterior_predict(best)
  # group <- best$data$Treatment
  
  
  SHAM <- mcmc[,1:CORES,"(Intercept)"]
  TreatmentGH <- mcmc[,1:CORES,"TreatmentGH"]
  TreatmentGHLD <- mcmc[,1:CORES,"TreatmentGH_LD"]
  TreatmentNF <- (mcmc[,1:CORES,"TreatmentNF"])
  # Difference between treatments
  GHvsGHLD <- TreatmentGH - TreatmentGHLD;
  GHvsNF <- TreatmentGH - TreatmentNF; 
  GHLDvsNF <- TreatmentGHLD - TreatmentNF;
  

  # PLOT
  
  # Precictive values
  nd <- data.frame(Time = c(rep("T1",12), rep("T2", 12)), Treatment = c(rep("SHAM",3),rep("GH",3),rep("GH_LD",3),rep("NF",3)), Position_Bis = rep(rev(levels(best$data$Position_Bis)),8))#,Channel = rep(c(rep("A",5),rep("B",5)),2))
  y_rep <- posterior_predict(best,newdata =nd, fun = NULL , re.form=NA)
  #boxplot(y_rep,col=rep(1:2,5))
  y.pred <- apply(y_rep,2, quantile, probs=c(0.025, 0.25, 0.5,0.75, 0.975))
  y.pred.tmp <- list()
  y.pred.tmp[[1]] <- y.pred[,1:12] # predicted values at T1
  y.pred.tmp[[2]] <- y.pred[,13:24] # predicted values at T2
  
  
  for (t in 1:2){
    j=j+1
    
    gap <- c(-.3, -.1, .1, .3)
    mycol <- c("black", "#009ACD", "#009ACD50", "black")
    if (t==1) xmain= expression(T[1]) else xmain=expression(T[2])
    plot(NULL,xlim=c(.5,3.5),ylim=c(min(y.pred),max(y.pred)),ylab=name[i],xlab="",xaxt="n",main=xmain,bty="n")
    axis(1, labels=rev(levels(best$data$Position_Bis)), at = 1:3)
    mtext(paste0("",letters[j],""), side = 3, line = 1, outer = FALSE, at = .5, cex=1.5)

    # SHAM
    points((1:3)+gap[1],y.pred.tmp[[t]]["50%",1:3],pch=16,col=mycol[1], cex=1.5)
    segments((1:3)+gap[1],y.pred.tmp[[t]]["25%",1:3],(1:3)+gap[1],y.pred.tmp[[t]]["75%",1:3],lwd=3,col=mycol[1])
    segments((1:3)+gap[1],y.pred.tmp[[t]]["2.5%",1:3],(1:3)+gap[1],y.pred.tmp[[t]]["97.5%",1:3],lwd=1,col=mycol[1])
    if (point == TRUE){
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="SHAM" & Position_Bis=="Upstream")
      points(rep(1,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="SHAM" & Position_Bis=="Middle")
      points(rep(2,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="SHAM" & Position_Bis=="Downstream")
      points(rep(3,5)+gap[1], get(paste0("tmp$",myvar)), col=mycol[1], pch=3)
    }
    
    # GH
    points((1:3)+gap[2],y.pred.tmp[[t]]["50%",4:6],pch=16,col=mycol[2], cex=1.5)
    segments((1:3)+gap[2],y.pred.tmp[[t]]["25%",4:6],(1:3)+gap[2],y.pred.tmp[[t]]["75%",4:6],lwd=3,col=mycol[2])
    segments((1:3)+gap[2],y.pred.tmp[[t]]["2.5%",4:6],(1:3)+gap[2],y.pred.tmp[[t]]["97.5%",4:6],lwd=1,col=mycol[2])
    if (point == TRUE){
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH" & Position_Bis=="Upstream")
      points(rep(1,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH" & Position_Bis=="Middle")
      points(rep(2,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH" & Position_Bis=="Downstream")
      points(rep(3,5)+gap[2], get(paste0("tmp$",myvar)), col=mycol[2], pch=3)
    }
    
    # GH_LD
    points((1:3)+gap[3],y.pred.tmp[[t]]["50%",7:9],pch=16,col=mycol[3], cex=1.5)
    segments((1:3)+gap[3],y.pred.tmp[[t]]["25%",7:9],(1:3)+gap[3],y.pred.tmp[[t]]["75%",7:9],lwd=3,col=mycol[3])
    segments((1:3)+gap[3],y.pred.tmp[[t]]["2.5%",7:9],(1:3)+gap[3],y.pred.tmp[[t]]["97.5%",7:9],lwd=1,col=mycol[3])
    if (point == TRUE){
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH_LD" & Position_Bis=="Upstream")
      points(rep(1,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH_LD" & Position_Bis=="Middle")
      points(rep(2,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="GH_LD" & Position_Bis=="Downstream")
      points(rep(3,5)+gap[3], get(paste0("tmp$",myvar)), col=mycol[3], pch=3)
    }
    
    # NF
    segments((1:3)+gap[4],y.pred.tmp[[t]]["25%",10:12],(1:3)+gap[4],y.pred.tmp[[t]]["75%",10:12],lwd=3,col=mycol[4])
    segments((1:3)+gap[4],y.pred.tmp[[t]]["2.5%",10:12],(1:3)+gap[4],y.pred.tmp[[t]]["97.5%",10:12],lwd=1,col=mycol[4])
    points((1:3)+gap[4],y.pred.tmp[[t]]["50%",10:12],pch=21,bg="white", cex=1.5)
    if (point == TRUE){
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="NF" & Position_Bis=="Upstream")
      points(rep(1,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="NF" & Position_Bis=="Middle")
      points(rep(2,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
      tmp <- subset(best$data, Time==paste0("T",t) & Treatment =="NF" & Position_Bis=="Downstream")
      points(rep(3,5)+gap[4], get(paste0("tmp$",myvar)), col=mycol[4], pch=3)
    }
    
   # legend("topright",legend=c("SHAM", "GH", "GH_LD", "NF"),col=mycol,pch=rep(16,4),bty="n", cex=.75)
  } #end loop t 
  
} # end loop
  
  
dev.off()





