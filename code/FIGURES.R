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
library(betareg)
library(bayesplot)
library(gridExtra)


#--------------- FUNCTIONS ---------------#
invlogit<-function(x) {1/(1+exp(-(x)))} # inverse logit
logit<-function(x) {log(x/(1-x))} # logisitic transformation

gf <- function(x){if(mean(x)>=0){mean(x>=0)}else{mean(x<0)}} # function to calculate the proportion of posterior with same sign as mean

plot.effect <- function(model, main, formula){
  fit <- as.array(model$stanfit)
  mcmc <- as.vector(fit[,1:3,"TreatmentGH"])
  res <- quantile(mcmc, probs=c(.025, .25, .5, .75, .975))
  
  plot(NULL, xlim=c(0,2),ylim=c(min(res), max(res)),xaxt="n",ylab=paste("Effect of GH (",expression(beta),")",sep=""),main= main,bty="n") #range(mcmc)
  abline(h=0,lty=2)
  points(1,res["50%"], pch=16,cex=1.5)
  segments(1,res["2.5%"],1,res["97.5%"])
  segments(1,res["25%"],1,res["75%"],lwd=3)
  
  text(1,1.5, paste0(formula),cex=0.75)
  text(1,1.4,  (paste0("beta = ",round(median(mcmc),2)," [",round(quantile(mcmc, probs=0.025),2),";",round(quantile(mcmc, probs=0.975),2),"]; P>0 = ",round(mean(mcmc>0),3))),cex=0.75)
  text(1,1.3,"P: confidence that the parameter is positive (P>0) or negative (P<0)",cex=.4)
  text(1,1.25,"(i.e., do not overlap with 0)",cex=.4)
  
}

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
  plot(NULL,xlim=c(1,5),ylim=c(min(y.pred),max(y.pred)),ylab=ylab,xlab=xlab,xaxt="n",main=main,bty="n")
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




#mycol <- c("#6B6B6B", "tomato", "#5CACEE")
mycol <- c("#6B6B6B", "#009ACD", "#5CACEE")

############################################
############Phenotypes###############
############################################

model.name <- c("Final_Mass_Indoor","Final_Mass_Outdoor"
                ,"Final_Length_Indoor","Final_Length_Outdoor"
                ,"Growth_Indoor","Growth_Outdoor_Global"
                #,"Final_K_Indoor" ,"Final_K_Outdoor"
                #,"Growth_Outdoor_TaggingRelease","Growth_Outdoor_ReleaseRecapture"
                #,"WARP1_indoor","WARP2_indoor","WARP3_indoor","WARP1_outdoor","WARP2_outdoor","WARP3_outdoor"
)

res=list()
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
}
  

pdf(paste("results/IMS2015_phenotypes.pdf",sep=""),width=8, height = 5)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(1,3))
par(mar=c(5.5, 4.5, 4.5, 1))
  

  ylim <- c(0,25)
  plot(NULL,xlim=c(.5,3.5),ylim=ylim,ylab="Mass (gr)",xlab="",xaxt="n",main="",bty="n")
  mtext(c("Indoor", "Outdoor"), side=1, line=1, at= c(.8,2.2), cex=0.8)
  abline(v=1.5,lty=3)
  mtext(paste0("(",letters[1],")"), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  ## INDOOR
  gap <- c(-.3, -.1)
  points(1+gap[1:2],res[[1]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(1+gap[1:2],res[[1]]["25%",1:2],1+gap[1:2],res[[1]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(1+gap[1:2],res[[1]]["2.5%",1:2],1+gap[1:2],res[[1]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  ## OUTDOOR
  gap <- c( .1, .3)
  points(2+gap[1:2],res[[2]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(2+gap[1:2],res[[2]]["25%",1:2],2+gap[1:2],res[[2]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(2+gap[1:2],res[[2]]["2.5%",1:2],2+gap[1:2],res[[2]]["97.5%",1:2],lwd=1,col=mycol[1:2])


  ylim <- c(60,130)
  plot(NULL,xlim=c(.5,3.5),ylim=ylim,ylab="Length (mm)",xlab="",xaxt="n",main="",bty="n")
  mtext(c("Indoor", "Outdoor"), side=1, line=1, at= c(.8,2.2), cex=0.8)
  abline(v=1.5,lty=3)
  mtext(paste0("(",letters[2],")"), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  ## INDOOR
  gap <- c(-.3, -.1)
  points(1+gap[1:2],res[[3]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(1+gap[1:2],res[[3]]["25%",1:2],1+gap[1:2],res[[3]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(1+gap[1:2],res[[3]]["2.5%",1:2],1+gap[1:2],res[[3]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  ## OUTDOOR
  gap <- c( .1, .3)
  points(2+gap[1:2],res[[4]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(2+gap[1:2],res[[4]]["25%",1:2],2+gap[1:2],res[[4]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(2+gap[1:2],res[[4]]["2.5%",1:2],2+gap[1:2],res[[4]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  
  
  ylim <- c(0,4)
  plot(NULL,xlim=c(.5,3.5),ylim=ylim,ylab="Growth Rate (SGRW)",xlab="",xaxt="n",main="",bty="n")
  mtext(c("Indoor", "Outdoor"), side=1, line=1, at= c(.8,2.2), cex=0.8)
  abline(v=1.5,lty=3)
  mtext(paste0("(",letters[3],")"), side = 3, line = 1, outer = FALSE, at = 0, cex=1)
  ## INDOOR
  gap <- c(-.3, -.1)
  points(1+gap[1:2],res[[5]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(1+gap[1:2],res[[5]]["25%",1:2],1+gap[1:2],res[[5]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(1+gap[1:2],res[[5]]["2.5%",1:2],1+gap[1:2],res[[5]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  ## OUTDOOR
  gap <- c( .1, .3)
  points(2+gap[1:2],res[[6]]["50%",1:2],pch=16,col=mycol[1:2], cex=1.5)
  segments(2+gap[1:2],res[[6]]["25%",1:2],2+gap[1:2],res[[6]]["75%",1:2],lwd=2,col=mycol[1:2])
  segments(2+gap[1:2],res[[6]]["2.5%",1:2],2+gap[1:2],res[[6]]["97.5%",1:2],lwd=1,col=mycol[1:2])
  
  dev.off()

############################################
##############Community / Ecosystem ################
############################################

model.name <- c("Rhyacophilidae","Polycentropodidae","Total_Prim_Cons","Chlo_Torch_All", "K_CM", "K_FM")

pdf(paste("results/IMS2015_fig3.pdf",sep=""),width=10, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,3))
par(mar=c(5.5, 4.5, 4.5, 1))



## Rhyacophilidae
load(paste0("results/IMS2015_Rhyacophilidae_",data.name,".Rdata"))
main = "Rhyacophilidae" 
plot.pred(best,ylab="Rhyacophilidae",xlab="",main=main, fun = NULL)
mtext(paste0("(",letters[1],")"), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(df$Position)), (df$Rhyacophilidae), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## Polycentropodidae
load(paste0("results/IMS2015_Polycentropodidae_",data.name,".Rdata"))
main = "Polycentropodidae"
plot.pred(best,ylab="Polycentropodidae",xlab="",main=main, fun = NULL)
mtext(paste0("(",letters[2],")"), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(df$Position)), (df$Polycentropodidae), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## Total_Prim_Cons
load(paste0("results/IMS2015_Total_Prim_Cons_",data.name,".Rdata"))
main = "Total_Prim_Cons"
plot.pred(best,ylab="Total_Prim_Cons (log)",xlab="",main=main, fun = NULL)
mtext(paste0("(",letters[3],")"), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(df$Position)), (df$Total_Prim_Cons), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## Chlo_Torch_All
load(paste0("results/IMS2015_Chlo_Torch_All_",data.name,".Rdata"))
main = "Chlo_Torch_All" 
plot.pred(best,ylab="Chlo_Torch_All",xlab="",main=main, fun = exp)
mtext(paste0("(",letters[4],")"), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(df$Position)), (df$Chlo_Torch_All), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## K_CM
load(paste0("results/IMS2015_K_CM_",data.name,".Rdata"))
main = "K_CM"
plot.pred(best,ylab="K_CM",xlab="",main=main, fun = NULL)
mtext(paste0("(",letters[5],")"), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(df$Position)), (df$K_CM), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## K_FM
load(paste0("results/IMS2015_K_FM_",data.name,".Rdata"))
main = "K_FM"
plot.pred(best,ylab="K_FM",xlab="",main=main, fun = NULL)
mtext(paste0("(",letters[6],")"), side = 3, line = 1, outer = FALSE, at = 1, cex=1.5)
points(as.numeric(as.factor(df$Position)), (df$K_FM), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

dev.off()



########################################################################
##################### PRIMARY CONSUMERS ################################
########################################################################
pdf(paste0("results/IMS2015_PrimaryConsumers.pdf"),width=10, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,3))
par(mar=c(5.5, 4.5, 4.5, 1))



## Rhyacophilidae
load(paste0("results/IMS2015_Rhyacophilidae_",data.name,".Rdata"))
main = "Rhyacophilidae" 
plot.pred(best,ylab="Rhyacophilidae",xlab="Position",main=main, fun = NULL)
points(as.numeric(as.factor(df$Position)), (df$Rhyacophilidae), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## Chironomidae
load(paste0("results/IMS2015_Chironomidae_",data.name,".Rdata"))
main = "Chironomidae"
plot.pred(best,ylab="Chironomidae",xlab="Position",main=main, fun = NULL)
points(as.numeric(as.factor(df$Position)), (df$Chironomidae), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

# ## Polycentropodidae
# load(paste0("results/IMS2015_Polycentropodidae_",data.name,".Rdata"))
# main = "Polycentropodidae"
# plot.pred(best,ylab="Polycentropodidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Polycentropodidae), col=col[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## Hydropsychidae
load(paste0("results/IMS2015_Hydropsychidae_",data.name,".Rdata"))
main = "Hydropsychidae"
plot.pred(best,ylab="Hydropsychidae",xlab="Position",main=main, fun = NULL)
points(as.numeric(as.factor(df$Position)), (df$Hydropsychidae), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

# ## Baetidae
# load(paste0("results/IMS2015_Baetidae_",data.name,".Rdata"))
# main =  "Baetidae"
# plot.pred(best,ylab="Baetidae (log)",xlab="Position",main=main, fun = log)
# points(as.numeric(as.factor(df$Position)), log(df$Baetidae), col=col[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## Simuliidae
load(paste0("results/IMS2015_Simuliidae_",data.name,".Rdata"))
main = "Simuliidae"
plot.pred(best,ylab="Simuliidae",xlab="Position",main=main, fun = NULL)
points(as.numeric(as.factor(df$Position)), (df$Simuliidae), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

## Planorbidae
load(paste0("results/IMS2015_Planorbidae_",data.name,".Rdata"))
main = "Planorbidae"
plot.pred(best,ylab="Planorbidae",xlab="Position",main=main, fun = NULL)
points(as.numeric(as.factor(df$Position)), (df$Planorbidae), col=mycol[as.numeric(as.factor(df$Treatment))], pch=as.numeric(as.factor(df$Channel)))

dev.off()



########################################################################
##################### Extended Data Figure 7 ###########################
########################################################################
#Extended Data Figure 7: Response variables measured at T1 and T2 in the stream mesocosms. Reported values are predicted from the models


pdf(paste0("results/IMS2016_ExtendedData_Figure7.pdf"),width=8, height = 8)
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mfrow=c(2,2))
#par(mar=c(3.5, 4.5, 4.5, 1.5))

point= FALSE # plot the data?

model.name <- c("All_Chloro", "Nb_Chironomidae")
#model.name <- c("K_CM", "K_FM","Cyanobacteria","Green_algae","Diatoms","All_Chloro", "Nb_Chironomidae")

CORES = 3

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
    gap <- c(-.3, -.1, .1, .3)
    mycol <- c("#6B6B6B", "#009ACD", "#009ACD50", "#6B6B6B50")
    plot(NULL,xlim=c(.5,3.5),ylim=c(min(y.pred),max(y.pred)),ylab=myvar,xlab="",xaxt="n",main=paste0("T",t),bty="n")
    axis(1, labels=rev(levels(best$data$Position_Bis)), at = 1:3)
    
    
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
    points((1:3)+gap[4],y.pred.tmp[[t]]["50%",10:12],pch=16,col=mycol[4], cex=1.5)
    segments((1:3)+gap[4],y.pred.tmp[[t]]["25%",10:12],(1:3)+gap[4],y.pred.tmp[[t]]["75%",10:12],lwd=3,col=mycol[4])
    segments((1:3)+gap[4],y.pred.tmp[[t]]["2.5%",10:12],(1:3)+gap[4],y.pred.tmp[[t]]["97.5%",10:12],lwd=1,col=mycol[4])
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




