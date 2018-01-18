####
#Top level script for manuscript:
#Title: Random measurement error: why worry? An example of cardiovascular risk factors
#Authors: Timo B Brakenhoff, Maarten van Smeden, Frank LJ Visseren, Rolf HH Groenwold
#Corresponding author: TB Brakenhoff; t.brakenhoff@gmail.com

#THis script is to load the clinical data, run the simulation,
#Plot the results from the simulation.
####
## Require PACKAGES
####

require(survival)
require(plyr)
require(dplyr)
require(devtools)
require(xtable)
require(DT)
require(ggplot2)
require(MASS)
require(tableone)

#####
## SET TOP WD
####

setwd("F:\\MIR5_SMART\\GitHubDataPackage\\")

#####
## Source Sim functions
#####

source("Analysis-R/R_Scripts/sim_smart_25.R") #Simulation function
source("Analysis-R/R_Scripts/plot_smart_28.R") #plot function

#####
## General Functions
#####

spec_dec <- function(x, k) format(round(x, k), nsmall=k)

####
## open data
####

com.c.m.sub <- readRDS("Data\\smart_com_c_m_sub.rds")


####
## BASELINE TABLE
####

#Select relevant vars and order
com.tab.dat <- com.c.m.sub[,c("age","sex","bp_s","bp_d","imt","abi",
                              "fu_vascom","vascom")]

#rename variables
names(com.tab.dat) <- c("Age in years","Sex","Systolic BP in mmHg",
                        "Diastolic BP in mmHg","Average IMT in mm",
                        "Average ABI (at rest)","Follow up in days",
                        "Vascular complication")

#make table1
tab1 <- CreateTableOne(data = com.tab.dat)

#Export to latex
xtable(print(tab1, nonnormal = "Follow up in days"))

#Fix further in latex and take inverse of sex variable. Sex = 2 = female

####
## HR TABLE (Now also for confounder relationships)
####

#HRs are taken from the coxph object and not cph() such that the HR is for one unit increase 
#and not for an IQR increase which is done by the summary(cph.object)

## MODEL1 : vascom _ bp s+bp_d+age+sex (DONE)

#crude models
bps.cru.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                             event = as.numeric(com.c.m.sub$vascom))~bp_s,
                        data=com.c.m.sub)

bpd.c.cru.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                             event = as.numeric(com.c.m.sub$vascom))~bp_d,
                        data=com.c.m.sub)

#adjusted models
bps.bpd.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                             event = as.numeric(com.c.m.sub$vascom))~
                          bp_s+bp_d+age+sex,data=com.c.m.sub)


#crude effects
bps.cru.ef <- summary(bps.cru.mod)$conf.int[1]
bps.cru.se <- exp(summary(bps.cru.mod)$coefficients[3])
bps.cru.ci <- summary(bps.cru.mod)$conf.int[3:4]

bpd.c.cru.ef <- summary(bpd.c.cru.mod)$conf.int[1]
bpd.c.cru.se <- exp(summary(bpd.c.cru.mod)$coefficients[3])
bpd.c.cru.ci <- summary(bpd.c.cru.mod)$conf.int[3:4]

#adjusted effects
bps.bpd.ef <- summary(bps.bpd.mod)$conf.int[1,1]
bps.bpd.se <- exp(summary(bps.bpd.mod)$coefficients[1,3])
bps.bpd.ci <- summary(bps.bpd.mod)$conf.int[1,3:4]

bpd.bpd.ef <- summary(bps.bpd.mod)$conf.int[2,1]
bpd.bpd.se <- exp(summary(bps.bpd.mod)$coefficients[2,3])
bpd.bpd.ci <- summary(bps.bpd.mod)$conf.int[2,3:4]


#crude On 10 mmHg scale
bps.cru.ef.10 <- spec_dec(bps.cru.ef^10,2)
bps.cru.se.10 <- spec_dec(bps.cru.se^10,2)
bps.cru.ci.10 <- spec_dec(bps.cru.ci^10,2)

bpd.c.cru.ef.10 <- spec_dec(bpd.c.cru.ef^10,2)
bpd.c.cru.se.10 <- spec_dec(bpd.c.cru.se^10,2)
bpd.c.cru.ci.10 <- spec_dec(bpd.c.cru.ci^10,2)

#adjusted On 10 mmHg scale
bps.bpd.ef.10 <- spec_dec(bps.bpd.ef^10,2)
bps.bpd.se.10 <- spec_dec(bps.bpd.se^10,2)
bps.bpd.ci.10 <- spec_dec(bps.bpd.ci^10,2)

bpd.bpd.ef.10 <- spec_dec(bpd.bpd.ef^10,2)
bpd.bpd.se.10 <- spec_dec(bpd.bpd.se^10,2)
bpd.bpd.ci.10 <- spec_dec(bpd.bpd.ci^10,2)


## MODEL2 : vascom _ bp s+ abi+ age+ sex 

#Crude models
abi.c.cru.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                               event = as.numeric(com.c.m.sub$vascom))~abi,
                          data=com.c.m.sub)

#adjusted models
bps.abi.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                             event = as.numeric(com.c.m.sub$vascom))~
                          bp_s+abi+age+sex,data=com.c.m.sub)


#crude effects
abi.c.cru.ef <- summary(abi.c.cru.mod)$conf.int[1]
abi.c.cru.se <- exp(summary(abi.c.cru.mod)$coefficients[3])
abi.c.cru.ci <- summary(abi.c.cru.mod)$conf.int[3:4]


#adjusted effects
bps.abi.ef <- summary(bps.abi.mod)$conf.int[1,1]
bps.abi.se <- exp(summary(bps.abi.mod)$coefficients[1,3])
bps.abi.ci <- summary(bps.abi.mod)$conf.int[1,3:4]

abi.abi.ef <- summary(bps.abi.mod)$conf.int[2,1]
abi.abi.se <- exp(summary(bps.abi.mod)$coefficients[2,3])
abi.abi.ci <- summary(bps.abi.mod)$conf.int[2,3:4]

#Crude effects (no transformation necessary)
spec_dec(abi.c.cru.ef,2)
spec_dec(abi.c.cru.ci,2)

#adjusted effects (transformation for bps on 10 mmHg scale)
bps.abi.ef.10 <- spec_dec(bps.abi.ef^10,2)
bps.abi.se.10 <- spec_dec(bps.abi.se^10,2)
bps.abi.ci.10 <- spec_dec(bps.abi.ci^10,2)

spec_dec(abi.abi.ef,2)
spec_dec(abi.abi.ci,2)


## MODEL3: vascom _ imt + bp s +age+sex (DONE)

#crude models
imt.cru.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                             event = as.numeric(com.c.m.sub$vascom))~imt,
                        data=com.c.m.sub)

bps.c.cru.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                             event = as.numeric(com.c.m.sub$vascom))~bp_s,
                        data=com.c.m.sub)

#adjusted models
imt.bps.mod    <- coxph(Surv(time  = com.c.m.sub$fu_vascom, 
                             event = as.numeric(com.c.m.sub$vascom))~
                          imt+bp_s+age+sex,data=com.c.m.sub)

#crude effects
imt.cru.ef <- summary(imt.cru.mod)$conf.int[1]
imt.cru.se <- exp(summary(imt.cru.mod)$coefficients[3])
imt.cru.ci <- summary(imt.cru.mod)$conf.int[3:4]

bps.c.cru.ef <- summary(bps.c.cru.mod)$conf.int[1]
bps.c.cru.se <- exp(summary(bps.c.cru.mod)$coefficients[3])
bps.c.cru.ci <- summary(bps.c.cru.mod)$conf.int[3:4]

#adjusted effects
imt.bps.ef <- summary(imt.bps.mod)$conf.int[1,1]
imt.bps.se <- exp(summary(imt.bps.mod)$coefficients[1,3])
imt.bps.ci <- summary(imt.bps.mod)$conf.int[1,3:4]

bps.bps.ef <- summary(imt.bps.mod)$conf.int[2,1]
bps.bps.se <- exp(summary(imt.bps.mod)$coefficients[2,3])
bps.bps.ci <- summary(imt.bps.mod)$conf.int[2,3:4]

#crude on 10mmHg scale
bps.c.cru.ef.10 <- spec_dec(bps.c.cru.ef^10,2)
bps.c.cru.se.10 <- spec_dec(bps.c.cru.se^10,2)
bps.c.cru.ci.10 <- spec_dec(bps.c.cru.ci^10,2)

#adjusted on 10mmHg scale
bps.bps.ef.10 <- spec_dec(bps.bps.ef^10,2)
bps.bps.se.10 <- spec_dec(bps.bps.se^10,2)
bps.bps.ci.10 <- spec_dec(bps.bps.ci^10,2)

####
## UNIVARIATE RELATIONSHIP CONFOUNDERS WITH EXPOSURE
####

### Model 1 corr (pearson) and LR of SBP ~ DBP
#cor
cor(com.c.m.sub$bp_s,com.c.m.sub$bp_d)

#lr
mod1.lr    <- lm(bp_s~bp_d,data=com.c.m.sub)
mod1.lr.ef <- summary(mod1.lr)$coefficients[2,1]
mod1.lr.se <- summary(mod1.lr)$coefficients[2,2]
mod1.lr.ci <- confint(mod1.lr,parm="bp_d")

mod1.lr.ef.10 <- 10*mod1.lr.ef

mod1.lr.ci.10 <- 10*mod1.lr.ci


### Model 2 corr (pearson) and LR of SBP ~ ABI
#cor
cor(com.c.m.sub$bp_s,com.c.m.sub$abi)

#lr
mod2.lr    <- lm(bp_s~abi,data=com.c.m.sub)
mod2.lr.ef <- summary(mod2.lr)$coefficients[2,1]
mod2.lr.se <- summary(mod2.lr)$coefficients[2,2]
mod2.lr.ci <- confint(mod2.lr,parm="abi")



### Model 3 corr (pearson) and LR of CIMT ~ SBP
#cor
cor(com.c.m.sub$imt,com.c.m.sub$bp_s)

#lr
mod3.lr    <- lm(imt~bp_s,data=com.c.m.sub)
mod3.lr.ef <- summary(mod3.lr)$coefficients[2,1]
mod3.lr.se <- summary(mod3.lr)$coefficients[2,2]
mod3.lr.ci <- confint(mod3.lr,parm="bp_s")

mod3.lr.ef.10 <- 10*mod3.lr.ef

mod3.lr.ci.10 <- 10*mod3.lr.ci


####
## COXPH ASSUMPTIONS OF ALL 3 MODELS
####

### ASSUMPTIONS

### MODEL 1

## martingale resids
fit1 <- coxph(Surv(com.c.m.sub$fu_vascom,as.numeric(com.c.m.sub$vascom))~bp_s,
      data=com.c.m.sub)
residM.fit1<-resid(fit1,type="martingale")

par(mfrow=c(1,1))
plot(com.c.m.sub$bp_d,residM.fit1)
lines(lowess(com.c.m.sub$bp_d,residM.fit1))
plot(com.c.m.sub$age,residM.fit1)
lines(lowess(com.c.m.sub$age,residM.fit1))
#LOESS line is straight indicating good functional form of bp_s and age (untransformed)

## PH assumption check

#crude model with only exposure
print(coxzph.fit1.I <- cox.zph(fit1, transform="identity"))
print(coxzph.fit1.KM <- cox.zph(fit1, transform="km"))
print(coxzph.fit1.L <- cox.zph(fit1, transform=log))
par(mfrow=c(1,3))	# graphs are plotted in 1 row, 3 next to each other
plot(coxzph.fit1.I);plot(coxzph.fit1.KM);plot(coxzph.fit1.L)
#In general, PH visually holds

#full adjusted model with all confs
print(coxzph.fit1.c.I <- cox.zph(bps.bpd.mod, transform="identity"))
print(coxzph.fit1.c.KM <- cox.zph(bps.bpd.mod, transform="km"))
print(coxzph.fit1.c.L <- cox.zph(bps.bpd.mod, transform=log))
par(mfrow=c(1,1))	# graphs are plotted in 1 row, 3 next to each other
plot(coxzph.fit1.c.I);plot(coxzph.fit1.c.KM);plot(coxzph.fit1.c.L)
#In general, PH visually holds

### MODEL 2

## Martingale resids
plot(com.c.m.sub$abi,residM.fit1)
lines(lowess(com.c.m.sub$abi,residM.fit1))
#LOESS line is straight indicating good functional form of abi (untransformed)

## PH assumption check

#full adjusted model with all confs
print(coxzph.fit2.c.I <- cox.zph(bps.abi.mod, transform="identity"))
print(coxzph.fit2.c.KM <- cox.zph(bps.abi.mod, transform="km"))
print(coxzph.fit2.c.L <- cox.zph(bps.abi.mod, transform=log))
par(mfrow=c(1,3))	# graphs are plotted in 1 row, 3 next to each other
plot(coxzph.fit2.c.I);plot(coxzph.fit2.c.KM);plot(coxzph.fit2.c.L)
#In general, PH visually holds

### MODEL 3

#martingale resids
fit3 <- coxph(Surv(com.c.m.sub$fu_vascom,as.numeric(com.c.m.sub$vascom))~imt,
              data=com.c.m.sub)
residM.fit3<-resid(fit3,type="martingale")

par(mfrow=c(1,2))
plot(com.c.m.sub$bp_s,residM.fit3)
lines(lowess(com.c.m.sub$bp_s,residM.fit3))
plot(com.c.m.sub$age,residM.fit3)
lines(lowess(com.c.m.sub$age,residM.fit3))
#LOESS line is straight indicating good functional form of bp_s and age (untransformed)

## PH Assumption check

#crude model with only exposure
print(coxzph.fit3.I <- cox.zph(fit3, transform="identity"))
print(coxzph.fit3.KM <- cox.zph(fit3, transform="km"))
print(coxzph.fit3.L <- cox.zph(fit3, transform=log))
par(mfrow=c(1,3))	# graphs are plotted in 1 row, 3 next to each other
plot(coxzph.fit3.I);plot(coxzph.fit3.KM);plot(coxzph.fit3.L)
#Straight line for all


#full adjusted model with all confs
print(coxzph.fit3.c.I <- cox.zph(imt.bps.mod, transform="identity"))
print(coxzph.fit3.c.KM <- cox.zph(imt.bps.mod, transform="km"))
print(coxzph.fit3.c.L <- cox.zph(imt.bps.mod, transform=log))
par(mfrow=c(1,3))	# graphs are plotted in 1 row, 3 next to each other
plot(coxzph.fit3.c.I);plot(coxzph.fit3.c.KM);plot(coxzph.fit3.c.L)
#In general, PH visually holds

####
## STANDARDIZE VARIABLES (Fixed)
####

#Standardize continuous variables 
com.c.m.st     <- cbind(scale(com.c.m.sub[,c(4:6,8:9)]),com.c.m.sub[,-c(4:6,8:9)])

#####################
## EXECUTE ME FUN ###
#####################

#Relationships to run
#MOD1: vascom ~ bp_s + bp_d + age + sex
#MOD2: vascom ~ bp_s + abi + age + sex 
#MOD3: vascom ~ imt + bp_s + age + sex 

#DO plotting manually after all simulation to make consistent heatmap limits
# So plt = sav = F

#Input parameters
df       <- com.c.m.st
out      <- "vascom"
ME_exp_sd  <- seq(0,0.50,0.10)
ME_conf_sd <- seq(0,0.50,0.10)
sims = 1000
seed = 333

#Start timer
start.t <- Sys.time()

##### BP_S #######
expo     <- "bp_s"

#MOD1: vascom ~ bp_d + plsprs + age + sex
confs    <- c("bp_d","age","sex")
ME_conf  <- "bp_d"

mod.1 <- me_mod(df=df,out=out,expo=expo,confs=confs,ME_conf=ME_conf,
                ME_exp_sd=ME_exp_sd,ME_conf_sd=ME_conf_sd,sims=sims,
                plt=F,sav=F,seed=seed)

mod.1$res.mat[1,1,1] #Reference logHR effect (left bottom cell)

#MOD2: vascom ~ bp_s + abi + age + sex 
confs    <- c("abi","age","sex")
ME_conf  <- "abi"

mod.2 <- me_mod(df=df,out=out,expo=expo,confs=confs,ME_conf=ME_conf,
                ME_exp_sd=ME_exp_sd,ME_conf_sd=ME_conf_sd,sims=sims,
                plt=F,sav=F,seed=seed)

mod.2$res.mat[1,1,1]   #Reference logHR effect (left bottom cell)
#mod.1$gheat.mod

##### IMT #######
expo     <- "imt"

#3: vascom ~ imt + bp_s + age + sex
confs    <- c("bp_s","age","sex")
ME_conf  <- "bp_s"

mod.3 <- me_mod(df=df,out=out,expo=expo,confs=confs,ME_conf=ME_conf,
                ME_exp_sd=ME_exp_sd,ME_conf_sd=ME_conf_sd,sims=sims,
                plt=F,sav=F,seed=seed)

mod.3$res.mat[1,1,1] #Reference logHR effect (left bottom cell)
#mod.3$gheat.mod

#End timer, compute total time taken
end.t <- Sys.time()
print(tot.t <- end.t-start.t)

####### SAVE RESULTS #########
sav.im <- paste0("Analysis-R\\R_Data\\","SMART_sim_",Sys.Date(),".rds")
saveRDS(list(mod.1,mod.2,mod.3,tot.t),file=sav.im)


##################################
### LOAD SIM DATA  if necessary ##
##################################

foo = readRDS("Analysis-R\\R_Data\\SMART_sim_2017-07-21.rds")

mod.1 <- foo[[1]]
mod.2 <- foo[[2]]
mod.3 <- foo[[3]]
tot.time   <- foo[[4]]

######################
### MANUAL PLOTTING ##
######################

#First find all maximum percentages for all simulations, combine in vector
lim.vec <- c(round_any(range(mod.1$tot.mat[,"exp_dif_rel"])[1], 1, f = floor),
             round_any(range(mod.1$tot.mat[,"exp_dif_rel"])[2], 1, f = ceiling),
             round_any(range(mod.2$tot.mat[,"exp_dif_rel"])[1], 1, f = floor),
             round_any(range(mod.2$tot.mat[,"exp_dif_rel"])[2], 1, f = ceiling),
             round_any(range(mod.3$tot.mat[,"exp_dif_rel"])[1], 1, f = floor),
             round_any(range(mod.3$tot.mat[,"exp_dif_rel"])[2], 1, f = ceiling))

#Select maximum limits for consistent heatmap colours across all plots
lims <- c(-max(abs(lim.vec)),max(abs(lim.vec)))
mdpt   <- 0

#Make plots (manually change mod.plt for each model run)
mod.plt <- mod.2

expo=mod.plt$in.list$expo
confs=mod.plt$in.list$confs
ME_conf=mod.plt$in.list$ME_conf
sims=mod.plt$in.list$sims
sav=F
tot.mat=mod.plt$tot.mat

plot_mod(expo=expo,confs=confs,ME_conf=ME_conf,tot.mat=tot.mat,sims=sims,sav=sav,lims=lims,mdpt=mdpt)


### LogHR Values for different cells
mod.plt$tot.mat