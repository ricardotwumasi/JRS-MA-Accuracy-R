###### This code is heavily guided and influenced by Chapter 10 of the Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy and its' supplementary material, 
available at https://training.cochrane.org/handbook-diagnostic-test-accuracy.
With only minor tweaks, both, the frequentist and Bayesian code, are almost exact replicas of Takwoingi et al. (2023), with our own data as input. #######

library(rjags)
library(DTAplots)
library(mcmcplots)


# Bayesian bivariate model
modelString =

"model {
#=== LIKELIHOOD ===#
for(i in 1:n) {		
TP[i] ~ dbin(se[i],pos[i])
TN[i] ~ dbin(sp[i],neg[i])

# === PRIOR FOR INDIVIDUAL LOGIT SENSITIVITY AND SPECIFICITY === #
logit(se[i]) <- l[i,1]
logit(sp[i]) <- l[i,2]
l[i,1:2] ~ dmnorm(mu[], T[,])
}

#=== HYPER PRIOR DISTRIBUTIONS POOLED LOGIT SENSITIVITY AND SPECIFICITY === #
mu[1] ~ dnorm(0,0.01)
mu[2] ~ dnorm(0,0.01) 
# Between-study variance-covariance matrix  

T[1:2,1:2]<-inverse(TAU[1:2,1:2])
TAU[1,1] <- tau[1]*tau[1]
TAU[2,2] <- tau[2]*tau[2]
TAU[1,2] <- rho*tau[1]*tau[2]	
TAU[2,1] <- rho*tau[1]*tau[2]	
#=== HYPER PRIOR DISTRIBUTIONS FOR PRECISION OF LOGIT SENSITIVITY ===# 
#=== AND LOGIT SPECIFICITY, AND CORRELATION BETWEEN THEM === #
prec[1] ~ dgamma(2,0.5)
prec[2] ~ dgamma(2,0.5)
rho ~ dunif(-1,1)

# ===  PARAMETERS OF INTEREST === #
# BETWEEN-STUDY STANDARD DEVIATION (tau) AND VARIANCE (tau.sq) OF LOGIT SENSITIVITY AND SPECIFICITY
tau[1]<-pow(prec[1],-0.5)
tau[2]<-pow(prec[2],-0.5)
tau.sq[1]<-pow(prec[1],-1)
tau.sq[2]<-pow(prec[2],-1)
# SUMMARY SENSITIVITY AND SPECIFICITY	
Summary_Se <- 1/(1+exp(-mu[1]))
Summary_Sp <- 1/(1+exp(-mu[2]))
# PREDICTED SENSITIVITY AND SPECIFICITY IN A NEW STUDY
l.predicted[1:2] ~ dmnorm(mu[],T[,])
Predicted_Se <- 1/(1+exp(-l.predicted[1]))
Predicted_Sp <- 1/(1+exp(-l.predicted[2]))
}

"

writeLines(modelString,con="model.txt")


# Data 

TP=c(21, 10, 1395, 127, 494)
FP=c(0, 3, 8, 54, 6)
FN=c(4, 6, 0, 21, 6)
TN=c(0, 46, 36, 13, 494)
pos=TP + FN
neg=TN + FP 

n=length(TP) #  Number of studies
dataList = list(TP=TP,TN=TN,n=n,pos=pos,neg=neg)

# Compile the model 
jagsModel = jags.model("model.txt",data=dataList,n.chains=3)

# Burn-in iterations 
update(jagsModel,n.iter=5000)
# Parameters to be monitored
parameters = c( "Summary_Se" , "Summary_Sp", "prec", "mu", "tau", "tau.sq", "rho",  "se","sp", "Predicted_Se", "Predicted_Sp") 

# Posterior samples
output = coda.samples(jagsModel,variable.names=parameters,n.iter=10000)

# Plots for evaluating convergence
png("Convergence.png",width = 23, height = 23, units = "cm", res=600)
par(oma=c(0,0,3,0))
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
denplot(output, parms=c("Summary_Se"), auto.layout=FALSE, main="(a)", xlab="Summary_Se", ylab="Density")
rmeanplot(output, parms=c("Summary_Se"), auto.layout=FALSE, main="(b)")
title(xlab="Iteration", ylab="Running mean")
traplot(output, parms=c("Summary_Se"), auto.layout=FALSE, main="(c)")
title(xlab="Iteration", ylab="Summary_Se")
mtext("Diagnostics for Summary_Se", side=3, line=1, outer=TRUE, cex=2)
dev.off()

gelman.diag(output)

# Summary statistics
summary(output) 

# Covariance between the posterior samples of the mean logit-transformed sensitivity and mean logit-transformed specificity.  This term is needed in order to create the SROC plot in RevMan
cov(as.matrix(output[, "mu[1]"]), as.matrix(output[, "mu[2]"]))

# Posterior density plots
denplot(output, parms=c("Summary_Se","Summary_Sp"))

# SROC plot
SROC_rjags(X=output, model="Bivariate",n=n, study_col1="blue", study_col2=rgb(0, 0, 1, 0.15), dataset=cbind(TP,FP,FN,TN), ref_std=FALSE,SROC_curve = F)





