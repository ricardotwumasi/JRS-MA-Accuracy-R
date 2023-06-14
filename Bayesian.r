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


# Data for the anti-CCP example
TP=c(115, 110, 40, 23, 236, 74, 89, 90, 31, 69, 25, 43, 70, 167, 26, 110, 26, 64, 71, 68, 38, 42, 149, 147, 47, 24, 40, 171, 72, 7, 481, 190, 82, 865, 139, 69, 90) 
FP=c(17, 24, 5, 0, 20, 11, 4, 2, 0, 8, 2, 1, 5, 8, 8, 3, 1, 14, 2, 14, 3, 2, 7, 10, 7, 3, 11, 26, 14, 2, 23, 12, 13, 79, 7, 5, 7) 
FN=c(16, 86, 58, 7, 88, 8, 29, 50, 22, 18, 10, 63, 17, 98, 15, 148, 20, 15, 58, 35, 0, 60, 109, 35, 20, 18, 46, 60, 77, 9, 68, 105, 71, 252, 101, 107, 101) 
TN=c(73, 215, 227, 39, 231, 130, 142, 129, 75, 38, 40, 120, 228, 88, 15, 118, 56, 293, 66, 132, 73, 96, 114, 106, 375, 79, 146, 443, 298, 51, 185, 408, 301, 2218, 464, 133, 313) 
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
tiff("Figure_11_3.tiff",width = 23, height = 23, units = "cm", res=600)
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

