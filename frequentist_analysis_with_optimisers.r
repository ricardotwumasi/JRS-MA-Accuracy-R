###### This code is heavily guided and influenced by Chapter 10 of the Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy and its' supplementary material, 
available at https://training.cochrane.org/handbook-diagnostic-test-accuracy.
With only minor tweaks, both, the frequentist and Bayesian code, are almost exact replicas of Takwoingi et al. (2023), with our own (Lange, Koutsouleris & Twumasi, in press) data as input. #######

#initialises data, creates new dataframe 'jobreco' for study identifiers and true positives (TP), false positives (FP), true negatives (TN), and false negatives (FN) 

Study_ID  <- c("Coelho 2015", "Chala 2018", "Mpela 2020", "Guleria 2023", "Rajathilagam 2019")

TP <- c(21, 10, 1395, 127, 114)
FP <- c(0, 3, 8, 54, 2)
TN <- c(0, 46, 36, 13, 374)
FN <- c(4, 6, 0, 21, 2)
jobreco <- data.frame(Study_ID  = Study_ID , TP = TP, FP = FP, TN = TN, FN = FN)

#Data transformation - modifies the original dataframe to include total numbers of disease positive (n1) and disease negative (n0) cases, as well as renaming TP and TN for clarity in modelling (true1 and true0).

X <- jobreco
X$n1 <- X$TP+X$FN
X$n0 <- X$FP+X$TN
X$true1 <- X$TP
X$true0 <- X$TN
X$recordid <- 1:5
Y = reshape(X, direction="long", varying=list(c("n1", "n0"), c("true1", "true0")), timevar="sens", times=c(1,0), v.names=c("n","true"))
Y = Y[order(Y$id),]
Y$spec<- 1-Y$sens

#model fitting + summary and output

(MA_Y = glmer(formula=cbind(true, n - true) ~ 0 + sens + spec + (0+sens + spec|Study_ID),data=Y, family=binomial, nAGQ=1, verbose=2))
(ma_Y = summary(MA_Y))
(summary(MA_Y))$vcov
labels( ma_Y)
ma_Y$coeff
(lsens = ma_Y$coeff[1,1])
(lspec = ma_Y$coeff[2,1])
se.lsens = ma_Y$coeff[1,2]
se.lspec = ma_Y$coeff[2,2]

#Sensitivity, Specificity, Diagnostic Odds Ratio, and Likelihood Ratios

Sens = c(lsens, lsens-qnorm(0.975)*se.lsens, lsens+qnorm(0.975)*se.lsens)
Spec = c(lspec, lspec-qnorm(0.975)*se.lspec, lspec+qnorm(0.975)*se.lspec)
logit_sesp = data.frame(estimate = c(lsens, lspec),
                        lci = c(lsens-qnorm(0.975)*se.lsens, lspec-qnorm(0.975)*se.lspec), uci = c(lsens+qnorm(0.975)*se.lsens, lspec+qnorm(0.975)*se.lspec), row.names = c("lSens", "lSpec"))

plogis(Sens)
plogis(Spec)

(DOR = exp(lsens+lspec))
(LRp = plogis(lsens)/(1-plogis(lspec)))
(LRn = ((1-plogis(lsens))/plogis(lspec)))

#error estimation
library(msm)
se.logDOR = deltamethod (~ (x1+x2), mean=c(lsens,lspec), cov=ma_Y$vcov)
se.logLRp = deltamethod (~ log((exp(x1)/(1+exp(x1)))/(1- (exp(x2)/(1+exp(x2))))), mean=c(lsens,lspec), cov=ma_Y$vcov)
se.logLRn = deltamethod (~ log((1- (exp(x1)/(1+exp(x1))))/(exp(x2)/(1+exp(x2)))), mean=c(lsens,lspec), cov=ma_Y$vcov)
data.frame(estimate = c(DOR, LRp, LRn),
           lci = c(exp(log(DOR)-qnorm(0.975)*se.logDOR), exp(log(LRp)- qnorm(0.975)*se.logLRp), exp(log(LRn)-qnorm(0.975)*se.logLRn)),
           uci = c(exp(log(DOR)+qnorm(0.975)*se.logDOR), exp(log(LRp)+ qnorm(0.975)*se.logLRp), exp(log(LRn)+qnorm(0.975)*se.logLRn)),
           row.names = c("DOR", "LR+", "LR-"))


(A = glmer(formula=cbind(true, n - true) ~ 0 + sens + spec + (0+sens + spec|Study_ID), data=Y, family=binomial))

library(optimx)
library(dfoptim)

allFit(show.meth.tab=TRUE)
A.all <- allFit(A)
ss <- summary(A.all)

### logical vector: which optimisers worked?
ss$which.OK

### vector of log-likelihoods
ss$llik

### table of fixed effect
ss$fixef

### table of random effect SDs and correlations
ss$sdcor

(A.all = summary(A.all))
(A = summary(A))
