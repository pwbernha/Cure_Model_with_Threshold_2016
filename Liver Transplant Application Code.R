##############Liver Transplant Application Code#####################

####################Necessary Packages######################
library(survival)  #For fitting standard survival distributions, creating Survival objects
library(smcure)	#For fitting semiparametric AFT and Cox cure models
library(numDeriv) #Needed to find Hessian matrix
library(flexsurv) #Needed for fitting generalized gamma AFT

#############Reading in the Data (and getting in appropriate form & defining variables)#############

#Reading in liver transplant data from survival package
data("transplant")
Data2 <- transplant #(used "Data2" to differentiate dataset name from bone marrow example)

d <- 100000	#defining "infinite" cure threshold (while some individuals are known cured, not because of cure threshold

#Putting data set in correct format
Data2 <- Data2[which(Data2[,1]!='NA'),]	#eliminating missing observations
Data2 <- Data2[which(Data2[,6]!='withdraw'),]	#eliminating those withdrawn from list
Data2[,2] <- as.numeric(Data2[,2])	#making male/female a numeric variable
Data2[which(Data2[,2]==2),2] <- 0	#Redefining females as 0 rather than 2
Data2[which(Data2[,5]==0),5] <- 0.01 #defining events at time 0 as slightly above 0 for numerical reasons

#Defining dummy variables for blood types
Data2 <- cbind(Data2, as.numeric(Data2[,3]=='A'), as.numeric(Data2[,3]=='B'), as.numeric(Data2[,3]=='AB'))
colnames(Data2)[7:9] <- c("A", "B", "AB")

#Reodering data points with those experiencing event first, then those known "cured", finally those with censored obs.
Data2 <- rbind(Data2[which(Data2[,6]=='ltx'),],Data2[which(Data2[,6]=='death'),],Data2[which(Data2[,6]=='censored'),])

#defining those experiencing event (dead, even though actually alive!), those known dead (labelled "alive here!), and censored 
dead <- dim(Data2[which(Data2[,6]=='ltx'),])[[1]]
alive <- dim(Data2[which(Data2[,6]=='death'),])[[1]]
unknown <- dim(Data2[which(Data2[,6]=='censored'),])[[1]]
n <- dead+alive+unknown

#Defining covariate and response matrices/vectors
Xnew <- as.matrix(cbind(1,Data2[,c(1,2,7,8,9)]))
ynew <- as.vector(Data2[,5])


###################Finding Model Parameter Estimates and SEs for Application###############


#Getting data into a format to get reasonable initial values for beta, gamma
datas2 <- data.frame(ynew2=ynew,delta=c(rep(1,dead),rep(0,n-dead)),Xnew2=Xnew)

#Obtaining initial values for gamma, beta;  also obtain semiparametric AFT estimates treating deaths as censored values
#SEs are obtained at bottom of this section
invisible(capture.output(sm <-smcure(Surv(ynew2, delta) ~ Xnew2.age + Xnew2.sex + Xnew2.A + Xnew2.B + Xnew2.AB , cureform=~ Xnew2.age + Xnew2.sex + Xnew2.A + Xnew2.B + Xnew2.AB, data=datas2, model = "aft", link="logit", Var = FALSE)))
Initials <- as.vector(c(-sm$b,sm$beta))

#Getting estimates based on fitting a generalized gamma AFT distribution
flexsurvreg(Surv(ynew2, delta) ~ Xnew2.age + Xnew2.sex + Xnew2.A + Xnew2.B + Xnew2.AB, data=datas2, dist="gengamma")

#Updating initial estimates by getting five sets of starting values
Initial <- InitVals(Xnew,ynew,dead,alive,n,d,9,Initials)

#####Model assuming dependent censoring#####

#initializing matrix of maximum likelihood estimate candidates
estimate <- matrix(0,5,20)

#Finding maximum likelihood estimate candidates for 5 sets of starting values	
#Note: EM algorithm not used in this particular application file
for(k in 1:5){
  values <-   optim(Initial[k,],MaxReg(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=100000))
  estimate[k,] <-  c(values$par, values$value)
}

#Eliminating any maximum likelihood estimate candidates that have infinite likelihood value at parameter estimates
if(any(estimate[,20]==-Inf)) estimate <- estimate[-which(estimate[,20]==-Inf),]

#Choosing the MLEs based on the set of estimates with highest likelihood value
FinalEstimate1 <-  t(estimate[which(estimate[,20]==min(estimate[,20])),])[1:19]

#Attempting continuing maximizations to update proposed MLEs from previous set of MLEs (in this application, with so many parameters, slight improvements were made)
for(k in 1:40){
  print(k)
  FinalEstimate1 <- optim(FinalEstimate1,MaxReg(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=100000))$par
}

#Maximum likelihood value at MLE
Maxi <- Maximum(FinalEstimate1[14:15],Xnew,ynew,dead, alive, n, d, c(FinalEstimate1[1:13],FinalEstimate1[16:19]))

#Estimated SEs at MLE (note: SEs are approximate here as large parameter values of cured inds. censoring dist. cause small issues)
SE <- sqrt(diag(solve(hessian(Maximumb, c(FinalEstimate1), X2=Xnew,y2=ynew,dead=dead,alive=alive,n=n,d=d))))
sqrt(diag(solve(optim(FinalEstimate1,MaxReg(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=1),hessian=TRUE)$hessian)))


#####Reduced model fit for purpose of LRT#####

#initializing matrix of maximum likelihood estimate candidates
estimate <- matrix(0,5,18)

#Finding maximum likelihood estimate candidates for 5 sets of starting values	
#Note: EM algorithm not used in this particular application file
for(k in 1:5){
  values <- tryCatch(optim(c(FinalEstimate1[1:13], Initial[k,14:15],FinalEstimate1[18:19]),MaxRegSimp(Xnew,ynew,dead,alive,n,d,2),method="BFGS"), error=function(...) optim(c(FinalEstimate1[1:13], Initial[k,14:15],FinalEstimate1[18:19]),MaxRegSimp(Xnew,ynew,dead,alive,n,d,1), control=list(maxit=100000, reltol = sqrt(.Machine$double.eps)/20))) # nlminb(Initial[k,1:8],MaxRegInd(Xnew,ynew,dead,alive,n,d,1)))
  estimate[k,] <-  c(values$par, c(values$objective,values$value))
}

#Eliminating any maximum likelihood estimate candidates that have infinite likelihood value at parameter estimates
if(any(estimate[,18]==-Inf)) estimate <- estimate[-which(estimate[,18]==-Inf),]

#Choosing the MLEs based on the set of estimates with highest likelihood value
FinalEstimate2 <-  t(estimate[which(estimate[,18]==min(estimate[,18])),])[1:17]

#Attempting continuing maximizations to update proposed MLEs from previous set of MLEs (in this application, with so many parameters, slight improvements were made)
for(k in 1:40){
  print(k)
  FinalEstimate2 <- optim(FinalEstimate2,MaxRegSimp(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=100000))$par
  print(Maximum(FinalEstimate2[14:15],Xnew,ynew,dead, alive, n, d, c(FinalEstimate2[1:13],FinalEstimate2[16:17],FinalEstimate2[16:17])))
  
}

#Maximum likelihood value at MLE
Maxi2 <- Maximum(FinalEstimate2[14:15],Xnew,ynew,dead, alive, n, d, c(FinalEstimate2[1:13],FinalEstimate2[16:17],FinalEstimate2[16:17]))

#####Likelihood ratio test for whether censoring is different for each cure status#####
2*(Maxi2-Maxi)


#####Model assuming independent censoring#####

#initializing matrix of maximum likelihood estimate candidates
estimate <- matrix(0,5,16)

#Finding maximum likelihood estimate candidates for 5 sets of starting values	
#Note: EM algorithm not used in this particular application file
for(k in 1:5){
  values <- tryCatch(optim(c(FinalEstimate2[1:13], Initial[k,14:15]),MaxRegInd(Xnew,ynew,dead,alive,n,d,2),method="BFGS"), error=function(...) optim(c(FinalEstimate1[1:13], Initial[k,14:15]),MaxRegInd(Xnew,ynew,dead,alive,n,d,1), control=list(maxit=100000, reltol = sqrt(.Machine$double.eps)/20))) # nlminb(Initial[k,1:8],MaxRegInd(Xnew,ynew,dead,alive,n,d,1)))
  estimate[k,] <-  c(values$par, values$value)
}

#Eliminating any maximum likelihood estimate candidates that have infinite likelihood value at parameter estimates
if(any(estimate[,16]==-Inf)) estimate <- estimate[-which(estimate[,16]==-Inf),]

#Choosing the MLEs based on the set of estimates with highest likelihood value
FinalEstimate3 <-  t(estimate[which(estimate[,16]==min(estimate[,16])),])[1:15]

#Attempting continuing maximizations to update proposed MLEs from previous set of MLEs (in this application, with so many parameters, slight improvements were made)
for(k in 1:40){
  print(k)
  FinalEstimate3 <- optim(FinalEstimate3,MaxRegInd(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=100000))$par
  print(Maximum2(FinalEstimate3[14:15],Xnew,ynew,dead, alive, n, d, c(FinalEstimate3[1:13])))
  
}

#Estimated SEs at MLE 
SE2 <- sqrt(diag(solve(optim(FinalEstimate3,MaxRegInd(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=1),hessian=TRUE)$hessian)))


#Obtaining standard errors for both semiparametric PH and AFT cure models using bootstrap 

B <- 1000 # number of bootstraps
Pars <- matrix(0,B,12)
t <- 0
while(t<B+1){
  print(t)
  err <- 0 #reset "err" to 0 (was previously defined as 1 if last bootstrapped data set led to numerical problems
  t <- t+1	#counter increase by one
  
  #Choose rows to include in bootstrapped data
  h <- sort(sample(1:n,replace=TRUE))
  
  #Defining covariate and response data based on bootstrapped rows
  Data3 <- Data2[h,]
  
  #Redefining data to be in proper form & defining required variables
  dead <- dim(Data3[which(Data2[,6]=='ltx'),])[[1]]
  alive <- dim(Data3[which(Data2[,6]=='death'),])[[1]]
  unknown <- dim(Data3[which(Data2[,6]=='censored'),])[[1]]
  n <- dead+alive+unknown
  Xnew <- as.matrix(cbind(1,Data3[,c(1,2,7,8,9)]))
  ynew <- as.vector(Data3[,5])
  
  #Defining dataframe in form needed for smcure function
  datas2 <- data.frame(ynew2=ynew,delta=c(rep(1,dead),rep(0,n-dead)),Xnew2=Xnew)
  
  #Obtaining estimates for semiparametric AFT
  invisible(capture.output(sm <-smcure(Surv(ynew2, delta) ~ Xnew2.age + Xnew2.sex + Xnew2.A + Xnew2.B + Xnew2.AB , cureform=~ Xnew2.age + Xnew2.sex + Xnew2.A + Xnew2.B + Xnew2.AB, data=datas2, model = "aft", link="logit", Var = FALSE)))
  Pars[t,] <- as.vector(c(-sm$b,sm$beta))
  
}

#Obtaining bootstrap SEs
apply(Pars[1:423,],2,sd)


#######################Program Functions#########################

##################SNP Functions####################


#A matrix in SNP fitting for normal 
Anorm <- matrix(c(1,0,1,0,1,0,1,0,3),3,3) 

#Gets "c" coefficients, then "a" coefficients (See Zhang and Davidian 2001)
acoef <- function(Phi,k) {
  B <- chol(Anorm[1:(k+1),1:(k+1)])
  Binv <- solve(B)
  c <- rep(1,k+1)
  if(k>0) { for(i in 2:(k+1)) c[i] <- c[i-1]*cos(Phi[i-1])
  for(i in 1:(k)) c[i] <- c[i]*sin(Phi[i])
  }   		   
  as.numeric(Binv%*%c)
}

#Finds P_K values
Pk <- function(T,X,t,k,a,n) {
  s <- matrix(0,n,(k+1))
  for(i in 1:(k+1)) s[,i] <- ((T-X%*%c(t[1:6]))/t[7])^(i-1)
  s%*%a
}

#CDF function for standard normal-SNP given "a" values
pSNPnorm <-function(x,a) 1-(a[1]^2*pnorm(x,lower.tail=FALSE)+2*a[1]*a[2]*dnorm(x)+(a[2]^2+2*a[1]*a[3])*(x*dnorm(x)+pnorm(x,lower.tail=FALSE))+2*a[2]*a[3]*(x^2+2)*dnorm(x)+a[3]^2*(x^3*dnorm(x)+3*x*dnorm(x)+3*pnorm(x,lower.tail=FALSE)))

#Density function for normal-SNP given parameter values and data
dSNPtNorm <- function(y,X,d,k,t) {
  fc <- rep(0,length(y))
  a <- acoef(t[(dim(X)[[2]]+2):length(t)],k)
  P_k <- Pk(log(y),X,t,k,a,length(y))
  if(k==1) fc <- 1/(y*t[7])*(P_k)^2*dnorm((log(y)-as.vector(X%*%c(t[1:6])))/t[7])/pSNPnorm((log(d)-as.vector(X%*%c(t[1:6])))/t[7],c(a,0))  
  if(k==2) fc <- 1/(y*t[7])*(P_k)^2*dnorm((log(y)-as.vector(X%*%c(t[1:6])))/t[7])/pSNPnorm((log(d)-as.vector(X%*%c(t[1:6])))/t[7],a)  
  fc
}

#Potentially truncated normal-SNP CDF given data, cure threshold, and parameters
pSNPtNorm <- function(y,X,d,k,t,lower=TRUE) {
  fc <- rep(0,length(y))
  a <- acoef(t[(dim(X)[[2]]+2):length(t)],k) 
  if(k==1) fc <- (pSNPnorm(((log(y)-as.vector(X%*%c(t[1:6])))/t[7]),c(a,0)))/(pSNPnorm(((log(d)-as.vector(X%*%c(t[1:6])))/t[7]),c(a,0)))
  if(k==2) fc <- (pSNPnorm(((log(y)-as.vector(X%*%c(t[1:6])))/t[7]),a))/(pSNPnorm(((log(d)-as.vector(X%*%c(t[1:6])))/t[7]),a))
  if(lower==FALSE) 1-fc else fc
}

##################Functions for Maximization##################

#Function to derive set of initial values for maximizing SNP-AFT cure model 
#Note:  Initial values for steps 1-3 in algorithm in paper were obtain using slightly alternative process
#Note:  Input data and "divs=number of grid points in each dimension" and get 5 sets of starting values as matrix
InitVals <- function(X,y,dead, alive,n, d,divs, initials) {
  Starting <-  optim(c(initials[7:12],initials[1:6],5,1,1000,1,1000),MaxInit(X,y,dead,alive,n,d))$par #tryCatch(nlminb(c(0.2,0.1,-0.74,0.03,0.004,4,-0.03,1.17,-0.07,-0.01,1,50,0.1,50,0.1),MaxInit(X,y,dead,alive,n,d))$par, error=function(...))
  
  PhiVals <- seq(-1.50,1.50,3/(divs-1))
  PhiValsM  <- PhiI <- NULL
  for(i in 1:2) PhiValsM <- cbind(rep(rep(PhiVals,each=divs^(i-1)),divs^(2-i)),PhiValsM)
  Grid <- array(0,rep(divs,2))
  Grid <- array(apply(PhiValsM,1,Maximum, X2=X, y2=y,dead=dead,alive=alive,n=n,d=d,Starting),rep(divs,2))
  
  for(i in 1:5){ g <- as.vector(which(Grid==min(Grid), arr.in=TRUE))
  Grid[(max(1,g[1]-round(divs/5,0)):min(length(PhiVals),g[1]+round(divs/5,0))),(max(1,g[2]-round(divs/5,0)):min(length(PhiVals),g[2]+round(divs/5,0)))] <- max(Grid)
  PhiI <- rbind(PhiI,c(PhiVals[g[2]],PhiVals[g[1]]))
  }
  cbind(matrix(Starting[1:13],5,13,byrow=TRUE),PhiI,matrix(Starting[14:17],5,4,byrow=TRUE))
}

#Likelihood function for normal-SNP; outputs likelihood value at parameters given data
Maximum <- function(Phi,X2,y2,dead,alive,n,d,par){
  par <- c(par[1:13],Phi,par[14:17])	
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:6])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[7:15])*(1-pweibull(y2[1:dead],par[18],par[19])))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y2[(dead+1):(dead+alive)], par[16],par[17]))*(1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
  llike[(dead+alive+1):n] <- -log(dweibull(y2[(dead+alive+1):n],par[16],par[17])*(1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(dweibull(y2[(dead+alive+1):n],par[18],par[19]))*(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[7:15],lower=FALSE))
  return(sum(llike))
}

#Likelihood function for normal-SNP; outputs likelihood value at parameters given data (modified input from previous function)
Maximumb <- function(par,X2,y2,dead,alive,n,d){
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:6])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[7:15])*(1-pweibull(y2[1:dead],par[18],par[19])))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y2[(dead+1):(dead+alive)], par[16],par[17]))*(1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
  llike[(dead+alive+1):n] <- -log(dweibull(y2[(dead+alive+1):n],par[16],par[17])*(1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(dweibull(y2[(dead+alive+1):n],par[18],par[19]))*(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[7:15],lower=FALSE))
  return(sum(llike))
}

#Likelihood function for normal-SNP not assuming dependent censoring; outputs likelihood value at parameters given data 
Maximum2 <- function(Phi,X2,y2,dead,alive,n,d,par){
  par <- c(par[1:13],Phi,par[14:17])	
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:6])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[7:15]))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
  llike[(dead+alive+1):n] <- -log((1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[7:15],lower=FALSE))
  return(sum(llike))
}

#Likelihood function for normal-SNP not assuming dependent censoring; outputs likelihood value at parameters given data (modified input from previous function)
Maximum2b <- function(par,X2,y2,dead,alive,n,d){
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:6])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[7:15]))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
  llike[(dead+alive+1):n] <- -log((1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[7:15],lower=FALSE))
  return(sum(llike))
}

#normal-AFT likelihood function for purpose of getting initial values for parameters
MaxInit <- function(X,y,dead,alive,n,d){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:6])))^(-1)*dlnorm(y[1:dead],X[1:dead,]%*%c(par[7:12]),par[13])/plnorm(d,X[1:dead,]%*%c(par[7:12]),par[13])*(1-pweibull(y[1:dead],par[16],par[17])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y[(dead+1):(dead+alive)], par[14], par[15]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
    llike[(dead+alive+1):n] <- -log(dweibull(y[(dead+alive+1):n],par[14],par[15])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(dweibull(y[(dead+alive+1):n],par[16],par[17]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*(1-plnorm(y[(dead+alive+1):n],X[(dead+alive+1):n,]%*%c(par[7:12]),par[13])/plnorm(d,X[(dead+alive+1):n,]%*%c(par[7:12]),par[13])))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for finding MLEs
MaxReg <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:6])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[7:15])*(1-pweibull(y[1:dead],par[18],par[19])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y[(dead+1):(dead+alive)],par[16],par[17]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
    llike[(dead+alive+1):n] <- -log(dweibull(y[(dead+alive+1):n],par[16],par[17])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(dweibull(y[(dead+alive+1):n],par[18],par[19]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[7:15],lower=FALSE))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for finding MLEs in reduced model (for LRT)
MaxRegSimp <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:6])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[7:15])*(1-pweibull(y[1:dead],par[16],par[17])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y[(dead+1):(dead+alive)],par[16], par[17]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
    llike[(dead+alive+1):n] <- -log(dweibull(y[(dead+alive+1):n],par[16],par[17])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(dweibull(y[(dead+alive+1):n],par[16],par[17]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[7:15],lower=FALSE))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for find MLEs for model assuming independent censoring and time-to-event variables
MaxRegInd <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:6])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[7:15]))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:6])))^(-1))
    llike[(dead+alive+1):n] <- -log((1+exp(-X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)+(1+exp(X[(dead+alive+1):n,]%*%c(par[1:6])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[7:15],lower=FALSE))
    return(sum(llike))
  }
  like
}
