##############Bone Marrow Transplant Application Code#####################


####################Necessary Packages######################
library(survival)  #For fitting standard survival distributions, creating Survival objects
library(smcure)	#For fitting semiparametric AFT and Cox cure models
library(numDeriv) #Needed to find Hessian matrix

#############Reading in the Data (and getting in appropriate form & defining variables)#############
#Reading in the Bone Marrow Data
Data <- read.table("Avalos.txt")


d <- 100000  #defining an "infinite" cure threshold since none exists in this application
n <- dim(Data)[[1]]	#defining the sample size
X <- as.matrix(cbind(1,Data[,c(1,2,5,6)]))  #Defining the covariate matrix
y <-  as.vector(Data[,3])	#Defining the time-to-event/censoring


#Reordering the data based on those individuals who experience the event, are known to survival long-term (none known here), censored
Xnew <- rbind(X[which(Data[,4]==1),], X[which(y==d),],X[which(Data[,4]==0),])
ynew <- c(y[which(Data[,4]==1)],y[which(y==d)],y[which(Data[,4]==0)])

#Defining number of individuals known to experience event, to be alive, to have censored event time
dead <- length(y[which(Data[,4]==1)])
alive <- length(y[which(y==d)])
unknown <- n-dead-alive


###################Finding Model Parameter Estimates and SEs for Application###############

#Finding 5x18 matrix of startinv values (each row corresponds to set of starting values for parameters
Initial <- InitVals(Xnew,ynew,dead,alive,n,d,9)


#####Model assuming dependent censoring#####

#initializing matrix of maximum likelihood estimate candidates
estimate <- matrix(0,5,18)

#Finding maximum likelihood estimate candidates for 5 sets of starting values	
#Note: EM algorithm not used in this particular application file
for(k in 1:5){
  values <- optim(Initial[k,],MaxReg(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=10000))
  estimate[k,] <-  c(values$par, values$value)
}

#Eliminating any maximum likelihood estimate candidates that have infinite likelihood value at parameter estimates
if(any(estimate[,18]==-Inf)) estimate <- estimate[-which(estimate[,18]==-Inf),]

#Choosing the MLEs based on the set of estimates with highest likelihood value
FinalEstimate1 <-  t(estimate[which(estimate[,18]==min(estimate[,18])),])[1:17]

#Attempting continuing maximizations to update proposed MLEs from previous set of MLEs (in this application, with so many parameters, slight improvements were made)
for(k in 1:40){
  print(k)
  FinalEstimate1 <- optim(c(FinalEstimate1),MaxReg(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=10000))$par
}

#Maximum likelihood value at MLE
Maxi <- Maximum(FinalEstimate1[12:13],Xnew,ynew,dead, alive, n, d, c(FinalEstimate1[1:11],FinalEstimate1[14:17]))

#Estimated SEs at MLE (note: SEs are approximate here as large parameter values of cured inds. censoring dist. cause small issues)
SE <- sqrt(diag(solve(hessian(Maximumb, c(FinalEstimate1), X2=Xnew,y2=ynew,dead=dead,alive=alive,n=n,d=d)[1:15,1:15])))


#####Reduced model fit for purpose of LRT#####

#initializing matrix of maximum likelihood estimate candidates
estimate <- matrix(0,5,16)

#Finding maximum likelihood estimate candidates for 5 sets of starting values	
#Note: EM algorithm not used in this particular application file
for(k in 1:5){
  values <- optim(c(FinalEstimate1[1:11], Initial[k,12:13],FinalEstimate1[14:15]),MaxRegSimp(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=10000))
  estimate[k,] <-  c(values$par, values$value)
}

#Eliminating any maximum likelihood estimate candidates that have infinite likelihood value at parameter estimates
if(any(estimate[,16]==-Inf)) estimate <- estimate[-which(estimate[,16]==-Inf),]

#Choosing the MLEs based on the set of estimates with highest likelihood value
FinalEstimate2 <-  t(estimate[which(estimate[,16]==min(estimate[,16])),])[1:15]

#Attempting continuing maximizations to update proposed MLEs from previous set of MLEs (in this application, with so many parameters, slight improvements were made)
for(k in 1:40){
  print(k)
  FinalEstimate2 <- optim(FinalEstimate2,MaxRegSimp(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=10000))$par
}

#Maximum likelihood value at MLE
Maxi2 <- Maximum(FinalEstimate2[12:13],Xnew,ynew,dead, alive, n, d, c(FinalEstimate2[1:11],FinalEstimate2[14:15],FinalEstimate2[14:15]))


#####Likelihood ratio test for whether censoring is different for each cure status#####
2*(Maxi2-Maxi)


#####Model assuming independent censoring#####

#initializing matrix of maximum likelihood estimate candidates
estimate <- matrix(0,5,14)

#Finding maximum likelihood estimate candidates for 5 sets of starting values	
#Note: EM algorithm not used in this particular application file
for(k in 1:5){
  values <- optim(c(FinalEstimate1[1:11], Initial[k,12:13]),MaxRegInd(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=10000))
  estimate[k,] <-  c(values$par, values$value)
}

#Eliminating any maximum likelihood estimate candidates that have infinite likelihood value at parameter estimates
if(any(estimate[,14]==-Inf)) estimate <- estimate[-which(estimate[,14]==-Inf),]

#Choosing the MLEs based on the set of estimates with highest likelihood value
FinalEstimate3 <-  t(estimate[which(estimate[,14]==min(estimate[,14])),])[1:13]

#Attempting continuing maximizations to update proposed MLEs from previous set of MLEs (in this application, with so many parameters, slight improvements were made)
for(k in 1:40){
  print(k)
  FinalEstimate3 <- optim(FinalEstimate3,MaxRegInd(Xnew,ynew,dead,alive,n,d,2),control=list(maxit=10000))$par
}

#Maximum likelihood value at MLE
Maxi3 <- Maximum2(FinalEstimate3[12:13],Xnew,ynew,dead, alive, n, d, c(FinalEstimate3[1:11]))

#Estimated SEs at MLE 
SE <- sqrt(diag(solve(hessian(Maximum2b, FinalEstimate3, X2=Xnew,y2=ynew,dead=dead,alive=alive,n=n,d=d))))


#####Finding estimates and SEs for Semiparametric models by Peng and Dear, 2000 (PH) and Zhang and Peng, 2007 (AFT)#####
#Note: Standard errors had to be found manually rather than use function's build in option due to numerical difficulties

#Putting data in form necessary for smcure function
datas <- data.frame(y=ynew,delta=c(rep(1,dead),rep(0,n-dead)),Xnew=Xnew)
invisible(capture.output(PH <-smcure(Surv(y, delta) ~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6 , cureform=~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6, data=datas, model = "ph", link="logit", Var=FALSE)))
invisible(capture.output(AFT <-smcure(Surv(y, delta) ~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6 , cureform=~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6, data=datas, model = "aft", link="logit", Var=FALSE)))
print(-PH$b)
print(c(-AFT$b,AFT$beta))

#Obtaining standard errors for both semiparametric PH and AFT cure models using bootstrap 

B <- 1000 # number of bootstraps

ParsPH <- matrix(0,B,9) #Initializing matrix for bootstrap estimates
ParsAFT <- matrix(0,B,10) #Initializing matrix for bootstrap estimates

t <- s <- 0 #counters

#continues looping until B bootstrap data sets obtained; 
#due to occasional numerical failures, some bootstrapped data sets need to be thrown out 
while(t<(B+1) | s<(B+1)){
  print(c(t,s))
  err1 <- err2 <-0	#reset "err1/err2" to 0 (was previously defined as 1 if last bootstrapped data set led to numerical problems
  t <- t+1	#counter increase by one
  s <- s +1	#counter increase by one
  
  #Define Data
  X <- as.matrix(cbind(1,Data[,c(1,2,4,5,6)]))
  y <-  as.vector(Data[,3])
  
  #Choose rows to include in bootstrapped data
  h <- sort(sample(1:n,replace=TRUE))
  
  #Defining covariate and response data based on bootstrapped rows
  X <- X[h,]
  y <- y[h]
  
  #Redefining data to be in proper form & defining required variables
  Xnew <- rbind(X[which(X[,4]==1),], X[which(y==d),],X[which(X[,4]==0),])
  ynew <- c(y[which(X[,4]==1)],y[which(y==d)],y[which(X[,4]==0)])
  dead <- length(y[which(X[,4]==1)])
  alive <- length(y[which(y==d)])
  unknown <- n-dead-alive
  Xnew <- Xnew[,-4]
  
  #Defining dataframe in form needed for smcure function
  datas <- data.frame(y=ynew,delta=c(rep(1,dead),rep(0,n-dead)),Xnew=Xnew)
  
  #tryCatch used where when numerical issues are experiences, err is set to 1 and loop is restarted to try a new bootstrap data set
  tryCatch(invisible(capture.output(smPH <-smcure(Surv(y, delta) ~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6 , cureform=~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6, data=datas, model = "ph", link="logit", Var = FALSE))), error=function(...) err1 <- 1)
  tryCatch(invisible(capture.output(smAFT <-smcure(Surv(y, delta) ~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6 , cureform=~ Xnew.V1 + Xnew.V2 + Xnew.V5 + Xnew.V6, data=datas, model = "aft", link="logit", Var = FALSE))), error=function(...) err2 <- 1)
  
  #Defining estimates for this bootstrapped data set
  if(err1==0 & t<(B+1)) ParsPH[t,] <- c(-smPH$b,smPH$beta)
  if(err2==0 & s<(B+1)) ParsAFT[s,] <- c(-smAFT$b,smAFT$beta)
  
  #a few bootstrap datasets led to estimates that were so bad and were thrown out;  this is very ad hoc (SEs estimates are already poor here)
  if(abs(smPH$b[1])>50 | abs(smPH$b[2])>12 | abs(smPH$b[3])>10 | err1==1) t <- t-1
  if(abs(smAFT$b[1])>50 | abs(smAFT$b[2])>12 | abs(smAFT$b[3])>10 | err2==1) s <- s-1
  
}

#Obtaining bootstrap SEs
apply(ParsPH,2,sd)
apply(ParsAFT,2,sd)



#######################Program Functions#########################

##############SNP Functions#################

#A matrix in SNP fitting for normal 
Anorm <- matrix(c(1,0,1,0,1,0,1,0,3),3,3) 

#Gets "c: coefficients, then "a" coefficients (See Zhang and Davidian 2001)
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
  for(i in 1:(k+1)) s[,i] <- ((T-X%*%c(t[1:5]))/t[6])^(i-1)
  s%*%a
}

#CDF function for standard normal-SNP given "a" values
pSNPnorm <-function(x,a) 1-(a[1]^2*pnorm(x,lower.tail=FALSE)+2*a[1]*a[2]*dnorm(x)+(a[2]^2+2*a[1]*a[3])*(x*dnorm(x)+pnorm(x,lower.tail=FALSE))+2*a[2]*a[3]*(x^2+2)*dnorm(x)+a[3]^2*(x^3*dnorm(x)+3*x*dnorm(x)+3*pnorm(x,lower.tail=FALSE)))

#Density function for normal-SNP given parameter values and data
dSNPtNorm <- function(y,X,d,k,t) {
  fc <- rep(0,length(y))
  a <- acoef(t[(dim(X)[[2]]+2):length(t)],k)
  P_k <- Pk(log(y),X,t,k,a,length(y))
  if(k==1) fc <- 1/(y*t[6])*(P_k)^2*dnorm((log(y)-as.vector(X%*%c(t[1:5])))/t[6])/pSNPnorm((log(d)-as.vector(X%*%c(t[1:5])))/t[6],c(a,0))  
  if(k==2) fc <- 1/(y*t[6])*(P_k)^2*dnorm((log(y)-as.vector(X%*%c(t[1:5])))/t[6])/pSNPnorm((log(d)-as.vector(X%*%c(t[1:5])))/t[6],a)  
  fc
}

#Potentially truncated normal-SNP CDF given data, cure threshold, and parameters
pSNPtNorm <- function(y,X,d,k,t,lower=TRUE) {
  fc <- rep(0,length(y))
  a <- acoef(t[(dim(X)[[2]]+2):length(t)],k) 
  if(k==1) fc <- (pSNPnorm(((log(y)-as.vector(X%*%c(t[1:5])))/t[6]),c(a,0)))/(pSNPnorm(((log(d)-as.vector(X%*%c(t[1:5])))/t[6]),c(a,0)))
  if(k==2) fc <- (pSNPnorm(((log(y)-as.vector(X%*%c(t[1:5])))/t[6]),a))/(pSNPnorm(((log(d)-as.vector(X%*%c(t[1:5])))/t[6]),a))
  if(lower==FALSE) 1-fc else fc
}


#Function to derive set of initial values for maximizing SNP-AFT cure model 
#Note:  Initial values for steps 1-3 in algorithm in paper were obtain using slightly alternative process
#Note:  Input data and "divs=number of grid points in each dimension" and get 5 sets of starting values as matrix
InitVals <- function(X,y,dead, alive,n, d,divs) {
  Starting <- nlminb(c(-14,-2,-0.03,.19,.02,.4,0.5,-1.57,0.018,0.02,1,1,500,1,500),MaxInit(X,y,dead,alive,n,d))$par # optim(c(-14,-2,-0.03,.19,.02,.4,0.5,-1.57,0.018,0.02,1,1,500,1,500),MaxInit(X,y,dead,alive,n,d),control=list(maxit=10000),method="BFGS")$par #tryCatch(nlminb(c(-14,-2,-0.03,.19,.02,.4,0.5,-1.57,0.018,0.02,1,1,500,1,500),MaxInit(X,y,dead,alive,n,d))$par, error=function(...))
  
  PhiVals <- seq(-1.50,1.50,3/(divs-1))
  PhiValsM  <- PhiI <- NULL
  for(i in 1:2) PhiValsM <- cbind(rep(rep(PhiVals,each=divs^(i-1)),divs^(2-i)),PhiValsM)
  Grid <- array(0,rep(divs,2))
  Grid <- array(apply(PhiValsM,1,Maximum, X2=X, y2=y,dead=dead,alive=alive,n=n,d=d,Starting),rep(divs,2))
  
  for(i in 1:5){ g <- as.vector(which(Grid==min(Grid), arr.in=TRUE))
  Grid[(max(1,g[1]-round(divs/5,0)):min(length(PhiVals),g[1]+round(divs/5,0))),(max(1,g[2]-round(divs/5,0)):min(length(PhiVals),g[2]+round(divs/5,0)))] <- max(Grid)
  PhiI <- rbind(PhiI,c(PhiVals[g[2]],PhiVals[g[1]]))
  }
  cbind(matrix(Starting[1:11],5,11,byrow=TRUE),PhiI,matrix(Starting[12:15],5,4,byrow=TRUE))
}

#Likelihood function for normal-SNP; outputs likelihood value at parameters given data
Maximum <- function(Phi,X2,y2,dead,alive,n,d,par){
  par <- c(par[1:11],Phi,par[12:15])	
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:5])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[6:13])*(1-pweibull(y2[1:dead],par[16],par[17])))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y2[(dead+1):(dead+alive)], par[14],par[15]))*(1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
  llike[(dead+alive+1):n] <- -log(dweibull(y2[(dead+alive+1):n],par[14],par[15])*(1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(dweibull(y2[(dead+alive+1):n],par[16],par[17]))*(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[6:13],lower=FALSE))
  return(sum(llike))
}

#Likelihood function for normal-SNP; outputs likelihood value at parameters given data (modified input from previous function)
Maximumb <- function(par,X2,y2,dead,alive,n,d){
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:5])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[6:13])*(1-pweibull(y2[1:dead],par[16],par[17])))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y2[(dead+1):(dead+alive)], par[14],par[15]))*(1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
  llike[(dead+alive+1):n] <- -log(dweibull(y2[(dead+alive+1):n],par[14],par[15])*(1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(dweibull(y2[(dead+alive+1):n],par[16],par[17]))*(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[6:13],lower=FALSE))
  return(sum(llike))
}

#Likelihood function for normal-SNP not assuming dependent censoring; outputs likelihood value at parameters given data 
Maximum2 <- function(Phi,X2,y2,dead,alive,n,d,par){
  par <- c(par[1:11],Phi)	
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:5])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[6:13]))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
  llike[(dead+alive+1):n] <- -log((1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[6:13],lower=FALSE))
  return(sum(llike))
  print(sum(llike))
}

#Likelihood function for normal-SNP not assuming dependent censoring; outputs likelihood value at parameters given data (modified input from previous function)
Maximum2b <- function(par,X2,y2,dead,alive,n,d){	
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1:5])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[6:13]))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
  llike[(dead+alive+1):n] <- -log((1+exp(-X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(1+exp(X2[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[6:13],lower=FALSE))
  return(sum(llike))
  print(sum(llike))
}


#normal-AFT likelihood function for purpose of getting initial values for parameters
MaxInit <- function(X,y,dead,alive,n,d){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:5])))^(-1)*dlnorm(y[1:dead],X[1:dead,]%*%c(par[6:10]),par[11])/plnorm(d,X[1:dead,]%*%c(par[6:10]),par[11])*(1-pweibull(y[1:dead],par[14],par[15])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y[(dead+1):(dead+alive)], par[12], par[13]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
    llike[(dead+alive+1):n] <- -log(dweibull(y[(dead+alive+1):n],par[12],par[13])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(dweibull(y[(dead+alive+1):n],par[14],par[15]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*(1-plnorm(y[(dead+alive+1):n],X[(dead+alive+1):n,]%*%c(par[6:10]),par[11])/plnorm(d,X[(dead+alive+1):n,]%*%c(par[6:10]),par[11])))
    return(sum(llike))
  }
  like
}


#normal-SNP Likelihood function which is used for find MLEs
MaxReg <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:5])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[6:13])*(1-pweibull(y[1:dead],par[16],par[17])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y[(dead+1):(dead+alive)],par[14],par[15]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
    llike[(dead+alive+1):n] <- -log(dweibull(y[(dead+alive+1):n],par[14],par[15])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(dweibull(y[(dead+alive+1):n],par[16],par[17]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[6:13],lower=FALSE))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for find MLEs for reduced model for LRT in paper
MaxRegSimp <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:5])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[6:13])*(1-pweibull(y[1:dead],par[14],par[15])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pweibull(y[(dead+1):(dead+alive)],par[14], par[15]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
    llike[(dead+alive+1):n] <- -log(dweibull(y[(dead+alive+1):n],par[14],par[15])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(dweibull(y[(dead+alive+1):n],par[14],par[15]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[6:13],lower=FALSE))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for find MLEs for model assuming independent censoring and time-to-event variables
MaxRegInd <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1:5])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[6:13]))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1:5])))^(-1))
    llike[(dead+alive+1):n] <- -log((1+exp(-X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)+(1+exp(X[(dead+alive+1):n,]%*%c(par[1:5])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[6:13],lower=FALSE))
    return(sum(llike))
  }
  like
}
