#####################################################################################
#NOTE: PLEASE DO NOT USE CODE FOR ANY PUBLICATION PURPOSES WITHOUT FIRST CONTACTING #
#      PAUL BERNHARDT AT PAUL.BERNHARDT@VILLANOVA.EDU                               #
#####################################################################################

####################Necessary Packages for Simulation######################
library(survival) #For fitting standard survival distributions, creating Survival objects
library(evd) #For generating from Gumbel distribution
library(smcure) #For fitting semiparametric AFT and Cox cure models
#library(actuar)
library(msm) #For generating from truncated normal


##################Simulation Parameters######################
set.seed(8089205)		#seed number
N <- 500			#number of simulated data sets
n <- 300			#number of observations per data set
d <- 1000			#cure threshold (at 1000, this is arbitrarily high for the generated data)
sigma <- 1			#variation in truncated normal to generating events
Gamma <- c(0.5,-0.2)	#parameters in logistic model for probability of cure
Beta <- c(0.3,0.15)	#parameters relating covariates to event in AFT model
EM <- 100 			#number of EM updates (preferrably repeat algorithm to convergence)


#################Simulation Functions###################

##################SNP Functions####################

Anorm <- matrix(c(1,0,1,0,1,0,1,0,3),3,3) #A matrix in SNP fitting for normal 

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
  for(i in 1:(k+1)) s[,i] <- ((T-X%*%c(t[1],t[2]))/t[3])^(i-1)
  s%*%a
}

#CDF function for standard normal-SNP given "a" values
pSNPnorm <-function(x,a) 1-(a[1]^2*pnorm(x,lower.tail=FALSE)+2*a[1]*a[2]*dnorm(x)+(a[2]^2+2*a[1]*a[3])*(x*dnorm(x)+pnorm(x,lower.tail=FALSE))+2*a[2]*a[3]*(x^2+2)*dnorm(x)+a[3]^2*(x^3*dnorm(x)+3*x*dnorm(x)+3*pnorm(x,lower.tail=FALSE)))

#Density function for normal-SNP given parameter values and data
dSNPtNorm <- function(y,X,d,k,t) {
  fc <- rep(0,length(y))
  a <- acoef(t[4:length(t)],k)
  P_k <- Pk(log(y),X,t,k,a,length(y))
  if(k==1) fc <- 1/(y*t[3])*(P_k)^2*dnorm((log(y)-as.vector(X%*%c(t[1],t[2])))/t[3])/pSNPnorm((log(d)-as.vector(X%*%c(t[1],t[2])))/t[3],c(a,0))  
  if(k==2) fc <- 1/(y*t[3])*(P_k)^2*dnorm((log(y)-as.vector(X%*%c(t[1],t[2])))/t[3])/pSNPnorm((log(d)-as.vector(X%*%c(t[1],t[2])))/t[3],a)  
  fc
}

#Potentially truncated normal-SNP CDF given data, cure threshold, and parameters
pSNPtNorm <- function(y,X,d,k,t,lower=TRUE) {
  fc <- rep(0,length(y))
  a <- acoef(t[4:length(t)],k) 
  if(k==1) fc <- (pSNPnorm(((log(y)-as.vector(X%*%c(t[1],t[2])))/t[3]),c(a,0)))/(pSNPnorm(((log(d)-as.vector(X%*%c(t[1],t[2])))/t[3]),c(a,0)))
  if(k==2) fc <- (pSNPnorm(((log(y)-as.vector(X%*%c(t[1],t[2])))/t[3]),a))/(pSNPnorm(((log(d)-as.vector(X%*%c(t[1],t[2])))/t[3]),a))
  if(lower==FALSE) 1-fc else fc
}


##################Functions for Maximization##################

#normal-SNP Likelihood function which is used for finding MLEs
MaxReg <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1],par[2])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[3:7])*(1-pgamma(y[1:dead],par[10],par[11])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pgamma(y[(dead+1):(dead+alive)],par[8],par[9]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1],par[2])))^(-1))
    llike[(dead+alive+1):n] <- -log(dgamma(y[(dead+alive+1):n],par[8],par[9])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)+(dgamma(y[(dead+alive+1):n],par[10],par[11]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[3:7],lower=FALSE))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for finding MLEs in reduced model (for LRT)
MaxRegSimp <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1],par[2])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[3:7])*(1-pgamma(y[1:dead],par[8],par[9])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pgamma(y[(dead+1):(dead+alive)],par[8],par[9]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1],par[2])))^(-1))
    llike[(dead+alive+1):n] <- -log(dgamma(y[(dead+alive+1):n],par[8],par[9])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)+(dgamma(y[(dead+alive+1):n],par[8],par[9]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[3:7],lower=FALSE))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function for finding MSEs assuming indep
MaxRegInd <- function(X,y,dead,alive,n,d,k){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1],par[2])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[3:7]))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1],par[2])))^(-1))
    llike[(dead+alive+1):n] <- -log((1+exp(-X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)+(1+exp(X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[3:7],lower=FALSE))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for finding MLEs based on EM algorithm
MaxEM <- function(X,y,dead,alive,n,d,k,PARS){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1],par[2])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[3:7])*(1-pgamma(y[1:dead],par[10],par[11])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pgamma(y[(dead+1):(dead+alive)], par[8], par[9]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1],par[2])))^(-1))
    llike[(dead+alive+1):n] <- -1/((1+exp(-X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[8],PARS[9])+(1+exp(X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[10],PARS[11])*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,PARS[3:7],lower=FALSE))*((1+exp(-X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[8],PARS[9])*log((1+exp(-X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*dgamma(y[(dead+alive+1):n],par[8],par[9]))+(1+exp(X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[10],PARS[11])*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,PARS[3:7],lower=FALSE)*log((1+exp(X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*dgamma(y[(dead+alive+1):n],par[10],par[11])*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[3:7],lower=FALSE)))
    return(sum(llike))
  }
  like
}

#normal-SNP Likelihood function which is used for finding MLEs based on EM algorithm in reduced model (for LRT)
MaxEM2 <- function(X,y,dead,alive,n,d,k,PARS){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1],par[2])))^(-1)*dSNPtNorm(y[1:dead],X[1:dead,],d,k,par[3:7])*(1-pgamma(y[1:dead],par[8],par[9])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pgamma(y[(dead+1):(dead+alive)], par[8], par[9]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1],par[2])))^(-1))
    llike[(dead+alive+1):n] <- -1/((1+exp(-X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[8],PARS[9])+(1+exp(X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[8],PARS[9])*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,PARS[3:7],lower=FALSE))*((1+exp(-X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[8],PARS[9])*log((1+exp(-X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*dgamma(y[(dead+alive+1):n],par[8],par[9]))+(1+exp(X[(dead+alive+1):n,]%*%c(PARS[1],PARS[2])))^(-1)*dgamma(y[(dead+alive+1):n],PARS[8],PARS[9])*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,PARS[3:7],lower=FALSE)*log((1+exp(X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*dgamma(y[(dead+alive+1):n],par[8],par[9])*pSNPtNorm(y[(dead+alive+1):n],X[(dead+alive+1):n,],d,k,par[3:7],lower=FALSE)))
    return(sum(llike))
  }
  like
}


#Function to derive set of initial values for maximizing SNP-AFT cure model 
#Note:  Initial values for steps 1-3 in algorithm in paper were obtain using slightly alternative process
#Note:  Input data and "divs=number of grid points in each dimension" and get 5 sets of starting values as matrix
InitVals <- function(X,y,dead, alive,n, d,divs) {
  delta <- c(rep(1,dead),rep(0,n-dead))
  datas <- data.frame(y=y,delta=delta,Xnew=X)
  invisible(capture.output(sm <-smcure(Surv(y, delta) ~ Xnew.2 , cureform=~Xnew.2, data=datas, model = "aft", link="logit", Var = FALSE)))
  
  Starting <- tryCatch(nlminb(c(-sm$b,sm$beta,1,3,0.15,3,0.3),MaxInit(X,y,dead,alive,n,d))$par, error=function(...) optim(c(0.5,-0.2,0.3,0.15,1,3,0.15,3,0.15),MaxInit(X,y,dead,alive,n,d))$par)
  
  PhiVals <- seq(-1.50,1.50,3/(divs-1))
  PhiValsM  <- PhiI <- NULL
  for(i in 1:2) PhiValsM <- cbind(rep(rep(PhiVals,each=divs^(i-1)),divs^(2-i)),PhiValsM)
  Grid <- array(0,rep(divs,2))
  Grid <- array(apply(PhiValsM,1,Maximum, X2=X, y2=y,dead=dead,alive=alive,n=n,d=d,Starting),rep(divs,2))
  
  for(i in 1:5){ g <- as.vector(which(Grid==min(Grid), arr.in=TRUE))
  Grid[(max(1,g[1]-round(divs/5,0)):min(length(PhiVals),g[1]+round(divs/5,0))),(max(1,g[2]-round(divs/5,0)):min(length(PhiVals),g[2]+round(divs/5,0)))] <- max(Grid)
  PhiI <- rbind(PhiI,c(PhiVals[g[2]],PhiVals[g[1]]))
  }
  cbind(matrix(Starting[1:5],5,5,byrow=TRUE),PhiI,matrix(Starting[6:9],5,4,byrow=TRUE))
}

#Likelihood function for normal-SNP; outputs likelihood value at parameters given data
Maximum <- function(Phi,X2,y2,dead,alive,n,d,par){
  par <- c(par[1:5],Phi,par[6:9])	
  llike <- rep(0,n)
  llike[1:dead] <- -log((1+exp(X2[1:dead,]%*%c(par[1],par[2])))^(-1)*dSNPtNorm(y2[1:dead],X2[1:dead,],d,2,par[3:7])*(1-pgamma(y2[1:dead],par[10], par[11])))
  if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pgamma(y2[(dead+1):(dead+alive)], par[8], par[9]))*(1+exp(-X2[(dead+1):(dead+alive),]%*%c(par[1],par[2])))^(-1))
  llike[(dead+alive+1):n] <- -log(dgamma(y2[(dead+alive+1):n],par[8],par[9])*(1+exp(-X2[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)+(dgamma(y2[(dead+alive+1):n],par[10], par[11]))*(1+exp(X2[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*pSNPtNorm(y2[(dead+alive+1):n],X2[(dead+alive+1):n,],d,2,par[3:7],lower=FALSE))
  return(sum(llike))
}

#normal-AFT likelihood function for purpose of getting initial values for parameters
MaxInit <- function(X,y,dead,alive,n,d){
  like <- function(par){
    llike <- rep(0,n)
    llike[1:dead] <- -log((1+exp(X[1:dead,]%*%c(par[1],par[2])))^(-1)*dlnorm(y[1:dead],X[1:dead,]%*%c(par[3],par[4]),par[5])/plnorm(d,X[1:dead,]%*%c(par[3],par[4]),par[5])*(1-pgamma(y[1:dead],par[8],par[9])))
    if(alive>0) llike[(dead+1):(dead+alive)] <- -log((1-pgamma(y[(dead+1):(dead+alive)], par[6], par[7]))*(1+exp(-X[(dead+1):(dead+alive),]%*%c(par[1],par[2])))^(-1))
    llike[(dead+alive+1):n] <- -log(dgamma(y[(dead+alive+1):n],par[6],par[7])*(1+exp(-X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)+(dgamma(y[(dead+alive+1):n],par[8],par[9]))*(1+exp(X[(dead+alive+1):n,]%*%c(par[1],par[2])))^(-1)*(1-plnorm(y[(dead+alive+1):n],X[(dead+alive+1):n,]%*%c(par[3],par[4]),par[5])/plnorm(d,X[(dead+alive+1):n,]%*%c(par[3],par[4]),par[5])))
    return(sum(llike))
  }
  like
}


###########################Simulation###################################


#initialized matrices for parameter estimates
SNPAFT <- matrix(0, N,11)
SNPAFT2 <- matrix(0,N,9)

#intialized vectors
Maxi <- Maxi2 <- rep(0,N)

#looping through all datasets
for(i in 1:N){
  print(i)
  
  X <- cbind(1,rnorm(n,6,4))	#generate covariate
  p <- rbinom(n,1,(1+exp(-as.vector(X%*%Gamma)))^(-1))	#generate cure status
  logt <- rtnorm(n, mean=as.vector(X%*%Beta), sd=sigma,upper=log(d))	#generate log time-to-event
  y <-  exp(logt)	#transforming to time-of-event
  
  #censoring distributions;  can be made different if desired
  a <- rgamma(n,3,0.15)
  b <- rgamma(n,3,0.15)
  
  #Mixture of censoring distributions (distribution a if uncured, b if cured)
  cens <- a*(1-p)+b*p
  
  #Defining those who survived to cure threshold (assuming no observations beyond that time)
  y[which(p==1)] <- d
  
  #Replacing time-to-event with censored value if censored value is less than time-to-event
  y[which(y>cens)] <- cens[which(y>cens)]
  
  #Defining matrix of covariates and vector of outcomes
  Xnew <- rbind(X[which(y!=d & y!=cens),], X[which(y==d),],X[which(y==cens),])
  ynew <- c(y[which(y!=d & y!=cens)],y[which(y==d)],y[which(y==cens)])
  
  #number who observed uncured (dead), observed cured (alive), censored below threshold (unknown)
  dead <- length(y[which(y!=d & y!=cens)])
  alive <- length(y[which(y==d)])
  unknown <- n-dead-alive
  
  #Obtaining starting values for ML methods
  Initial <- InitVals(Xnew,ynew,dead,alive,n,d,9)
  
  
  ######EM Algorithm Method for obtaining parameter estimates in proposed SNP-AFT model#####
  estvals <- matrix(0,5,11)
  Maxvals <- rep(0,5)
  for(j in 1:5){
    estimates <- matrix(0,EM,11)
    estimates[1,] <- Initial[j,]
    ##Looping through EM updates of EM algorithm
    for(k in 2:EM){
      estimates[k,] <- nlm(MaxEM(Xnew,ynew,dead,alive,n,d,2,estimates[(k-1),]),estimates[(k-1),])$estimate
    }
    estvals[j,] <- estimates[EM,]
    Maxvals[j] <- Maximum(estvals[j,6:7],Xnew,ynew,dead, alive, n, d, c(estvals[j,1:5],estvals[j,8:11]))
  }
  SNPAFT[i,] <- estvals[which(Maxvals==min(Maxvals)),]	#Definining estimate for current data set
  
  
  ######Alternatively running N-R maximum likelihood methods (or similar vein of algorithms) to get SNP-AFT MLEs######
  #estvals <- matrix(0,5,12) #initialing matrix of proposed estimates
  
  #Looping through 5 initial data sets to get potential ML values
  #for(k in 1:5){
  #	values <- tryCatch(optim(Initial[k,],MaxReg(Xnew,ynew,dead,alive,n,d,2)), error=function(...) nlminb(Initial[k,],MaxReg(Xnew,ynew,dead,alive,n,d,2)))
  #	estvals[k,] <-  c(values$par, c(values$objective,values$value))
  #}
  
  #Removing any set of estimates from possible MLE list if its maximized likelihood is "infinite"
  #if(any(estvals[,12]==-Inf)) estvals <- estvals[-which(estvals[,12]==-Inf),]
  
  #Defining estimate for current data set
  #SNPAFT[i,] <-  t(estvals[which(estvals[,12]==min(estvals[,12])),])[1:11]
  
  
  #Maximized likelihood value at parameter estimates
  Maxi[i] <- Maximum(SNPAFT[i,6:7],Xnew,ynew,dead, alive, n, d, c(SNPAFT[i,1:5],SNPAFT[i,8:11]))
  
  
  #######EM Algorithm Method for obtaining parameter estimates in reduced SNP-AFT model (for obtaining LRT value)######
  estvals <- matrix(0,5,9)
  Maxvals <- rep(0,5)
  for(j in 1:5){
    estimates <- matrix(0,EM,9)
    estimates[1,] <- Initial[j,1:9]
    ##Looping through EM updates of EM algorithm
    for(k in 2:EM){
      estimates[k,] <- nlm(MaxEM2(Xnew,ynew,dead,alive,n,d,2,estimates[(k-1),]),estimates[(k-1),])$estimate
    }
    estvals[j,] <- estimates[EM,]
    Maxvals[j] <- Maximum(estvals[j,6:7],Xnew,ynew,dead, alive, n, d, c(estvals[j,1:5],estvals[j,8:9],estvals[j,8:9]))
  }
  SNPAFT2[i,] <- estvals[which(Maxvals==min(Maxvals)),]	#Definining estimate for current data set
  
  
  ######Alternatively running N-R maximum likelihood methods (or similar vein of algorithms) to get SNP-AFT MLEs######
  #estvals <- matrix(0,5,10) #initialing matrix of proposed estimates
  
  #Looping through 5 initial data sets to get potential ML values
  #for(k in 1:5){
  #	values <- tryCatch(optim(Initial[k,1:9],MaxReg(Xnew,ynew,dead,alive,n,d,2)), error=function(...) nlminb(Initial[k,1:9],MaxReg(Xnew,ynew,dead,alive,n,d,2)))
  #	estvals[k,] <-  c(values$par, c(values$objective,values$value))
  #}
  
  #Removing any set of estimates from possible MLE list if its maximized likelihood is "infinite"
  #if(any(estvals[,10]==-Inf)) estvals <- estvals[-which(estvals[,10]==-Inf),]
  
  #Defining estimate for current data set
  #SNPAFT2[i,] <-  t(estvals[which(estvals[,10]==min(estvals[,10])),])[1:9]
  
  
  #Maximized likelihood value at parameter estimates
  Maxi2[i] <- Maximum(SNPAFT2[i,6:7],Xnew,ynew,dead, alive, n, d, c(SNPAFT2[i,1:5],SNPAFT2[i,8:9],SNPAFT2[i,8:9]))
}
}

LRTstats <- 2*(Maxi2-Maxi)



