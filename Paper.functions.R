# Copyright 2017, Andrey Skripnikov, All rights reserved.

rm(list=ls())
library(plotrix)
library(glmnet)
library(MASS)
library(glasso)
library(methods)


##########################################
## Creates block-diagonal matrix       ###
## with A.list[[i]] as its i-th block  ###
##########################################

block.diag <- function(A.list){
  p <- nrow(A.list[[1]])
  m <- ncol(A.list[[1]])
  K <- length(A.list)
  
  block_matrix <- matrix(0,K*p,K*m)
  for (i in 1:K){
    block_matrix[1:p + (i-1)*p,1:m + (i-1)*m] <- A.list[[i]]
  }  
  return(block_matrix)
}




#####################
### GENERATES TRANSITION MATRIX FOR THE TIME SERIES
###
### Function provided by Sumanta Basu
#####################

gen_A = function(  # returns a list of d matrices A_1, ..., A_d
  p,    # generate p x p VAR
  d=1,   # generate VAR(d)
  max_eig,   # spectral norm of A
  edge_density = 0.1, # if different for different lags, provide a vector of length d
  nonzero_diag = 1,  # ensure all diagonal entries are non-zero; if different for different lags, provide a d-vector of 0 and 1
  stationary = 1,  # ensure stationarity
  network.family = "random",  # A matrix filled randomly or with specific structure
  structure = NULL  # a pre-determined pattern of non-zero elements in generated matrix
){
  
  A = list()
  if (is.null(structure)){
    for (lag in 1:d){
      e_d = ceiling(p^2 * ifelse(length(edge_density) == 1, edge_density, edge_density[lag]))
      temp = sample(p^2, e_d)
      temp_v = rep(0, p^2); temp_v[temp] = 1
      A[[lag]] = array(temp_v, c(p,p))
      if (ifelse(length(nonzero_diag) == 1, nonzero_diag, nonzero_diag[lag]))
        diag(A[[lag]]) = 1
    }
  }
  if(!is.null(structure)){
    for (lag in 1:d){
      A[[lag]]=structure
    }
  }
  
  if (stationary){
    A.big = array(0, c(p*d, p*d))
    
    for (i in 1:d)
      A.big[1:p, ((i-1)*p+1):(i*p)] = max_eig*A[[i]]
    
    if (d > 1)
      diag(A.big[(p+1):(p*d), 1:((d-1)*p)]) = 1
    
    temp = max(abs(eigen(A.big)$values))
    count = 0
    while (temp > max_eig){
      count = count+1
      #  print(paste("count:", count, "max_eigen_value:", round(temp, 2), "signal:", round(max(abs(A.big[1:p,])), 2)))
      A.big[1:p,] = A.big[1:p,]*0.95
      temp = max(abs(eigen(A.big)$values))
    }
    
    for (i in 1:d)
      A[[i]] = A.big[1:p, ((i-1)*p+1):(i*p)]
    #  print(paste("signal reduced to", round(max(abs(A.big[1:p,])), 2)))
    
  }
  return(list(A=A,
              Signal=round(max(abs(A.big[1:p,])), 2))
  )
}




####################################
####  Makes a list of matrices to be transition matrices
####  of a stationary time series
####################################

make.stationary <- function(#returns the stationary matrix
  #        the signal of the matrix
  A,             #list of matrices to be made stationary
  max_eig,   #max matrix eigenvalue allowed
  d=1            #number of matrices in the list(VAR(d))
){
  
  A.big = array(0, c(p*d, p*d))
  
  for (i in 1:d)
    A.big[1:p, ((i-1)*p+1):(i*p)] = A
  
  temp = max(abs(eigen(A.big)$values))
  count = 0
  
  while (temp > max_eig){
    count = count+1
    # print(paste("count:", count, "max_eigen_value:", round(temp, 2), "signal:", round(max(abs(A.big[1:p,])), 2)))
    A.big[1:p,] = A.big[1:p,]*0.95
    temp = max(abs(eigen(A.big)$values))
  }
  
  for (i in 1:d)
    A = A.big[1:p, ((i-1)*p+1):(i*p)] 
  
  return(list(A=A,
              Signal=round(max(abs(A.big[1:p,])), 2)))
}



################
##### Function generating transition matrices for related VAR(1) models
################

A.setup <- function( # returns two transition matrices,
  p,          # number of variables per entity
  thresh,     # non-zero matrix elements should be greater than 'thresh'
  diff=FALSE, # "different" - whether to generate identical or different matrices
  ed=0.05,    #  edge density - proportion of non-zero off-diagonal elements
  comm=0.01,  #  for the case of different matrices - proportion of off-diagonal elements that don't share same position in two matrices
  max_eig #  maximum eigenvalue of transition matrices
){
  
  ### For the case of two different transition matrices:
  ###    - if comm=0: we generate them separately;
  ###    - if comm>0: we generate a common matrix and then randomly add
  ###                 non-zero elements randomly to both transition matrices
  ###                 (that way there will be some common off-diagonal zeros as well)
  
  if (diff==TRUE){  
    if (comm == 0){
      A.true <- list()
      for (i in 1:K){
        n.iter <- 0
        repeat{
          n.iter <- n.iter + 1
          A.obj <- gen_A(p,ed=ed,max_eig=max_eig)
          Signal <- A.obj$Signal
          if (Signal > thresh) break; 
        }
        A.true[[i]] <- A.obj$A[[1]];
      }
    }
    
    if (comm != 0){   
      repeat{
        A.Gen <- gen_A(p,edge_density=ed,max_eig=max_eig)
        if (A.Gen$Signal > thresh) break;
      }
      
      A.true <- list()
      
      for(i in 1:K){
        repeat{
          Ahet.Gen <- gen_A(p,edge_density=comm,max_eig=max_eig)
          A.Obj <- make.stationary(pmax(Ahet.Gen$A[[1]],A.Gen$A[[1]]),max_eig=max_eig)
          if (A.Obj$Signal > thresh) break;
        }
        A.true[[i]] <- A.Obj$A;
      }
    }
    return(A.true=A.true)
  }
  
  ### If we want identical matrices: just generate a matrix once
  ### and copy it for the second entity.
  
  if (diff == FALSE){ 
    A.true <- list()
    repeat{
      A.Gen <- gen_A(p,edge_density=ed,max_eig=max_eig)
      if (A.Gen$Signal > thresh) break;
    }  
    A.true[[1]] <- A.Gen$A[[1]]
    for (i in 2:K){
      A.true[[i]] <- A.true[[1]]
    }
    return(A.true=A.true)
  }
}


##################
### GENERATES DATA FROM A PRE-SPECIFIED VAR PROCESS:
### X_t = A[[1]]X_t-1 + ... + A[[d]]X_t-d + eps, eps ~ N(O,Sigma_error)
###
### Function provided by Sumanta Basu
###################

require(MASS)
gen_dat = function(  	# returns a p x T matrix of observations {X^1, ..., X^T}
  T = NULL,
  error_sd = NULL, # optional, a p-dim vector of the individual sdevs of the p time series
  Sigma_error = NULL, # input a p x p correlation matrix; otherwise it is set to diag(error_sd) * identity
  SNR = 2,	# signal-to-noise ratio, used to determine error_sd (= abs(A[i,j])/SNR)
  A = NULL,	# a list of d adjacency matrices, as returned by gen_A
  cut_T = 500*length(A) # time to wait before process reaches stability 
){
  d = length(A)
  p = dim(A[[1]])[1]
  X = array(0, c(p, T+cut_T))
  if(is.null(error_sd))
    error_sd = rep(max(abs(A[[1]]))/SNR, p)
  else if (length(error_sd) == 1)
    error_sd = rep(error_sd, p)
  
  if (is.null(Sigma_error))
    Sigma_error = diag(p)
  Sigma_error = diag(error_sd) %*% Sigma_error %*% diag(error_sd); #print(round(Sigma_error, 4))
  
  X = t(mvrnorm(n = T+cut_T, mu = rep(0, p), Sigma = Sigma_error))
  for (tt in (d+1):(T+cut_T))
    for (lg in 1:d)
      X[,tt] = X[,tt]+A[[lg]] %*% X[,tt-lg]
  
  return(X[,-seq(cut_T)])
}


########################
### Function that generates error covariance matrix that
### follows a L-factor model with a diagonal Sigma_U matrix
### Generation approach was taken from Fan et al., 2011
########################

Sigma.Gen <- function(# returns error covariance matrix Sigma,
                      #         inverse of that matrix,
                      #         Sigma_U estimate
                      p, # number of variables
                      t, # number of time points
                      L, # number of factors
                      Sigma_Ux=NULL # Sigma_U matrix
                      ){
  
  repeat{
    if (is.null(Sigma_Ux)){
      Sigma_Ux.Inv <- diag(runif(p,1,3),p)
      Sigma_Ux <- solve(Sigma_Ux.Inv)
    }
    
    U <- mvrnorm(t,rep(0,p),Sigma_Ux)
    
    f <- matrix(rnorm(t*L,0,1),t,L)
    b <- matrix(rnorm(p*L,0,1),p,L)
    
    Common <- b %*% t(b)
    SnR <- sum(eigen(Common)$values)/sum(eigen(Sigma_Ux)$values)  
    if ((SnR >= 2)) break;  
  }
  Sigma <- b %*% t(b) + Sigma_Ux
  Sigma.Inv <- solve(Sigma)
  
  return(list(Sigma=Sigma,
              Sigma.Inv=Sigma.Inv,
              Sigma_Ux=Sigma_Ux))
}


#################
#### Uses L-factor model to calculate error covariance and its inverse
#### Input - time series of error vectors
#################

factor.glasso <- function(# returns estimate of error covariance,
                          #         estimate of inverse error covariance,
                          #         estimate of Sigma_U,
                          #         number of factors,
                          #         sample error covariance
                          eps,     # error vector
                          lambda,  # graphical lasso parameter for Sigma_U estimation
                          L        # in case we know number L of factos in advance
                          ){
  
  p <- nrow(eps)
  Cov.Y <- cov(t(eps))
  
  if (L==0){
    L.obj <- L.est(Cov.Y,p-1)
    L.hat <- L.obj$L1
    TVE <- L.obj$TVE
  } else {eig <- eigen(Cov.Y)$values; L.hat <- L; TVE <- sum(eig[1:L.hat])/sum(eig)}
  
  SVD <- eigen(Cov.Y)
  Phi_l <- diag(SVD$values[1:L.hat],L.hat)
  Lambda_est <- SVD$vectors[,1:L.hat] %*% sqrt(Phi_l)
  Sigma_Y <- Lambda_est %*% t(Lambda_est)
  Sigma_Ux_est <- Cov.Y - Sigma_Y
  Sigma_Ux_est1 <- nearPD(Sigma_Ux_est)$mat
  
  gl <- glasso(as.matrix(Sigma_Ux_est1),lambda)
  Sigma_Ux_est2 <- gl$w
  Sigma_Ux_est2.Inv <- gl$wi
  
  Sigma.est <- Sigma_Y + Sigma_Ux_est2
  Sigma.est.Inv <- solve(Sigma_Y + Sigma_Ux_est2)
  
  return(list(Sigma.est = Sigma.est,
              Sigma.est.Inv = Sigma.est.Inv,
              Sigma_Ux.est = Sigma_Ux_est2,
              L=L.hat,
              Cov.Mat=Cov.Y,
              TVE=TVE))
}


########################
########### Function that:
###########  - uses the sparse estimate of transition matrix to calculate error vector;
###########  - uses factor models to calculate error covariance estimate via factor.glasso() function
#########################

L1.1step <- function(# returns estimate of error covariance,
                     #         estimate of inverse error covariance,
                     #         estimate of Sigma_U,
                     #         number of factors,
                     #         sample error covariance,
                     Data, # Data matrix
                     L=0,  # number of factors(if known in advance)
                     A.L1.est # sparse transition matrix estimate
                     ){ 
  
  t <- ncol(Data)
  p <- nrow(Data)

  Response <- matrix(c(t(Data[,-1])))
  B <- matrix(Data[,-t],t-1,p,byrow=TRUE)
  X <- diag(1,p) %x% B

  ## RESIDUALS, DOING FACTOR MODELS, SVD+GLASSO 
  
  eps.L1 <- Response - X %*% (A.L1.est)
  eps.L1 <- matrix(eps.L1,p,t-1,byrow=TRUE)  
  eps_c.L1 <- t(scale(t(eps.L1),scale=FALSE))
  
  diag.est <- apply(eps_c.L1,1,sd)
  
  Fac.glasso <- factor.glasso(eps_c.L1,lambda=sqrt(log(p)/t),L=L)
  Sigma.est <- Fac.glasso$Sigma.est
  Sigma.est.Inv <- Fac.glasso$Sigma.est.Inv
  Sigma_Ux.est <- Fac.glasso$Sigma_Ux.est
  Cov.Mat <- Fac.glasso$Cov.Mat
  L <- Fac.glasso$L
  TVE <- Fac.glasso$TVE
  
  return(list(Sigma.est=Sigma.est,
              Sigma.est.Inv=Sigma.est.Inv,
              Sigma_Ux.est=Sigma_Ux.est,
              L=L,
              Cov.Mat=Cov.Mat,
              diag.est=diag.est,
              TVE=TVE))
}


########################
### Factor estimation function:
### figures out number of common underlying factors for error covariance matrix
### for one entity.
########################

L.est <- function(Cov.y,M){

  criterion <- rep(0,M)

  SVD <- eigen(Cov.y)
  eig <- SVD$values
  eig.total <- sum(eig)
  
  eig.med <- mean(eig)
  for(L1 in 1:M){
       if (eig[L1] < eig.med){
         TVE = (sum(eig[1:(L1-1)])/sum(eig))
         break;
       }
     }

  return(list(L1=L1,
              TVE=TVE))
}




######################
##### Function performs our estimation procedure for one entity
##### (takes INVERSES OF SIGMA as input)
######################

sep.lasso.Inv <- function(# returns set of estimates,
                          # data matrix,
                          # response vector,
                          # lambda.path - the solution path produced by glmnet
                          Data,          # Data matrix
                          Sigma.Inv  # Inverse error covariance matrix
                          ){
  
  Data.P <- Data.Prep.Sep.Inv(Data,Sigma.Inv) 
  
  Glmnet.obj <- glmnet(Data.P$X,
                       Data.P$Y,
                       family="gaussian",
                       intercept=intercept,
                       standardize=standardize)
  
  beta <- as.matrix(coef(Glmnet.obj)[-1,])
    
  ### We exclude replicates that provide solution paths with highest density 
  ### of estimate being (fail.thresh2*100%) or lower, which is not detailed enough to select from
  if (max(Glmnet.obj$df) <= (fail.thresh*(p^2))) fail <<- 1
    
    return(list(Est=beta,
                X=Data.P$X,
                Y=Data.P$Y,
                lambda.path=Glmnet.obj$lambda))
}



##############
#### Function calculates all the performance measurements of the estimates:
####             FP,FN,TP,TN,Frobenius difference
##############

Measures.Vec <- function(# returns performance measurements of estimates
                         A.est,   # estimates
                         A.true   # true matrix
                         ){ 
  
  FP <- sum(!!A.est  & !A.true)/(sum(!A.true))
  FN <- sum(!A.est  & !!A.true)/(sum(!!A.true))
  TN <- 1 - FP
  TP <- 1 - FN
  Frob <- norm(A.est - A.true,type="F")/norm(A.true,type="F")
  return(list(FP=FP,
              FN=FN,
              TP=TP,
              TN=TN,
              Frob=Frob))
}



############################################################################
#### CALCULATES STANDARDIZED FROBENIUS DIFFERENCE BETWEEN TWO MATRICES  ####
############################################################################

Frob.comp <- function(Sigma.est,Sigma.true){   
  return(norm(Sigma.est - Sigma.true,type="F")/norm(Sigma.true,type="F")) 
}



#################
### Function takes square root of a matrix
##  Take mA as input, returns B: B'B = mA
#################

fnMatSqrt <- function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = sqrt(d)
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}



#############
### Function takes inverse of a matrix
##  Take mA as input, returns mA^(-1)
#############

fnMatInverse <- function(mA) {
  ei = eigen(mA)
  d = ei$values
  d = (d+abs(d))/2
  d2 = 1/d
  d2[d == 0] = 0
  return(ei$vectors %*% diag(d2) %*% t(ei$vectors))
}



################
### Soft-thresholding operator
### a - argument
### b - threshold
################

softthresh <- function(a,b){
  r <- abs(a) - b
  return(ifelse(r>=0,sign(a)*r,0))
}

#######################################
#### Hard-thresholding function #######
#######################################

sparsify <- function(m,a){   
  m1 <- ifelse(abs(m)<a,0,m)
  return(m1)  
}



####################
#### Converts a vectorized matrix
#### (matrix stretched into a 1-dimensional vector)
#### back into its original form
####################


ConvertToMatrix.Full <- function(# returns a list of matrices(one matrix per entity)
  vec,        #vectorized matrix
  p           #number of variables
){
  K <- length(vec)/p^2
  A <- list()
  for (i in 1:K){
    A[[i]] <- matrix(vec[(1+(i-1)*(p^2)):(i*p^2)],p,p,byrow=TRUE)
  }
  return(A)
}




######################################
### Function prepares response vector and design matrix
### for a standard regression setup for SINGLE model(one entity)
### (takes INVERSE OF SIGMA as input)
######################################

Data.Prep.Sep.Inv <- function(# returns response vector,
                              #         design matrix,
                              #         number of variables and time points
                              Data, # time series matrix, rows - variables, columns - time points
                              Sigma.Inv # inverse error covariance matrix
){   
  p <- nrow(Data)
  t <- ncol(Data)
  
  Sigma.Sqrt.Inv.new <- fnMatSqrt(Sigma.Inv) %x% diag(1,t-1)
 # print(sum(abs(Sigma.Sqrt.Inv1.new - fnMatSqrt(Sigma.Inv) %x% fnMatSqrt(diag(1,t-1)))))
 # Sigma.Sqrt.Inv1.new <- sqrtm(Sigma.Inv %x% diag(1,t-1))
 # Entity <- Data[1:p,]
 
  ### Setting up all the matrices for regression problem, accounting for Sigma estimate
  Response <- matrix(c(t(Data[,-1])))
  X <- matrix(Data[,-t],t-1,p,byrow=TRUE)
  X_Sigma <- optim.product(Sigma.Sqrt.Inv.new,X,p)
  Response_Sigma <- Sigma.Sqrt.Inv.new %*% Response
  
  X.full <- X_Sigma
  Y.full <- Response_Sigma
  
  return(list(X=X.full,
              Y=Y.full,
              p=p,
              t=t))
}


############################
### Function prepares response vector and design matrix
### for a standard regression setup for JOINT model(multiple entities)
############################


Data.Prep.Inv <- function(# returns response vector,
                          #         design matrix,
                          #         number of variables and time points
  Data, # time series matrix, rows - variables, columns - time points
  Sigma.Inv
){  
  
  ### K1 - number of entities
  if (is.list(Sigma.Inv)) K1 <- length(Sigma.Inv)
  if (!is.list(Sigma.Inv)) K1 <- 1
  
  p <- nrow(Data)/K1
  t <- ncol(Data)
  X.full <- list()
  Y.full <- NULL
  
  for (k in 1:K1){
    Sigma.Sqrt.Inv.new <- fnMatSqrt(Sigma.Inv[[k]]) %x% diag(1,t-1)
    Entity <- Data[(k-1)*p + 1:p,]
    Response_Entity <- matrix(c(t(Entity[,-1])))
    X_Entity <- matrix(Entity[,-t],t-1,p,byrow=TRUE)
    X_Entity_Sigma <- optim.product(Sigma.Sqrt.Inv.new,X_Entity,p)
    Response_Entity_Sigma <- Sigma.Sqrt.Inv.new %*% Response_Entity
    
    X.full[[k]] <- X_Entity_Sigma
    Y.full <- c(Y.full,Response_Entity_Sigma)
  }

  X.full <- block.diag(X.full)

  return(list(X=X.full,
              Y=Y.full,
              p=p,
              t=t,
              connList=connList))
}



############################
### Function prepares response vector and design matrix
### for a standard regression setup for GENERAL NUMBER OF ENTITIES K
############################
# 
# Data.Prep.Inv.Full <- function(# returns response vector,
#   #         design matrix,
#   #         number of variables and time points
#   Data, # time series matrix, rows - variables, columns - time points
#   Sigma.Inv# list of Sigma Inverse matrices
# ){  
#   K <- length(Sigma.Inv)
#   p <- nrow(Data)/K
#   t <- ncol(Data)
#   
#   X_Sigma <- list()
#   Response_Sigma <- list()
#   
#   for (i in 1:K){
#     Sigma.Sqrt.Inv1.new <- fnMatSqrt(Sigma.Inv[[i]]) %x% diag(1,t-1)
#     
#     USA <- Data[1:p,]
#     
#     Response_USA <- matrix(c(t(USA[,-1])))
#     X_USA <- matrix(USA[,-t],t-1,p,byrow=TRUE)
#     X_Sigma[[i]] <- optim.product(Sigma.Sqrt.Inv1.new,X_USA,p)
#     Response_Sigma[[i]] <- Sigma.Sqrt.Inv1.new %*% Response_USA
#     
#   }
#   
#   
#   
#   
#   X.full <- block.diag(X_Sigma)
#   Y.full <- c(unlist(Response_Sigma))
#   
#   return(list(X=X.full,
#               Y=Y.full,
#               p=p,
#               t=t,
#               connList=connList))
# }


##########################################################
## function connListCalc: calculates a connection list ###
## for generalized lasso for p vars per entity  
## It matches corresponding elements of two transition matrices
## for the fusion penalty.
##########################################################

connListCalc <- function(p){  # returns the connection list
  
  connList <- vector("list",2*p^2)
  class(connList) <- "connListObj"
  for (i in 1:(p^2)) connList[[i]] <- as.integer(i+p^2-1)
  for (i in (p^2+1):(2*(p^2))) connList[[i]] <- as.integer(i- p^2-1)
  
  return(connList=connList)
}




#################################################################################
#################################################################################
### IN THAT SECTION I HAVE ALL THE POSSIBLE CRITERIONS                     ######
### FOR PICKING SPARSITY AND FUSION PARAMETERS                             ######
###                                                                        ######
### The parameters for all of them will be:                                ######
###                                                                        ######
###  Est - set of estimates that we calculate the criterion for,           ######
###  X - data matrix,                                                      ######
###  Y - response vector,                                                  ######
###  lambda.path - vector of sparsity parameter values                     ######
###                that estimates from 'Est' correspond to,                ######
###  df.coef - the coefficient for degrees of freedom to be multiplied by  ######
###                                                                        ######
### The functions will return:                                             ######  
###                                                                        ######
###  the estimate that minimizes the criterion,                            ######
###  the corresponding minimum value of the criterion,                     ######
###  set of criterion values for all possible sparsity parameter values,   ######
###  sparsity parameter value corresponding to minimum criterion value,    ######
###  index of that sparsity parameter value in lambda.path,                ######
###  log-likelihood part of the criterion(for the whole lambda path)       ######
###  degrees of freedom of the criterion(for the whole lambda path)        ######
#################################################################################
#################################################################################



################################################
### AIC criterion ##############################
###############################################

AIC <- function(Est,X,Y,lambda.path,df.coef=2){
  
  AIC <- rep(0,length(lambda.path))
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  n <- nrow(Y)
  t <- n/p + 1
  
  for(i in 1:length(lambda.path)){
    df <- sum(Est[,i] != 0)
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- df.coef*df
    AIC[i] <- loglik.part[i] + df.part[i]  
  }
  
  min <- which.min(AIC)
  
  return(list(Est=Est[,min],
              Criter.min=AIC[min],
              Criter=AIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}



################################################
### AICc criterion ##############################
###############################################

AICc <- function(Est,X,Y,lambda.path){
  
  AIC <- rep(0,length(lambda.path))
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  n <- nrow(Y)
  t <- n/p + 1
  
  for(i in 1:length(lambda.path)){
    df <- sum(Est[,i] != 0)
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- 2*df + 2*df*(df+1)/(n-df-1)
    AIC[i] <- loglik.part[i] + df.part[i]  
  }
  
  min <- which.min(AIC)
  
  return(list(Est=Est[,min],
              Criter.min=AIC[min],
              Criter=AIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}


#############################################
#### AIC with DISTINCT DEGREES OF FREEDOM 
#### (we only count DISTINCT non-zero elements as degrees of freedom)
#############################################

AIC.dist <- function(Est,X,Y,lambda.path,df.coef=2){
  AIC <- rep(0,ncol(Est))
  n <- nrow(Y)
  t <- n/p + 1
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  
  for(i in 1:ncol(Est)){
    df <- sum(Est[1:p^2,i] != 0) + sum((Est[(p^2+1):(2*p^2),i] != 0) & (Est[(p^2+1):(2*p^2),i] != Est[1:p^2,i]))  
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- df.coef*df
    AIC[i] <- loglik.part[i] + df.part[i]
  }

  min <- which.min(AIC)
  
  return(list(Est=Est[,min],
              Criter.min=AIC[min],
              Criter=AIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}


#############################################
#### AIC CORRECTED with DISTINCT DEGREES OF FREEDOM 
#### (we only count DISTINCT non-zero elements as degrees of freedom)
#############################################

AICc.dist <- function(Est,X,Y,lambda.path){
  AIC <- rep(0,ncol(Est))
  n <- nrow(Y)
  t <- n/p + 1
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  
  for(i in 1:ncol(Est)){
    df <- sum(Est[1:p^2,i] != 0) + sum((Est[(p^2+1):(2*p^2),i] != 0) & (Est[(p^2+1):(2*p^2),i] != Est[1:p^2,i]))  
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- 2*df + 2*df*(df+1)/(n-df-1)
    AIC[i] <- loglik.part[i] + df.part[i]
  }
  
  min <- which.min(AIC)
  
  return(list(Est=Est[,min],
              Criter.min=AIC[min],
              Criter=AIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}



################################################
### BIC criterion ##############################
###############################################

BIC <- function(Est,X,Y,lambda.path){
  
  BIC <- rep(0,ncol(Est))
  n <- nrow(Y)
  t <- n/p + 1
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  
  for(i in 1:ncol(Est)){ 
    df <- sum(Est[,i] != 0)
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- df*log(n)
    BIC[i] <- loglik.part[i] + df*log(n) 
  }
 
  min <- which.min(BIC)
  
  return(list(Est=Est[,min],
              Criter.min=BIC[min],
              Criter=BIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}


################################################
### BIC criterion with DISTINCT DEGREES OF FREEDOM 
### (we only count DISTINCT non-zero elements as degrees of freedom)
#################################################


BIC.dist <- function(Est,X,Y,lambda.path){
  
  BIC <- rep(0,ncol(Est))
  n <- nrow(Y)
  t <- n/p + 1
  loglik.part <- rep(0,length(lambda.path))
  df.part <- rep(0,length(lambda.path))
  
  for(i in 1:ncol(Est)){
    df <- sum(Est[1:p^2,i] != 0) + sum((Est[(p^2+1):(2*p^2),i] != 0) & (Est[(p^2+1):(2*p^2),i] != Est[1:p^2,i]))  
    loglik.part[i] <- 2*n*log(norm(Y - X %*% Est[,i],type="F")/sqrt(n))
    df.part[i] <- df*log(n)
    BIC[i] <- loglik.part[i] + df.part[i]
  }

  min <- which.min(BIC)
  #print(min)
  
  return(list(Est=Est[,min],
              Criter.min=BIC[min],
              Criter=BIC,
              lambda1=lambda.path[min],
              ind=min,
              loglik.part=loglik.part,
              df.part=df.part))
}






###################################################################
#### Function calculating H-step forecast values and MSFE's #######
###################################################################



forecast.mine <- function( # returns list of forecasted values,
                           #         list of corresponding MSFE
                           #         (Mean Squared Forecast Error) values
  Data,              # Data matrix
  Test=Data[,t - H], # Test data points
  K,                 # number of entities
  Est,               # List of length K of transition matrix estimates for each entity
  H=1                # Number of forecasting steps
){
  
  p <- nrow(Data)/K
  t <- ncol(Data)
  l <- 1        
  
  forecast.list  <- array(0, c(l, K*p, H))
  
  ### when we have more than 1 estimate, we run a loop for all of them
  # if (l>1){
  #   for (i in 1:l){
  #     forecast <- cbind(Data, array(0, c(K*p, H)))
  #     A <- block.diag(Est)
  #     
  #     for(k in 1:H){
  #       forecast[,1+k] <- A %*% forecast[,(1+k)-1]
  #     }
  #     forecast.list[i,,] <- forecast[,1+seq(H)]
  #   }
  #   
  #   forecast.err <- array(0, c(l, K*p, H))
  #   MSFE <- array(0,c(l,H))
  #   
  #   ## USING SWEEP TO SUBTRACT A ROW WITH TEST SUBSET VALUES FROM ALL THE PREDICTION ROWS
  #   for (k in 1:H){
  #     forecast.err[,,k] <- sweep(forecast.list[,,k],MARGIN=2,Test[,k],"-")^2
  #     MSFE[,k] <- apply(forecast.err[,,k],1,mean)
  #   }
  # }
  
   ### when we have just 1 estimate, we run the forecast just for this estimate 
  
  if (l==1){ 
    forecast <- cbind(Data, array(0, c(K*p, H)))
    A <- block.diag(Est)

    for(k in 1:H){
      forecast[,1+k] <- A %*% forecast[,(1+k)-1]
    }

    forecast.list[1,,] <- forecast[,1+seq(H)]
    forecast.err <- array(0, c(l, K*p, H))
    MSFE <- array(0,c(l,H))
  
  for (k in 1:H){
    forecast.err[,,k] <- (forecast.list[,,k] - as.matrix(Test[,k]))^2
    MSFE[,k] <- mean(forecast.err[,,k])
    }
  }
  return(list(val=forecast.list,
              MSFE=MSFE))
}


##########################################
#### Calculating Matthews coefficient based on TP,FP,TN,FN rates
#### of the estimate
##########################################

Matthews.Coef <- function(TP,FP,TN,FN){
  return(ifelse(TP*TN - FP*FN == 0,0,(TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
}




#########################################################
## Product of matrix M1 with block-diagonal matrix M2  ##
## Matrix M2 consists of identical blocks "Block"      ##
## Optimal calculation procedure through calculating   ##
## products of smaller submatrices                     ##
#########################################################

optim.product <- function(#returns product of M1 by M2
                          M1,
                          Block,    # block of M2 matrix
                          nblocks){ # number of diagonal blocks in M2
  dim1 <- nrow(M1)
  Prod <- matrix(0,dim1,ncol(Block)*nblocks)
  
  p <- nblocks
  t <- dim1/nblocks
  
  for (i in 1:nblocks) 
    for (j in 1:nblocks)
    Prod[(i-1)*t + 1:t,(j-1)*p + 1:p] <- M1[(i-1)*t + 1:t,(j-1)*t + 1:t] %*% Block
  
  return(Prod)
}




###########################################################
#### Function that initializes the ADMM algorithm #########
###########################################################

ADMM.precalc <- function(# returns number of columns of data matrix X,
                         #         rho - Lipshitz constant for the algorithm,
                         #         Z1,Z2,W - matrices needed for completing a step of ADMM algorithm, 
                         Y, # response vector 
                         X, # data matrix
                         Lm,  # matrix for the Lm*beta - gamma = 0 constraint
                         rho
                         ){                             
  
  ncol_X <- ncol(X)
  
  A <- Y
  B <- X
  D <- sqrt(rho/2)*Lm
  W <- fnMatSqrt(t(B)%*%B + t(D)%*%D)
  
  W.inv <- fnMatInverse(W)
  Z1 <- W.inv %*% t(B) %*% A
  Z2 <- W.inv %*% t(D)
  
  return(list(ncol_X=ncol_X,
              rho=rho,
              Z1=Z1,
              Z2=Z2,
              W=W))
}



################################# 
#### ADMM function ##############
####
#### Performs ADMM optimization for our generalized sparse fused lasso criterion. 
#### In particular, for a fixed value of lambda2 (fused lasso parameter) and a grid of values for lambda1 (sparsity parameter)
#### Update rules are as discussed in the paper
#### If grid for sparsity parameter is not submitted(lambda1.path), calculates it automatically.
#################################

fused.ADMM <- function(# returns the fused estimates for our grid,
                       #         set of u values(from ADMM algorithm),       
                       #         set of objective function values for the grid,
                       #         summary for number of iterations,
                       #         lambda1.path, calculated(if needed) for the sparsity parameter
                       Y,                  # Y - response vector          
                       X,                  # X - data matrix
                       Z1,Z2,W,rho,ncol_X, # output of ADMM.precalc function
                       Lm,                  # Lm - matrix for the ADMM constraint
                       lambda1.path=NULL,        # grid of sparsity parameter values
                       l_lambda1.path=NULL,       # expected length of grid for lambda1
                       lambda2,            # fusion parameter value
                       A.init,         # initial value for beta
                       gamma.start,        # initial value for gamma
                       u.start,            # initial value for u 
                       eps,                # stopping criterion
                       beta.true,          # true transition matrix value
                       iter.thresh,        # max number of iterations of the algorithm
                       fus.thresh          # hard threshold that shrinks differences between coefficients to 0
                      ){
  
  beta.start <- A.init
  gamma.start <- Lm %*% as.matrix(beta.start)
  beta.prev <- beta.start
  gamma.prev <- gamma.start
  u.prev <- u.start
 
  obj.prev.1 <- sum((Y - X %*% beta.start)^2) + lambda2*sum(abs(gamma.start))
  

    #################################
    #### INITIALIZING THE PATH ######
    #################################
  
  if (is.null(lambda1.path) == TRUE){
    C <- sqrt(rho/2)*(gamma.prev - u.prev)
    Z <- Z1 + Z2 %*% C
    glm.obj <- glmnet(x=X,y=Y,family="gaussian",standardize=standardize,intercept=intercept)
    glm.coef <- coef(glm.obj)
    lambda.path <- glm.obj$lambda
    lambda1.path <- lambda.path
    #print(paste(c("RANGE OF INITIAL PATH OF lambda1's:",head(lambda1.path,1),tail(lambda1.path,1),length(lambda1.path))))
  }
  
  l_lambda1.path <- length(lambda1.path)
  beta.est <- matrix(0,ncol_X,l_lambda1.path)
  obj.est <- matrix(0,1,l_lambda1.path)
  u.est <- matrix(0,ncol_X/2,l_lambda1.path)
  iter <- rep(0,l_lambda1.path)
  sum.shrunk <- rep(0,l_lambda1.path)

  
  ###########
  ### RUNNING THE ADMM FOR ALL THE lambda1 VALUES
  ###########
  
  ind <- 0
  
  for (lambda1 in lambda1.path){
    ind <- ind + 1
    beta.prev <- beta.start
    gamma.prev <- gamma.start
    u.prev <- u.start
    obj.prev <- lambda1*sum(abs(beta.start)) + obj.prev.1
    
    obj.val <- list()
    obj.val.diff <- list()
    
    iter[ind] <- 0
    
    obj.val[[iter[ind]+1]] <- obj.prev
    obj.val.diff[[iter[ind]+1]] <- obj.prev
    
    repeat{
      iter[ind] <- iter[ind] + 1
      C <- sqrt(rho/2)*(gamma.prev - u.prev)
      Z <- Z1 + Z2 %*% C
 
      ## beta update: solution to l1-optimization problem
      ## ||Z - W*beta||_2^2 + lambda1*||beta||_1

      glm.obj <- glmnet(W,Z,lambda=lambda1,family="gaussian",standardize=standardize,intercept=intercept)
      beta.next <- as.matrix(coef(glm.obj)[-1])
      Lm.b <- Lm %*% beta.next
   
      ## gamma update: proximal operator for lasso signal approximator
      ## (rho/2)*||Y - gamma||_2^2 + lambda2*||gamma|_1
      
      gamma.next <- softthresh(Lm.b+u.prev,(lambda2/rho))
      
      ## u update:
      
      u.next <- u.prev + Lm.b - gamma.next
      
      ## objective function update:
      
      obj.next <- sum((Y - X %*% beta.next)^2) + 
        lambda1*sum(abs(beta.next)) +
        lambda2*sum(abs(gamma.next))
      
      obj.val[[iter[ind]+1]] <- obj.next
      obj.val.diff[[iter[ind]+1]] <- (abs(obj.next - obj.prev)/abs(obj.prev))
         
      # hops <- ConvertToMatrix.Full(beta.next,p)
      # par(mfrow=c(1,2))
      # color2D.matplot(hops[[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
      #                 show.legend=TRUE,
      #                 show.values=TRUE,
      #                 main = paste(c("First entity(SHRINKING ITER)",lambda1,lambda2,iter[ind])))
      # color2D.matplot(hops[[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
      #                 show.legend=TRUE,
      #                 show.values=TRUE,
      #                 main = paste(c("Second entity(SHRINKING ITER)",lambda1,lambda2,iter[ind])))
      # 
      
      if (((abs(obj.next - obj.prev)/abs(obj.prev))<eps) | (iter[ind]>iter.thresh)) break;
      
      beta.prev <- beta.next
      gamma.prev <- gamma.next
      u.prev <- u.next  
      obj.prev <- obj.next      
    }
    
    
    ################################################
    ### shrinking small differences to zero   ######
    ################################################

    ### finding the indices of CORRESPONDING elements(of two matrices) that are within fus.thresh of each other
    ### taking average of two matrices(to use for finalizing fusion), BUT ALSO SETTING TINY AVERAGES(<fus.thresh) TO ZERO
    ### (in case we had one element at 0, corresponding element at 0.001 => 
    ### instead of setting them both to (0+0.001)/2 = 0.0005, we just shrink that to 0, ENCOURAGING SPARSITY)
    
    fus.shrink <- ifelse(abs(beta.next[1:(ncol_X/2)] - beta.next[((ncol_X/2)+1):ncol_X])<fus.thresh,1,0)
    beta.next.avg <- sparsify((head(beta.next,ncol_X/2) + tail(beta.next,ncol_X/2))/2,fus.thresh)
    
    beta.next[1:(ncol_X/2)] <- fus.shrink*beta.next.avg + (1-fus.shrink)*head(beta.next,(ncol_X/2))
    beta.next[((ncol_X/2)+1):ncol_X] <- fus.shrink*beta.next.avg + (1-fus.shrink)*tail(beta.next,(ncol_X/2))
    
    sum.shrunk[ind] <- sum(head(beta.next,(ncol_X/2)) == tail(beta.next,(ncol_X/2)))

    beta.est[,ind] <- beta.next
    u.est[,ind] <- u.next
    obj.est[,ind] <- obj.next

  }

  #### Here I have commented the plots of matrices for the case of a single submitted value of sparsity parameter,
  #### helped me see the effects of increasing the fusion parameter(lambda2), while having lambda1 fixed.
  
  if (l_lambda1.path == 1){ 
    if (lambda2>0){
      #print("Totals of iterations:")
      #print(iter)
      #print("Total of shrunk differences:")
      #print(sum.shrunk)
    }
      if (lambda2==0){
        #print("Totals of iterations:")
        #print(iter)
      }
    
      # hops <- ConvertToMatrix.Full(beta.est[,1],p)
      # par(mfrow=c(1,2))
      # color2D.matplot(hops[[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
      #                 show.legend=FALSE,
      #                 show.values=TRUE,
      #                 main="A(AFTER fusion)",
      #                 xlab="",
      #                 ylab="",
      #                 axes=FALSE)
      #                 #show.legend=TRUE,
      #                 #main = paste(c("First entity(lambda2 run)",lambda1,lambda2)))
      # color2D.matplot(hops[[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
      #                 show.legend=FALSE,
      #                 show.values=TRUE,
      #                 main="B(AFTER fusion)",
      #                 xlab="",
      #                 ylab="",
      #                 axes=FALSE)
      #                 #show.legend=TRUE,
      #                 #main = paste(c("Second entity(lambda2 run)",lambda1,lambda2)))


  }
  
  ### here I just had plots and printouts for the lambda1.grid with fixed lambda2.

  if (l_lambda1.path > 1){

    if (lambda2>0){
      #print("Totals of iterations:")
      #print(iter)
      #print("Total of shrunk differences:")
      #print(sum.shrunk)
    }
    
    if (lambda2==0){
      #print("Totals of iterations:")
      #print(iter)
    }
  #print(paste(c("Totals of iterations:",iter)))
 
    
  ##############################################
  ### PRINTOUT OF LAMBDA1 RUNS FOR JOINT  ######
  ##############################################
    
  #   for(l in 1:l_lambda1.path){
  # hops <- ConvertToMatrix.Full(beta.est[,l],p)
  # 
  # par(mfrow=c(1,2))
  # 
  # color2D.matplot(hops[[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
  #                 show.legend=TRUE,
  #                 show.values=TRUE,
  #                 main = paste(c("First entity(lambda1 run)",lambda1.path[l],lambda2)))
  # color2D.matplot(hops[[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
  #                 show.legend=TRUE,
  #                 show.values=TRUE,
  #                 main = paste(c("Second entity(lambda1 run)",lambda1.path[l],lambda2)))
  # }
    
  } 
  return(list(beta.est=as.matrix(beta.est),
              u.est=u.est,
              iter=iter,
              obj.est = obj.est,
              lambda.path = lambda1.path))
}



###################################################
### Function performing the grid search:
###  fix lambda2 at lambda2.init, optimize with respect to lambda1, get lambda1.est
###  fix lambda1 at lambda1.est, optimize with respect to lambda2
###################################################

seq.step <- function(# returns lambda1 and lambda2 picked by the criterion,
                     #         their indices in lambda1.path and lambda2.path respectively,
                     #         joint estimate corresponding to these parameter values
  Y,                  # response vector
  X,                  # data matrix
  Z1,Z2,W,rho,ncol_X, # ADMM precalculations output
  p,                  # number of variables per entity
  Lm,                  # matrix for the Lm*beta - gamma=0 constraint
  lambda1.path=NULL,        #path of lambda1 values
  l_lambda1.path=NULL,        # length of path of lambda1 values
  lambda2.path,       # path of lambda2 values
  lambda2.init,       # initial value for lambda2
  eps,                # stopping criterion
  A.true.vec,         # true value of transition matrices(for simulated data)
  iter.thresh,        # max number of iterations for ADMM
  fused.thresh,       # hard threshold that shrinks differences between coefficients to 0
  df=3,               # degrees of freedom multiplier for the criterion
  criter.joint.1="AICc.dist",  # the criterion
  criter.joint.2="BIC.dist",
  A.init=c(matrix(diag(1,p),byrow=TRUE),matrix(diag(1,p),byrow=TRUE)) #initializes vector beta
){
 
  ####Initialization part

  gamma.start <- Lm %*% as.matrix(A.init)
  u.start <- rep(0,p^2)
  u.Est <- list()

  #### Running the ADMM algorithm
  
  Joint.Result <- fused.ADMM(Y=Y,
                             X=X,
                             Z1=Z1,
                             Z2=Z2,
                             W=W,
                             rho=rho,
                             ncol_X=ncol_X,
                             Lm=Lm,
                             lambda1.path = lambda1.path,
                             l_lambda1.path=l_lambda1.path,
                             lambda2=lambda2.init,
                             gamma.start=gamma.start,
                             A.init=A.init,
                             u.start=u.start,
                             eps=eps,
                             beta.true=A.true.vec,
                             iter.thresh=iter.thresh, 
                             fus.thresh=fused.thresh)
  
  Joint.Result.Est <- Joint.Result$beta.est
  Joint.Result.lambda.path <- Joint.Result$lambda.path
  
  #### Looking for optimal estimate depending on the criterion
  #### Getting the optimal sparsity parameter value lambda1.est
  
  if (criter.joint.1 == "AIC.dist"){
    Criter.out <- AIC.dist(Joint.Result.Est,X,as.matrix(Y),lambda.path=Joint.Result.lambda.path,df.coef=df)  
  }
  if (criter.joint.1 == "AICc.dist"){
    Criter.out <- AICc.dist(Joint.Result.Est, X,as.matrix(Y),lambda.path=Joint.Result.lambda.path)  
  }
  if (criter.joint.1 == "AIC"){
    Criter.out <- AIC(Joint.Result.Est,X,as.matrix(Y),lambda.path=Joint.Result.lambda.path)  
  }
  if (criter.joint.1 == "AICc"){
    Criter.out <- AICc(Joint.Result.Est,X,as.matrix(Y),lambda.path=Joint.Result.lambda.path)  
  }
  if (criter.joint.1 == "BIC.dist"){
    Criter.out <- BIC.dist(Joint.Result.Est,X,as.matrix(Y),lambda.path=Joint.Result.lambda.path)  
  }
  if (criter.joint.1 == "BIC"){
    Criter.out <- BIC(Joint.Result.Est,X,as.matrix(Y),lambda.path=Joint.Resul.lambda.path)  
  }
  
  
  Criter.J <- Criter.out$Criter.min
  Criter.val <- Criter.out$Criter
  J.Est <- Criter.out$Est 
  lambda1.est <- Criter.out$lambda1
  ind1.est <- Criter.out$ind
  
  
  # hops <- ConvertToMatrix.Full(Joint.USACANADA.Est[,ind1.est],p)
  # par(mfrow=c(1,2))
  # color2D.matplot(hops[[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
  #                 #show.legend=FALSE,
  #                 show.values=TRUE,
  #                 xlab="",
  #                 ylab="",
  #                 axes=FALSE,
  #                 main = paste(c("First entity, sparse selected",lambda1.est,lambda2.init)),
  #                 show.legend=TRUE)
  # color2D.matplot(hops[[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
  #                 #show.legend=FALSE,
  #                 show.values=TRUE,
  #                 xlab="",
  #                 ylab="",
  #                 axes=FALSE,
  #                 main = paste(c("Second entity, sparse selected",lambda1.est,lambda2.init)),
  #                 show.legend=TRUE)
  # 
  
  
  #####  Fixing lambda1=lambda1.est and running ADMM algorithm for 
  #####  fusion parameter lambda2
  
  Joint.Result.Est <- list()
  Criter.J <- rep(0,L2)
  Criter.val <- rep(0,L2)
  J.Est <- list()
  
  it <- 0
  
  for (lambda2 in lambda2.path){    
    it <- it+1
    #print(c("lambda2:",lambda2))
    
    Joint.Result <- fused.ADMM(Y=Y,
                               X=X,
                               Z1=Z1,
                               Z2=Z2,
                               W=W,
                               rho=rho,
                               ncol_X=ncol_X,
                               Lm=Lm,
                               lambda1.path=lambda1.est,
                               l_lambda1.path=NULL,
                               lambda2=lambda2,
                               gamma.start=gamma.start,
                               A.init=A.init,
                               u.start=u.start,
                               eps=eps,
                               beta.true=A.true.vec,
                               iter.thresh=iter.thresh,
                               fus.thresh=fused.thresh)
    
    Joint.Result.Est[[it]] <- Joint.Result$beta.est
    
    #### Picking estimate and value of lambda2 that minimizes corresponding criterion
    
    if (criter.joint.2 == "AIC.dist"){
      Criter.out <- AIC.dist(Joint.Result.Est[[it]],X,as.matrix(Y),lambda.path=lambda1.est,df.coef=df)  
    }
    if (criter.joint.2 == "AICc.dist"){
      Criter.out <- AICc.dist(Joint.Result.Est[[it]],X,as.matrix(Y),lambda.path=lambda1.est) 
    }
    if (criter.joint.2 == "AIC"){
      Criter.out <- AIC(Joint.Result.Est[[it]],X,as.matrix(Y),lambda.path=lambda1.est)  
    }
    if (criter.joint.2 == "BIC.dist"){
      Criter.out <- BIC.dist(Joint.Result.Est[[it]],X,as.matrix(Y),lambda.path=lambda1.est) 
    }
    if (criter.joint.2 == "BIC"){
      Criter.out <- BIC(Joint.Result.Est[[it]],X,as.matrix(Y),lambda.path=lambda1.est)  
    }
 
    Criter.J[it] <- Criter.out$Criter.min
    Criter.val[it] <- Criter.out$Criter
    J.Est[[it]] <- Criter.out$Est 
  } 
  
  ind2.est <- which.min(Criter.val)
  lambda2.est <- lambda2.path[ind2.est]
  Joint.Est <- J.Est[[ind2.est]]
  
  return(list(lambda1.est=lambda1.est,
              ind1.est=ind1.est,
              lambda2.est=lambda2.est,
              ind2.est=ind2.est,
              Joint.Est=Joint.Est,
              lambda1.path=Joint.Result.lambda.path))
}




##########################################################
#### GENERATING SYNTHETIC DATA FOR TWO ENTITIES   ########
##########################################################

Gener.Simul.Data <- function(# USES A LOT OF GLOBAL VARIABLES
                             # returns both generated transition matrices,
                             #         both generated Sigma matrices and their inverses,
                             #         generated training and test datasets
  ){

  #####################################################
  ### Generating Sigma matrices from L-factor model ###
  #####################################################
  
  Sigma.true <- list()
  Sigma.U <- list()
  Sigma.Inv.true <- list()
  
  for(j in 1:K){
   Gen <- Sigma.Gen(p,t,L)
   Sigma.true[[j]] <- Gen$Sigma
   Sigma.U[[j]] <- Gen$Sigma_Ux
   Sigma.Inv.true[[j]] <- solve(Sigma.true[[j]])
  }
  
   Sigma <- block.diag(Sigma.true)
   
  ####################################################
  ### Generating stationary transition matrices for K related VAR(1) models  
  ####################################################
  
  A.object <- A.setup(p=p,thresh=max_eig-0.2,diff=diff,ed=ed,comm=comm,max_eig=max_eig)
  
  A.true.mat <- list()
  A.true.mat[[1]] <- A.object[[1]]
  
  ## getting vectorized versions of matrices, by row
  A.true.vec <- matrix(t(A.object[[1]]),1,p^2)
  
  for(j in 2:K){
  A.true.mat[[j]] <- A.object[[j]]
  A.true.vec <- cbind(A.true.vec,matrix(t(A.object[[j]]),1,p^2))
  }
  
  
  ## Full block-diagonal matrix A with blocks A1, A2, ..., AK
  
  A.full <- block.diag(A.true.mat)
  A.list <- list()
  A.list[[1]] <- A.full
  
  
  ##########
  ####### DATA GENERATION 
  ##########

  GenData <- gen_dat(T=t,Sigma_error=Sigma,A=A.list,SNR=SNR)
  GenData.train <- GenData[,1:train]
  GenData.test <- GenData[,train+1:h]
  
  return(list(A.true=A.true.mat,
              Sigma.true=Sigma.true,
              Sigma.Inv.true=Sigma.Inv.true,
              GenData.train=GenData.train,
              GenData.test=GenData.test))
  
}





################################################
#### SEPARATE ESTIMATION FOR SIMULATED DATA ####
################################################


Simul.Data.Separate.Main <- function(# USES A LOT OF GLOBAL VARIABLES
                                     # returns the picked tuning parameter values(lambda1),
                                     #         the transition matrix estimates,
                                     #         all performance measures of these estimates
                                     #         (FP,FN,TP,TN,Matthews,Frob diff,MSFE),
                                     #         error covariance and its inverse estimates,
                                     #         performance measures of these estimates
                                     #         (Frobenius differences)
  ){
  
  ##################################################################
  ### Calculating l1-estimates of transition matrices               ##
  ### with identity error covariance matrices                     ##
  ##################################################################
  
  Entity.Sep.Est <- list()
  lambda1.picked <- numeric(K)
  ind1.picked <- numeric(K)
  
  for(k in 1:K){
    Entity <- GenData.train[(k-1)*p + 1:p,]
    Entity.Sep <- sep.lasso.Inv(Entity,diag(1,p))
    
    ## in case generated data leads to overly sparse glmnet solution path,
    ## return 0 and regenerate the data
    
    if (fail == 1) return(0);
    
    ## tuning parameter selection for separate method
    if (criter.sep == "AIC") criter.out <- AIC(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path,df.coef=df.sep)
    if (criter.sep == "AICc") criter.out <- AICc(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path)
    
    Entity.Sep.Est[[k]] <- criter.out$Est
    lambda1.picked[k]  <- criter.out$lambda1
    ind1.picked[k]  <- criter.out$ind
    
    if (ConstThresh == TRUE){
      Entity.Sep.Est[[k]] <- sparsify(criter.out$Est,Thresh)
    }
  }

  ##Using residuals off these initial l1-estimates as data for error covariance estimation
  ##Apply L-factor model to estimate off-diagonal elements of inverse error covariance
  ##Apply straightforward graphical lasso to estimate diagonals of inverse error covariance
  
  L1.est <- numeric(K)
  TVE1.est <- numeric(K)
  Sigma.Inv.est <- list()
  Frob <- numeric(K)
  Frob.Inv <- numeric(K)
  
  for (j in 1:sigma.iter){
   for(k in 1:K){
   Entity <- GenData.train[(k-1)*p + 1:p,]
   Entity_1step <- L1.1step(Entity, A.L1.est=Entity.Sep.Est[[k]]) 
   L1.est[k]  <- Entity_1step$L
   TVE1.est[k] <- Entity_1step$TVE
   
   Sigma.Inv.est[[k]] <- Entity_1step$Sigma.est.Inv
   Frob[k]  <- Frob.comp(solve(Sigma.Inv.est[[k]]),Sigma.true[[k]])
   Frob.Inv[k]  <- Frob.comp(Sigma.Inv.est[[k]],Sigma.Inv.true[[k]])
   
   Entity.Sep <- sep.lasso.Inv(Entity,Sigma.Inv.est[[k]])
   
   if (fail == 1) return(0);
   
   if (criter.sep == "AIC") criter.out <- AIC(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path,df.coef=df.sep)
   if (criter.sep == "AICc") criter.out <- AICc(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path)
   
   Entity.Sep.Est[[k]] <- criter.out$Est
   lambda1.picked[k]  <- criter.out$lambda1
   ind1.picked[k]  <- criter.out$ind
   
   if (ConstThresh == TRUE){
     Entity.Sep.Est[[k]] <- sparsify(criter.out$Est,Thresh)
   }
   }
  }
    
    Measures.Sep <- Measures.Vec(c(unlist(Entity.Sep.Est)),A.true.vec)
    FP.Sep  <- Measures.Sep$FP
    FN.Sep  <- Measures.Sep$FN
    TP.Sep  <- Measures.Sep$TP
    TN.Sep  <- Measures.Sep$TN
    Matt.Coef.Sep  <- Matthews.Coef(TP.Sep,
                                    FP.Sep,
                                    TN.Sep,
                                    FN.Sep)
    Frob.Sep  <- Measures.Sep$Frob
    
    Sep.Est  <- ConvertToMatrix.Full(c(unlist(Entity.Sep.Est)),p)
    
    f.mine <- forecast.mine(as.matrix(GenData.train[,train]), Test = as.matrix(GenData.test),Est = Sep.Est, H=h,K=2)
    Pred.Err.Sep  <- f.mine$MSFE[1]
    
    Results <- c(Pred.Err.Sep,
                 FP.Sep,
                 FN.Sep,
                 Matt.Coef.Sep,
                 Frob.Sep)
    
    return(list(lambda1.picked=lambda1.picked,
                ind1.picked=ind1.picked,
                L1.est=L1.est,
                TVE1.est=TVE1.est,
                Frob=Frob,
                Frob.Inv=Frob.Inv,
                Pred.Err.Sep=Pred.Err.Sep,
                FP.Sep=FP.Sep,
                FN.Sep=FN.Sep,
                TP.Sep=TP.Sep,
                TN.Sep=TN.Sep,
                Matt.Coef.Sep=Matt.Coef.Sep,
                Frob.Sep=Frob.Sep,
                Sep.Est=Sep.Est,
                Sigma.Inv.est=Sigma.Inv.est,
                Results=Results))
    

}




################################################
#### JOINT ESTIMATION FOR SIMULATED DATA    ####
################################################



Simul.Data.Joint.Main <- function( # USES A LOT OF GLOBAL VARIABLES
                                   # returns the picked tuning parameter values(lambda1 and lambda2),
                                   #         the transition matrix estimates,
                                   #         all performance measures of these estimates
                                   #         (FP,FN,TP,TN,Matthews,Frob diff,MSFE),
                                   #         error covariance and its inverse estimates,
                                   #         performance measures of these estimates
                                   #         (Frobenius differences) 
  ){
  
  #################################
  #### ADMM function ##############
  #################################
 
  #### Running joint method
  #### as multiple runs of the seq.step function
  
  Data.P <- Data.Prep.Inv(GenData.train,Sigma.Inv.est)
  
  ## setting up the L matrix for the gamma - Lbeta = 0 constraint in ADMM task
  
  Lm <- NULL
  for (j in 1:(K-1)){
    for (l in (j+1):K){
      Lvec <- numeric(K)
      Lvec[j] <- 1
      Lvec[l] <- -1
      Lsub <- diag(Lvec[1],p^2)
      for (o in 2:K) Lsub <- cbind(Lsub,diag(Lvec[o],p^2))
      if (!is.null(Lm)) Lm <- rbind(Lm,Lsub);
      if (is.null(Lm)) Lm <- Lsub;
    }
  }  
  
  precalc <- ADMM.precalc(Data.P$Y,Data.P$X,Lm,rho=rho)
  
  ########
  ## Running the grid search
  ########
  
  lambda2.init <- 0
  for (i in 1:n.iter){ 
    Result <- seq.step(Y=Data.P$Y,
                       X=Data.P$X,
                       Z1=precalc$Z1,
                       Z2=precalc$Z2,
                       W=precalc$W,
                       rho=precalc$rho,
                       ncol_X=precalc$ncol_X,
                       p=p,
                       Lm=Lm,
                       lambda1.path=NULL,
                       l_lambda1.path=L1,
                       lambda2.path=lambda2.path,
                       lambda2.init=lambda2.init,
                       eps=eps,
                       A.true.vec=A.true.vec,
                       iter.thresh=iter.thresh,
                       fused.thresh=fused.thresh,
                       df=df,
                       criter.joint.1=criter.joint.1,
                       criter.joint.2=criter.joint.2,
                       A.init=A.init)
    
    lambda2.init <- Result$lambda2.est
    lambda1.path_est <- Result$lambda1.path
  }
  
  if (ConstThresh == TRUE){
    Result$Joint.Est <- sparsify(Result$Joint.Est,Thresh)
  }
  
  Joint.Est <- ConvertToMatrix.Full(Result$Joint.Est,p)
  
  f.mine <- forecast.mine(as.matrix(GenData.train[,train]), Test = as.matrix(GenData.test),Est = Joint.Est, H=h,K=2)
  Pred.Err.Joint  <- f.mine$MSFE[1]
  
  Measures.Joint <- Measures.Vec(Result$Joint.Est,A.true.vec)
  FP.Joint  <- Measures.Joint$FP
  FN.Joint  <- Measures.Joint$FN
  TP.Joint  <- Measures.Joint$TP
  TN.Joint  <- Measures.Joint$TN
  Matt.Coef.Joint  <- Matthews.Coef(TP.Joint,
                                    FP.Joint,
                                    TN.Joint,
                                    FN.Joint)
  Frob.Joint  <- Measures.Joint$Frob
  
  
  Joint.lambda1.picked  <- Result$lambda1.est
  Joint.lambda2.picked  <- Result$lambda2.est
  Joint.ind1.picked  <- Result$ind1.est
  Joint.ind2.picked  <- Result$ind2.est 
  
  Results <- c(Pred.Err.Joint,
               FP.Joint,
               FN.Joint,
               Matt.Coef.Joint,
               Frob.Joint)
  
  return(list(Pred.Err.Joint=Pred.Err.Joint,
              FP.Joint=FP.Joint,
              FN.Joint=FN.Joint,
              TP.Joint=TP.Joint,
              TN.Joint=TN.Joint,
              Matt.Coef.Joint=Matt.Coef.Joint,
              Frob.Joint=Frob.Joint,
              Joint.Est=Joint.Est,
              Joint.lambda1.picked=Joint.lambda1.picked,
              Joint.lambda2.picked=Joint.lambda2.picked,
              Joint.ind1.picked=Joint.ind1.picked,
              Joint.ind2.picked=Joint.ind2.picked,
              Results=Results, 
              lambda.path=lambda1.path_est))
}



##########################################
#### LOADING THE REAL DATA             ###
##########################################


Load.Data <- function(#returns the time series matrix(K*p X T),
                      #        number of variables
                      sit,        # situation number(1,2,3,4)
                      state,      # states under consideration
                      K           # number of entities
                      ){ 
  Data <- list()
  
  for (i in 1:K){
    Data[[i]] <- read.csv(paste(state[i],".csv",sep=""),header=TRUE)
    Data[[i]] <- na.omit(Data[[i]][,-1])
  }
  
  Data.final <- list()
  
  if (sit==1){
    p <- 8
    start <- 35
    period <- 70
    good.ind <- c(1,2,3,5,9,6,7,8)
    for(i in 1:K){
      Data.final[[i]] <- Data[[i]][start:(start+period),good.ind]
    }
  }
  
  if (sit==2){
    p <- 13
    start <- 35
    period <- 70
    good.ind <- c(1,2,3,5,9,12,18,16,14,10,6,7,8)
    for(i in 1:K){
      Data.final[[i]] <- Data[[i]][start:(start+period),good.ind]
    }
  }
  
  if (sit==3){
    p <- 13
    start <- 35
    period <- 70
    good.ind <- c(1,2,3,5,9,13,19,17,15,11,6,7,8)
    for(i in 1:K){
      Data.final[[i]] <- Data[[i]][start:(start+period),good.ind]
    }
  }
  
  if (sit==4){
    p <- 18
    start <- 35
    period <- 70
    good.ind <- c(1,2,3,5,9,12,18,16,14,10,13,19,17,15,11,6,7,8) 
    for(i in 1:K){
      Data.final[[i]] <- Data[[i]][start:(start+period),good.ind]
    }
  }
  
  Data.final.m <- t(Data.final[[1]])
  
  for(i in 2:K){
    Data.final.m <- rbind(Data.final.m,t(Data.final[[i]]))
  }
  
  return(list(Data.final.m=Data.final.m,
              p=p))
}



#########################################
### SEPARATE ESTIMATION FOR REAL DATA         
#########################################

## Each replicate inside the function corresponds to a "shifting subsets procedure":
##  first replicate corresponds to training subset from time point 'start' to point 'start+train'
##                                 testing time point is 'start+train+1'
##  second replicate corresponds to subset shifted by 1 to the right:
##    it goes from time point 'start+1' to 'start+train+1',
##    testing time point is 'start+train+2'
##  etc 

Real.Data.Separate.Main <- function(# USES A LOT OF GLOBAL VARIABLES
                                    # returns separate estimates,
                                    #         prediction errors,
                                    #         stability of estimates
                                    #         (how often are particular elements non-zero),
                                    #         inverse of Sigma estimate for each entity,
                                    #         estimate of number of factors L per entity
                                    ){

t0 <- start0-2

for (run in 1:rep){
  print(c("Run:",run))
  
  t0 <- t0 + 1
  
  ## standardizing the whole subset of size train+1(train - size of training set, 1 - size of testing set, as we make 1-step forecasts) 
  ## starting at point t0+1
  
  Data.st <- t(scale(t(Data.final.m[,(t0+1):(t0+train+1)]),center=TRUE))
  
  ## breaking into training & testing subsets for each entity
  Data.train <- list()
  Data.test <- list()
  
  for (k in 1:K){
    Data.train[[k]] <- Data.st[1:p + (k-1)*p,1:train]
    Data.test[[k]] <- Data.st[1:p + (k-1)*p,train+1:h]
  }
  
  Full.data.train <- Data.st[,1:train]
  Full.data.test <- Data.st[,train+1:h]
  
  ############
  ## Running the whole singly estimation process
  ## for each entity(US state, in our case)
  #############
  
  A.Sep.Est <- list()
  
  for(k in 1:K){
    Entity.Sep <- sep.lasso.Inv(Data.train[[k]],diag(1,p))

    if (criter.sep == "AIC") AIC.out <- AIC(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path,df.coef=df.sep)
    if (criter.sep == "AICc") AIC.out <- AICc(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path)
   
    A.Sep.Est[[k]] <- AIC.out$Est

    if (ConstThresh == TRUE){
      A.Sep.Est[[k]] <- sparsify(AIC.out$Est,Thresh)
    }

    for (j in 1:sigma.iter){
    Entity_1step <- L1.1step(Data.train[[k]],A.L1.est=A.Sep.Est[[k]])
    L1.est[k,run] <- Entity_1step$L
    TVE1.est[k,run] <- Entity_1step$TVE
    Sigma.Inv.est[[k]] <- Entity_1step$Sigma.est.Inv
    Entity.Sep <- sep.lasso.Inv(Data.train[[k]], Sigma.Inv=Sigma.Inv.est[[k]])

    if (criter.sep == "AIC") AIC.out <- AIC(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path,df.coef=df.sep)
    if (criter.sep == "AICc") AIC.out <- AICc(Entity.Sep$Est,Entity.Sep$X,as.matrix(Entity.Sep$Y),lambda.path=Entity.Sep$lambda.path)

    A.Sep.Est[[k]] <- AIC.out$Est

    if (ConstThresh == TRUE){
      A.Sep.Est[[k]] <- sparsify(AIC.out$Est,Thresh)
    }

  }
  }
  
  ## calculating performance measures: MSFE and Stability of estimates

  Sep.Est[[run]] <- ConvertToMatrix.Full(unlist(A.Sep.Est),p)
  
  f.mine <- forecast.mine(as.matrix(Full.data.train[,train]), 
                          Test = as.matrix(Full.data.test),
                          K=K,
                          Est = Sep.Est[[run]], 
                          H=h)
  
  Pred.Err.Sep[run] <- f.mine$MSFE[1]
  
  for (k in 1:K){
    Sep.Stability[[k]] <- Sep.Stability[[k]] + unlist(ConvertToMatrix.Full(ifelse(A.Sep.Est[[k]]==0,0,1),1))
  }
}

return(list(Sep.Est = Sep.Est,
            Pred.Err.Sep = Pred.Err.Sep,
            Sep.Stability = Sep.Stability,
            Sigma.Inv.est = Sigma.Inv.est,
            L1.est = L1.est,
            TVE1.est = TVE1.est))
}



####################################################
#### JOINT ESTIMATION(PAIRWISE APPROACH):
###############################################

##############################################
#### FOR THE CASE OF K=4:
#### Perform pairwise joint estimation of 
####   1st and 2nd entity(US state, in our case); 
####   1st and 3rd;
####   1st and 4th;
####   2nd and 3rd;
####   2nd and 4th;
####   3rd and 4th.
####
####  Afterwards combine all the estimates for each entity for further analysis.
####################################################

Real.Data.Joint.Main <- function(# returns joint estimates,
                                 #         prediction errors,
                                 #         stability of estimates
                                 #         (how often particular elements are non-zero)
                                 ){
  
  t0 <- start0-2
  
  for (run in 1:rep){
    print(c("Run:",run))
    t0 <- t0 + 1
    
    ### initializing all transition matrices with identities
    A.init <- list()
    for (k in 1:K){
      A.init[[k]] <- list()
    }
    
    ### Making pairwise snake-like joint estimation of K states
    ### Each pairwise estimation follows the sequential approach
    
    for (q in 1:(K-1)){
      for(l in (q+1):K){
      Two_Entities <- t(scale(t(Data.final.m[c((q-1)*p + 1:p, (l-1)*p + 1:p),(t0+1):(t0+train+1)]),center=TRUE))
      
      Two_Entities.train <- Two_Entities[,1:train]
      Two_Entities.test <- Two_Entities[,train+1:h]
      
      connList <- connListCalc(p)
      Data.P <- Data.Prep.Inv(Two_Entities.train,list(Sigma.Inv.est[[q]],Sigma.Inv.est[[l]]))
      
      ## setting up the L matrix for two entities optimization criterion, s.t L*(beta1,beta) = beta1 - beta2
      Lm <- cbind(diag(1,p^2),diag(-1,p^2))
      
      ## setting up the rest of the matrices for ADMM algorithm
      precalc <- ADMM.precalc(Data.P$Y,Data.P$X,Lm,rho=rho)
      lambda2.init <- 0
      
      for (i in 1:n.iter){
        Result <- seq.step(Y=Data.P$Y,
                           X=Data.P$X,
                           Z1=precalc$Z1, 
                           Z2=precalc$Z2, 
                           W=precalc$W, 
                           rho=precalc$rho, 
                           ncol_X=precalc$ncol_X,
                           p=p,
                           Lm=Lm,
                           lambda1.path=NULL,
                           l_lambda1.path=L1,
                           lambda2.path=lambda2.path,
                           lambda2.init=lambda2.init,
                           eps=eps,
                           A.true.vec=A.true.vec,
                           iter.thresh=iter.thresh,
                           fused.thresh=fused.thresh,
                           df=df,
                           criter.joint.1=criter.joint.1,
                           criter.joint.2=criter.joint.2,
                           A.init=c(matrix(diag(1,p),byrow=TRUE),matrix(diag(1,p),byrow=TRUE))
        )
        lambda2.init <- Result$lambda2.est
      }
      A.init[[q]] <- c(A.init[[q]],Result$Joint.Est[1:p^2])
      A.init[[l]] <- c(A.init[[l]],Result$Joint.Est[1:p^2 + p^2])
      
      }
    }
    
    Joint.Est[[run]] <- list()
      
    for(l in 1:(K-1)){
    A.init.med <- list()
    for(q in 1:K){
      A.init.med <- c(A.init.med,A.init[[q]][(l-1)*p^2 + 1:p^2])
    }
    Joint.Est[[run]][[l]] <- ConvertToMatrix.Full(sparsify(unlist(A.init.med),Thresh),p)
    }
    
    ### keeping track of stability measures
    for (k in 1:K){
      for(l in 1:(K-1)){
      Joint.Stability[[k]] <- Joint.Stability[[k]] + unlist(ConvertToMatrix.Full(ifelse(sparsify(unlist(A.init[[k]])[(l-1)*p^2 + 1:p^2],Thresh)==0,0,1),1))
      }
    }
    
    ## standardizing the whole subset of size "train" starting at point "t0+1"
    Data.st <- t(scale(t(Data.final.m[,(t0+1):(t0+train+1)]),center=TRUE))
    Full.data.train <- Data.st[,1:train]
    Full.data.test <- Data.st[,train+1:h]
    
    
    ### calculating the MSFE
    for(l in 1:(K-1)){
     f.mine <- forecast.mine(as.matrix(Full.data.train[,train]), 
                             Test = as.matrix(Full.data.test),
                             K=K,
                             Est = Joint.Est[[run]][[l]], 
                             H=h)
     Pred.Err.Joint[run] <- Pred.Err.Joint[run] + f.mine$MSFE[1]
    }
    Pred.Err.Joint[run] <- Pred.Err.Joint[run]/(K-1)
  }
    return(list(Joint.Est=Joint.Est,
                Joint.Stability=Joint.Stability,
                Pred.Err.Joint=Pred.Err.Joint))
}


############################
##### FUNCTION OUTPUTTING NAMES OF VARIABLES UNDER CONSIDERATION
##### IN THE CORRESPONDING SITUATION OF REAL DATA APPLICATION
###########################

situation.names <- function(sit){
  if (sit==1){
    return(c("Constr Total","Edu/Health Total","Finance Total","Manuf Total","GProd Total",
                             "Total Nonfarm","Lead Ind","UnempR"))
  }
  if (sit==2){ 
    return(c("Constr Total","Edu/Health Total","Finance Total","Manuf Total","GProd Total",
             "Constr WeeklyH","Edu/Health WeeklyH","Finance WeeklyH","Manuf WeeklyH","GProd WeeklyH",
              "Total Nonfarm","Lead Ind","UnempR"))
  }
  if (sit==3){
    return(c("Constr Total","Edu/Health Total","Finance Total","Manuf Total","GProd Total",
             "Consr HEarn","Edu/Health HEarn","Finance HEarn","Manuf HEarn","GProd HEarn",
             "Total Nonfarm","Lead Ind","UnempR"))
  }
  if (sit==4){
    return(c("Constr Total","Edu/Health Total","Finance Total","Manuf Total","GProd Total",
             "Constr WeeklyH","Edu/Health WeeklyH","Finance WeeklyH","Manuf WeeklyH","GProd WeeklyH",
             "Consr HEarn","Edu/Health HEarn","Finance HEarn","Manuf HEarn","GProd HEarn",
             "Total Nonfarm","Lead Ind","UnempR"))
  }
}


#### Plotting the heatmaps for stability matrices

simple.plot <- function(mat1,           # stability matrix
                        rown,          # names of variables under consideration
                        state,         # state under consideration
                        meth="Joint"   # method under consideration
                        ){
  par(mar=c(9,9,3,3))
  rownames(mat1) <- rown
  colnames(mat1) <- rown
  
  color2D.matplot(mat1,cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                  show.values=TRUE,
                  axes=FALSE,
                  xlab="",
                  ylab="",
                  show.legend=F,
                  main=paste(state,meth,sep=","))
  
  axis(1,at=0.5:(length(rown)-.5),las=2,labels=colnames(mat1))
  axis(2,at=0.5:(length(rown)-.5),las=2,labels=rev(rownames(mat1)))
}

