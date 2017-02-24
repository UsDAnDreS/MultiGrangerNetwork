#setwd("/home/usdandres/Documents/Study_stuff/George/Joint Granger Project/Research20/NEW STUFF, FEBRUARY 13TH")
# Copyright 2017, Andrey Skripnikov, All rights reserved.
source("Paper.functions.R")

K <- 4                             # number of entities(states)
state <- c("PA","IL","OH","MI")     # state abbreviations
sit <- 1                            # situation number(either 2,3 or 4)

########
########
#### Out of 70 available time points:
########
########

train <- 30                    # number of training points
start0 <- 1                     # the time point from which the training starts
rep <- 2                       # number of runs
sigma.iter <- 2                 # number of iterations to get Sigma estimate(before I always just used 1)
h <- 1                          # H-step forecasts will be made
test <- train + 1:h             # test time points(to measure forecasting precision)


#### criterions for tuning parameter selection
criter.joint.1 <- "AICc"  # criterion to select sparsity parameter for joint estimates
criter.joint.2 <- "BIC.dist" # criterion to select fusion parameter for joint estimates
criter.sep <- "AICc"       # criterion to select sparsity parameter for separate estimates
df.sep <- 2                     # the multiplier for "degrees of freedom" part of the criterion for SEPARATE


### parameters for glmnet() function performing solution path calculation for standard lasso problem
intercept <- FALSE
standardize <- TRUE
fail <- 0
fail.thresh <- 0



#####
### ADMM algorithm parameters
#####

rho <- 10                   # ADMM rho penalty parameter
init <- "0"                 # initialize ADMM with Identity("1"),Zero-matrix("0"), or with separate estimates("Sep")
eps <- 0.001           # stopping criterion
iter.thresh <- 200     # max number of iterations per one run of ADMM algorithm
fused.thresh <- 0.01  # threshold for differences between elements
n.iter <- 1            # number of sequential iterations

#####
### Thresholding parameters
#####

ConstThresh <- TRUE         # whether to use the hard threshold for final estimates
Thresh <- 0.1               # the value of hard threshold
df <- 2                     # for AIC.dist

######
## Printing out all the main parameter values
######

#print(c("p:",p))
print(c("start0:",start0))
print(c("rep",rep))
print(c("train:",train))
print(c("ConstThresh",ConstThresh))
print(c("Thresh",Thresh))
print(c("Init:",init))
print(c("eps:",eps))
print(c("iter.thresh",iter.thresh))
#print(c("lambda1:",lambda.path))
print(c("n.iter:",n.iter))
#print(c("lambda2:",lambda2.path))
print(c("df:",df))
print(c("situation:",sit))
print(c("df.sep",df.sep))
print(c("criter.joint.1:",criter.joint.1))
print(c("criter.joint.2:",criter.joint.2))



#####
### Starting the code timer
#####

ptm <- proc.time()


#####
## INITIALIZING VARIABLES
#####

Pred.Err.Sep <- rep(0,rep)
Pred.Err.Joint <- rep(0,rep)

Pred.Err.Sep.Refit <- rep(0,rep)
Pred.Err.Joint.Refit <- rep(0,rep)

Sep.Est <- list()
Sep.Est.Refit <- list()


L1.est <- matrix(0,K,rep)    # number of factors picked for all states and all replicates
TVE1.est <- matrix(0,K,rep)  # total variance explained by the picked factors for all states and all replicates

Joint.Est <- list()
Joint.Est.Refit <- list()
Sigma.Inv.est <- list()


#########################################
### LOADING THE DATA ####################
#########################################

### function Load.Data is described in Paper.functions.R
LData <- Load.Data(sit=sit,state=state,K=K)
Data.final.m <- LData$Data.final.m            ### full data matrix of size (K*p x T) with all K econometric time series
p <- LData$p                                  ### number of variables per state



#####
## Grid parameters
#####

knots1 <- 15                                           # controls lower bound of lambda1 path
knots2 <- 20                                           # controls lower bound of lambda2 path
L1 <- NULL                                              # length of lambda1 grid
L2 <- 20                                                # length of lambda2 grid
lambda2.path <- c(0,(2*p)*rho*(1/2)^seq(0, knots1, length = L2))   # lambda2 grid
L2 <- length(lambda2.path)

Sep.Stability <- ConvertToMatrix.Full(numeric(K*p^2),p)
Joint.Stability <- ConvertToMatrix.Full(numeric(K*p^2),p)


#########################################
### SEPARATE ESTIMATION              ####
#########################################

t0 <- start0-2
connList <- connListCalc(p)

Real.Data.Sep.Obj <- Real.Data.Separate.Main()

Sep.Est <- Real.Data.Sep.Obj$Sep.Est
Pred.Err.Sep <- Real.Data.Sep.Obj$Pred.Err.Sep
Sep.Stability <- Real.Data.Sep.Obj$Sep.Stability
Sigma.Inv.est <- Real.Data.Sep.Obj$Sigma.Inv.est
L1.est <- Real.Data.Sep.Obj$L1.est
TVE1.est <- Real.Data.Sep.Obj$TVE1.est


####################################################
#### JOINT ESTIMATION(PAIRWISE APPROACH):
#### ESTIMATING ALL 6 PAIRS OF STATES JOINTLY,
#### COMBINING THE ESTIMATES
####################################################


t0 <- start0-2

Real.Data.Joint.Obj <- Real.Data.Joint.Main()
Joint.Est <- Real.Data.Joint.Obj$Joint.Est
Joint.Stability <- Real.Data.Joint.Obj$Joint.Stability
Pred.Err.Joint <- Real.Data.Joint.Obj$Pred.Err.Joint


### Plotting all the estimates

pdf(paste("All_Four_Together_ESTIMATES_","p=",p,"t=",train,"rep=",rep,"start0=",start0,"situation=",sit,"_K=",K,".pdf"),
      width=12)
  par(mfrow=c(1,2))
  for(run in 1:rep){
    for(i in 1:K){
      color2D.matplot(Sep.Est[[run]][[i]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                      show.legend=TRUE,
                      show.values=TRUE,
                      main = paste(c("Separate Estimate",state[i])))

      color2D.matplot(Joint.Est[[run]][[1]][[i]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                      show.legend=TRUE,
                      show.values=TRUE,
                      main = paste(c("Joint Estimate",state[i])))
    }
  }

  par(mfrow=c(1,(K-1)))

  for(run in 1:rep){
   for(i in 1:K){
    for(l in 1:(K-1)){
      color2D.matplot(Joint.Est[[run]][[l]][[i]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                      show.legend=TRUE,
                      show.values=TRUE,
                      main = paste(c("Joint Estimate",state[i])))
    }
   }
  }

  dev.off()


  ### Plotting all the stability matrices
  
pdf(paste("All_Four_Together_STABILITY_","p=",p,"t=",train,"rep=",rep,"start0=",start0,"situation=",sit,"_K=,",K,"df_sep=",df.sep,".pdf"),
    width=12)
par(mfrow=c(1,2))
for(i in 1:K){
  simple.plot(Sep.Stability[[i]]/rep,rown = situation.names(sit),state=state[i],meth="Separate")
  simple.plot(Joint.Stability[[i]]/(rep*(K-1)),rown = situation.names(sit),state=state[i],meth="Joint")
}
dev.off()


par(mfrow=c(1,1))


## Estimated number of factors and total variance explained for each state for each replicate
print("Number of factors:")
print(L1.est)

print("Total variance explained:")
print(TVE1.est)

print("MSFE for separate:")
print(Pred.Err.Sep)
print(mean(Pred.Err.Sep))
print(sd(Pred.Err.Sep))

print("MSFE for joint:")
print(Pred.Err.Joint)
print(mean(Pred.Err.Joint))
print(sd(Pred.Err.Joint))



print("MSFE for separate:")
cat(round(Pred.Err.Sep,3),sep=",")

print("MSFE for joint:")
cat(round(Pred.Err.Joint,3),sep=",")

