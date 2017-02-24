#setwd("/home/usdandres/Documents/Study_stuff/George/Joint Granger Project/Research20/NEW STUFF, FEBRUARY 13TH")
# Copyright 2017, Andrey Skripnikov, All rights reserved.
source("Paper.functions.R")

#set.seed(2)

###########
## Code for simulated experiment: 
## comparing joint and separate approaches for estimating TWO simulated Granger networks.
## (code needs some adjustments for higher number K of networks)
##
## Here is the layout of the code:
##   1. User-defined parameters(all variable descriptions are given next to parameters)
##   2. Data generation(both stationary transition matrices, 
##                      both error covariance matrices that follow K-factor model, 
##                      both simulated time series following the VAR model)
##   3. Separate estimation.
##   4. Joint estimation.
##   5. Printing all the results(details of each performance measure are given
##                               in the end)
##
## Separate estimation is done in Simul.Data.Separate.Main() function
## that takes no arguments but uses global variables.
## Summary of what it does: 
##     - initializes both error covariance matrices with identities,
##     - estimates transition matrix for each entity separately via l1-estimation 
##   and AIC criterion for tuning parameter selection,
##     - uses the transition matrix estimates to get residuals,
##     - use residuals as data to estimate error inverse covariance matrices 
##   help of L-factor model and graphical lasso separately with,
##     - plug-in the resulting inverse error covariance estimates and re-estimate
##   transition matrices separately for each entity via l1-estimation and AIC
##   criterion for tuning parameter selection,
##   For more details on the code - see Paper.functions.R
##
## Joint estimation is done in Simul.Data.Joint.Main() function
## that takes no arguments but uses global variables.
## 
## Summary of what it does: 
##     - uses the estimates of inverse error covariances from the separate estimation process
##       (just plugs them in)
##     - combines both entites into one standardized dataset, sets up the generalized
##       fused lasso as discussed in paper with two tuning parameters lambda1 and lambda2
##     - does sequential search of tuning parameter values: fix lambda2, grid search for lambda1.est,
##       fix lambda1=lambda1.est, grid search for lambda2.est, fix lambda2=lambda2.est etc etc
##     - criterion used for search: AIC distinct(described in Paper.functions.R)
##     - the fusion of parameters is done with help of ADMM algorithm, described
##       in detail in the paper and in Paper.functions.R,
##     - refitting the estimate using the structure(skeleton) of the original
##     estimate: get rid of columns of data matrix(that COMBINES entities) 
##     that correspond to zero elements in the original estimate,
##     fit a regression with a tiny sparsity penalty
##     (to shrink elements a little bit, but so that none of them are 0)


## For more details on how the code functions work see Paper.function.R.


########################################
####  USER-DEFINED PARAMETERS         ##
########################################


### General parameters

K <- 2                     # number of entities
p <- 10                     # number of variables per entity
train <- 30                  # number of training points per entity
rep <- 1                       # number of replicates
max_eig <- 0.6            # max eigenvalue of VAR transition matrices
diff = TRUE            # whether to generate same or different matrices for two entities
h <- 1                      # h-step forecasts will be made
test <- train + 1:h         # test time points(to measure forecasting precision) per entity
t <- train + h              # total number of time points per entity
SNR <- 2                    # Signal-to-Noise ratio
L <- 2                     # number of factors for the L-factor model per entity
sigma.iter <- 2           # number of iterations for error covariance estimation procedure
ed <- 0.02                # off-diagonal edge density shared for both entities
comm <- 0.02               # off-diagonal edge density specific for each entity
Results <- matrix(0,rep,10) # performance measures put into one vector:
# Pred.Err.Joint, FP.Joint, FN.Joint, Matt.Coef.Joint,
# Frob.Joint, Pred.Err.Sep, FP.Sep, FN.Sep, 
# Matt.Coef.Sep, Frob.Sep



#### criterions for tuning parameter selection
criter.joint.1 <- "AICc"  # criterion to select sparsity parameter for joint estimates
criter.joint.2 <- "BIC.dist" # criterion to select fusion parameter for joint estimates
criter.sep <- "AICc"       # criterion to select sparsity parameter for separate estimates
df.sep <- 2                     # the multiplier for "degrees of freedom" part of the criterion for SEPARATE


### parameters for glmnet() function performing solution path calculation for standard lasso problem
fail <- 0
fail.thresh <- 0.7
intercept <- FALSE
standardize <- TRUE


### ADMM algorithm parameters ###

rho <- 10                   # ADMM rho penalty parameter
init <- "0"                 # initialize ADMM with Identity("1"),Zero-matrix("0"), or with separate estimates
eps <- 0.001           # stopping criterion
iter.thresh <- 200     # max number of iterations per one run of ADMM algorithm
fused.thresh <- 0.01  # threshold for differences between elements
n.iter <- 1            # number of sequential iterations


### lambda grid parameters  
knots1 <- 15                                           # controls lower bound of lambda1 path
knots2 <- 20                                           # controls lower bound of lambda2 path
L1 <- NULL                                              # length of lambda1 grid(NULL means we calculate lambda1 path automatically)
L2 <- 20                                                # length of lambda2 grid
#lambda.path <- 1*(1/2)^seq(0, knots1, length = L1)     # lambda1 grid
#L1 <- length(lambda.path)
lambda2.path <- c((2*p)*rho*(1/2)^seq(0, knots1, length = L2),0)   # lambda2 grid, making sure we capture the whole specter of fusion:
                                                                   # from identical A11=A22(for high values of lambda2) to dissimilar A11 and A22(low values of lambda2)
L2 <- length(lambda2.path)

### Thresholding parameters

ConstThresh <- TRUE         # whether to use the hard threshold for final estimates
Thresh <- 0.1               # the value of hard threshold for the final estimate
#criter.sep <- "AIC"        # the criterion to use for tuning parameter selection
df <- 2                     # for AIC.dist 


################################################
## Printing out all the main parameter values ##
################################################


print(c("p:",p))
print(c("t:",t))
print(c("Init:",init))
print(c("criter.joint.1:",criter.joint.1))
print(c("criter.joint.2:",criter.joint.2))
print(c("criter.sep:",criter.sep))
print(c("df.sep",df.sep))
print(c("rho:",rho))
print(c("fail.thresh=",fail.thresh))


print(c("L:",L))
print(c("SNR:",SNR))
print(c("diff:",diff))
print(c("rep",rep))
print(c("ConstThresh",ConstThresh))
print(c("Thresh",Thresh))
print(c("l_lambda1.path",L1))
#print(c("lambda1:",lambda.path))
print(c("lambda2:",lambda2.path))
print(c("eps:",eps))
print(c("iter.thresh",iter.thresh))
print(c("ed:",ed))
print(c("comm:",comm))
print(c("max_eig:",max_eig))
print(c("n.iter:",n.iter))



##############################
### Starting the code timer ##
##############################

ptm <- proc.time()


############################
## INITIALIZING VARIABLES ##
############################


A.true <- list()
A11.true <- list()
A22.true <- list()

FP.Sep <- rep(0,rep)
FN.Sep <- rep(0,rep)
TP.Sep <- rep(0,rep)
TN.Sep <- rep(0,rep)
Matt.Coef.Sep <- rep(0,rep)
Frob.Sep <- rep(0,rep)
Pred.Err.Sep <- rep(0,rep)

FP.Joint <- rep(0,rep)
FN.Joint <- rep(0,rep)
TP.Joint <- rep(0,rep)
TN.Joint <- rep(0,rep)
Matt.Coef.Joint <- rep(0,rep)
Frob.Joint <- rep(0,rep)
Pred.Err.Joint <- rep(0,rep)

lambda1.picked.A11 <- rep(0,rep)
lambda1.picked.A22 <- rep(0,rep)
ind1.picked.A11 <- rep(0,rep)
ind1.picked.A22 <- rep(0,rep)


Joint.lambda1.picked <- rep(0,rep)
Joint.lambda2.picked <- rep(0,rep)
Joint.ind1.picked <- rep(0,rep)
Joint.ind2.picked <- rep(0,rep)

lambda1.path <- list()

L.11.est <- rep(0,rep)
L.22.est <- rep(0,rep)
TVE.11.est <- rep(0,rep)
TVE.22.est <- rep(0,rep)

Frob.11 <- rep(0,rep)
Frob.22 <- rep(0,rep)

Frob.11.Inv <- rep(0,rep)
Frob.22.Inv <- rep(0,rep)

Sep.Est <- list()
Joint.Est <- list()

connList <- connListCalc(p) #the list of correspondence for joint estimation

Pred.Err.Sep.True <- numeric(rep)

for(run in 1:rep){
  print(c("run",run))
  
  
  ##########################################
  #### GENERATING SIMULATED DATA    ########
  ##########################################
  repeat{
  fail <- 0
  Simul.Obj <- Gener.Simul.Data()

  A11.true <- Simul.Obj$A.true[[1]]
  A22.true <- Simul.Obj$A.true[[2]]
  A.true[[run]] <- list(A11.true,A22.true)
  
  ## vectorized versions of matrices, by row
  A11.true.vec <- matrix(t(A11.true),1,p^2)
  A22.true.vec <- matrix(t(A22.true),1,p^2)
  A.true.vec <- cbind(A11.true.vec,A22.true.vec)
  
  Sigma.true <- Simul.Obj$Sigma.true
  Sigma.Inv.true <- Simul.Obj$Sigma.Inv.true
  
  GenData.train <- Simul.Obj$GenData.train
  GenData.test <- Simul.Obj$GenData.test
  
  GenData.train[,train]
  GenData.test
  #mean((block.diag(A.true) %*% GenData.train[,train] - GenData.test)^2)

  ##########################################
  ##### SEPARATE ESTIMATION  ###############
  ##########################################
  
  Separate.Est.Obj <- Simul.Data.Separate.Main()
  if (fail == 0) break;
  }
  
  
  lambda1.picked.A11[run] <- Separate.Est.Obj$lambda1.picked[1]
  ind1.picked.A11[run] <- Separate.Est.Obj$ind1.picked[1]
  lambda1.picked.A22[run] <- Separate.Est.Obj$lambda1.picked[2]
  ind1.picked.A22[run] <- Separate.Est.Obj$ind1.picked[2]
  L.11.est[run] <- Separate.Est.Obj$L1.est[1]
  L.22.est[run] <- Separate.Est.Obj$L1.est[2]
  TVE.11.est[run] <- Separate.Est.Obj$TVE1.est[1]
  TVE.22.est[run] <- Separate.Est.Obj$TVE1.est[2]
  Frob.11[run] <- Separate.Est.Obj$Frob[1]
  Frob.22[run] <- Separate.Est.Obj$Frob[2]
  Frob.11.Inv[run] <- Separate.Est.Obj$Frob.Inv[1]
  Frob.22.Inv[run] <- Separate.Est.Obj$Frob.Inv[2]
  Pred.Err.Sep[run] <- Separate.Est.Obj$Pred.Err.Sep
  FP.Sep[run] <- Separate.Est.Obj$FP.Sep
  FN.Sep[run] <- Separate.Est.Obj$FN.Sep
  TP.Sep[run] <- Separate.Est.Obj$TP.Sep
  TN.Sep[run] <- Separate.Est.Obj$TN.Sep
  Matt.Coef.Sep[run] <- Separate.Est.Obj$Matt.Coef.Sep
  Frob.Sep[run] <- Separate.Est.Obj$Frob.Sep
  Sigma.Inv.est <- Separate.Est.Obj$Sigma.Inv.est
  Results[run,6:10] <-   c(Pred.Err.Sep[run],
                           FP.Sep[run],
                           FN.Sep[run],
                           Matt.Coef.Sep[run],
                           Frob.Sep[run])
  
  Sep.Est[[run]] <- Separate.Est.Obj$Sep.Est
  
  
  
   #############################################
  ######   JOINT ESTIMATION  ##################
  #############################################
  
  if (init == "sep") A.init <- c(matrix(Sep.Est[[run]][[1]],byrow=TRUE),matrix(Sep.Est[[run]][[2]],byrow=TRUE))
  if (init == "0") A.init <- rep(0,2*(p^2))
  if (init == "1") A.init <- c(matrix(diag(1,p),byrow=TRUE),matrix(diag(1,p),byrow=TRUE))

  Joint.Est.Obj <- Simul.Data.Joint.Main()


  Pred.Err.Joint[run] <- Joint.Est.Obj$Pred.Err.Joint
  FP.Joint[run] <- Joint.Est.Obj$FP.Joint
  FN.Joint[run] <- Joint.Est.Obj$FN.Joint
  TP.Joint[run] <- Joint.Est.Obj$TP.Joint
  TN.Joint[run] <- Joint.Est.Obj$TN.Joint
  Matt.Coef.Joint[run] <- Joint.Est.Obj$Matt.Coef.Joint
  Frob.Joint[run] <- Joint.Est.Obj$Frob.Joint
  Joint.lambda1.picked[run] <- Joint.Est.Obj$Joint.lambda1.picked
  Joint.lambda2.picked[run] <- Joint.Est.Obj$Joint.lambda2.picked
  Joint.ind1.picked[run] <- Joint.Est.Obj$Joint.ind1.picked
  Joint.ind2.picked[run] <- Joint.Est.Obj$Joint.ind2.picked
  Results[run,1:5] <- Joint.Est.Obj$Results
  lambda1.path[[run]] <- Joint.Est.Obj$lambda.path


  Joint.Est[[run]] <- Joint.Est.Obj$Joint.Est

  A11.true <- Simul.Obj$A.true[[1]]
  A22.true <- Simul.Obj$A.true[[2]]

}


pdf(paste("SIMUL_ESTIMATES_","p=",p,"_t=",train,"_rep=",rep,"_diff=", diff,"_ed=",ed,"_comm=",comm,"_max_eig=",max_eig,"_K=",K,".pdf"),
    width=12)
par(mfrow=c(1,2))
for(run in 1:rep){
color2D.matplot(A.true[[run]][[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("First entity, TRUE")))
color2D.matplot(A.true[[run]][[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("First entity, TRUE")))
color2D.matplot(Sep.Est[[run]][[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("First entity, Separate Estimate")))
color2D.matplot(Joint.Est[[run]][[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("First entity, Joint Estimate")))
color2D.matplot(A.true[[run]][[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("Second entity, TRUE")))
color2D.matplot(A.true[[run]][[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("Second entity, TRUE")))
color2D.matplot(Sep.Est[[run]][[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("Second entity, Separate Estimate")))
color2D.matplot(Joint.Est[[run]][[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("Second entity, Joint Estimate")))
color2D.matplot(Joint.Est[[run]][[1]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("First J",Joint.lambda2.picked[run],
                               Joint.ind2.picked[run])))
color2D.matplot(Joint.Est[[run]][[2]],cs1=c(0,1),cs2=c(1,0),cs3=c(1,0),
                show.legend=TRUE,
                show.values=TRUE,
                main = paste(c("Second J",Joint.lambda2.picked[run],
                               Joint.ind2.picked[run])))

}

par(mfrow=c(1,1))

dev.off()

#### How long it took for code to run

time.taken <- proc.time() - ptm
print(time.taken)

###############################################
#### These results are measurements for each run in comma-separated form,
#### so that I could compile results from different jobs into one vector
#### for each type of measurement(Prediction error(MSFE),FP,FN,Matthews,
#### Frobenius difference, Frobenius difference refitted)
###################################################


print("JOINT ESTIMATION COMMA-SEPARATED RESULTS:")
print("Prediction error:")
paste(Pred.Err.Joint,collapse=",")

print("FP:")
paste(FP.Joint,collapse=",")

print("FN:")
paste(FN.Joint,collapse=",")

print("Matt:")
paste(Matt.Coef.Joint,collapse=",")

print("Frob:")
paste(Frob.Joint,collapse=",")


print("SEPARATE ESTIMATION COMMA-SEPARATED RESULTS:")
print("Prediction error:")
paste(Pred.Err.Sep,collapse=",")

print("FP:")
paste(FP.Sep,collapse=",")

print("FN:")
paste(FN.Sep,collapse=",")

print("Matt:")
paste(Matt.Coef.Sep,collapse=",")

print("Frob:")
paste(Frob.Sep,collapse=",")


##########################################################
#### Mean and SD for each type of measurement(Prediction error(MSFE),FP,
#### FN,Matthews,Frobenius difference, Frobenius difference refitted)
##########################################################

print("MEAN AND SD OF RESULTS FOR JOINT:")
print("Forecast:")
print(mean(Pred.Err.Joint))
print(sd(Pred.Err.Joint))

print("FP")
print(mean(FP.Joint))
print(sd(FP.Joint))

print("FN")
print(mean(FN.Joint))
print(sd(FN.Joint))

print("Matthews coeff:")
print(mean(Matt.Coef.Joint))
print(sd(Matt.Coef.Joint))

print("Frob")
print(mean(Frob.Joint))
print(sd(Frob.Joint))


print("MEAN AND SD OF RESULTS FOR SEPARATE:")
print("Forecast:")
print(mean(Pred.Err.Sep))
print(sd(Pred.Err.Sep))

print("FP")
print(mean(FP.Sep))
print(sd(FP.Sep))

print("FN")
print(mean(FN.Sep))
print(sd(FN.Sep))

print("Matthews coeff:")
print(mean(Matt.Coef.Sep))
print(sd(Matt.Coef.Sep))

print("Frob")
print(mean(Frob.Sep))
print(sd(Frob.Sep))


#####################################################################################
### ALL THE PICKED TUNING PARAMETER VALUES                                         ##
###                                                                                ##
### FOR JOINT: SPARSITY(lambda1) and FUSION(lambda2)                               ##
###                                                                                ##
### FOR SEPARATE: SPARSITY FOR ENTITY 1(lambda1.A11) AND FOR ENTITY 2(lambda1.A22) ##
#####################################################################################


print("JOINT TUNING PARAMETER VALUES:")
print("lambda1's picked")
print(Joint.lambda1.picked)

print("lambda2's picked")
print(Joint.lambda2.picked)

print("JOINT TUNING INDEX SPARSITY:")
print(Joint.ind1.picked)
print("JOINT TUNING INDEX FUSION:")
print(Joint.ind2.picked)

print("SEPARATE TUNING PARAMETER VALUES:")
print("A11 lambda1's picked")
print(lambda1.picked.A11)

print("A22 lambda1's picked")
print(lambda1.picked.A22)



###################################################################
#### ESTIMATED NUMBER OF FACTORS AND TOTAL VARIANCE EXPLAINED  ####
#### MEASURES OF ERROR COVARIANCE ESTIMATION AND ITS INVERSE    ###
#### MEAN and SD of FROBENIUS DIFFERENCE WITH THE TRUE MATRICES ###
###################################################################

print("Number of factors:")
print(c(L.11.est,L.22.est))
print(table(c(L.11.est,L.22.est)))

print("Total variance explained:")
print(mean(c(TVE.11.est,TVE.22.est)))
print(sd(c(TVE.11.est,TVE.22.est)))

print("Frobenius norm difference for Sigma estimates:")
print(mean(Frob.11 + Frob.22)/2)
print(sd((Frob.11+Frob.22)/2))


print("Frobenius norm difference for Inverse Sigma estimates:")
print(mean(Frob.11.Inv + Frob.22.Inv)/2)
print(sd((Frob.11.Inv+Frob.22.Inv)/2))


########################
#######################

print("lambda1's picked")
print(Joint.lambda1.picked)
print("lambda2's picked")
print(Joint.lambda2.picked)

print("JOINT TUNING INDEX SPARSITY:")
print(Joint.ind1.picked)
print("JOINT TUNING INDEX FUSION:")
print(Joint.ind2.picked)


############################
###########################

cat("Output for direct entry into Latex table:")
cat("\n",sep="")
cat("\n",sep="")
cat(" \\begin{tabular} {l} ","p=",p," t=",train," \\\\ ","diff=",diff," \\end{tabular}"," & ",sep="")
cat("\n",sep="")
cat(" \\begin{tabular} {l} ","S"," \\\\ ","J"," \\end{tabular}"," & ",sep="")
cat("\n",sep="")
cat(" \\begin{tabular} {l} ",round(mean(Pred.Err.Sep),2)," (",round(sd(Pred.Err.Sep),2),") \\\\ ",round(mean(Pred.Err.Joint),2)," (",round(sd(Pred.Err.Joint),2),")"," \\end{tabular}"," & ",sep="")
cat("\n",sep="")
cat(" \\begin{tabular} {l} ",round(mean(FP.Sep),2)," (",round(sd(FP.Sep),2),") \\\\ ",round(mean(FP.Joint),2)," (",round(sd(FP.Joint),2),")"," \\end{tabular}"," & ",sep="") 
cat("\n",sep="")
cat(" \\begin{tabular} {l} ",round(mean(FN.Sep),2)," (",round(sd(FN.Sep),2),") \\\\ ",round(mean(FN.Joint),2)," (",round(sd(FN.Joint),2),") "," \\end{tabular}"," & ",sep="") 
cat("\n",sep="") 
cat(" \\begin{tabular} {l} ",round(mean(Matt.Coef.Sep),2)," (",round(sd(Matt.Coef.Sep),2),")  \\\\ ",round(mean(Matt.Coef.Joint),2)," (",round(sd(Matt.Coef.Joint),2),") "," \\end{tabular}"," & ",sep="") 
cat("\n",sep="") 
cat(" \\begin{tabular} {l} ",round(mean(Frob.Sep),2)," (",round(sd(Frob.Sep),2),")  \\\\ ",round(mean(Frob.Joint),2)," (",round(sd(Frob.Joint),2),") "," \\end{tabular}"," & ",sep="") 
cat("\n",sep="") 





