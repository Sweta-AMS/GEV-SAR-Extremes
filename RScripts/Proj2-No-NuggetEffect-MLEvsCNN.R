rm(list=ls())

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")

source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")

## Location
M <- 16
sGrid <- list(x=1:M, y=1:M)
s <- make.surface.grid(sGrid)
plot(s, pch=19)


n <- nrow(s)
print(paste0('Field size: ', n))

# ## -- Generating training set --: 
# n_grid <- 6
# n_sims <- n_grid*n_grid
# print(paste0('full size: ', n_sims))
 
# # ### -- Parameter configuration --:
# # # Test parameter configuration
# AWghtGrid <- 4 + exp(seq(log(0.005), log(1.8), length.out=n_grid))
# summary(AWghtGrid)
# length(AWghtGrid)
# 
# shapeGrid <- seq(0.2, 0.8, length.out=n_grid)
# summary(shapeGrid)
# length(shapeGrid)
# 
# lambdaGrid <- rep(0.00001, length.out=n_sims)
# length(lambdaGrid)
# 
# parameterComb <- expand.grid(shapeGrid, AWghtGrid)
# dim(parameterComb)
# 
# parameterComb <- cbind(parameterComb, lambdaGrid)
# dim(parameterComb)
# 
# 
# ## -- Generation of test sample, withtout replication --:
# m <- 30 # no of replication in the spatial field
# R <- 60
# extremalFields <- array(NA,
#                         dim=c(n_sims, R, n, m))
# dim(extremalFields)
# 
# set.seed(123)
# tic()
# for (i in 1:n_sims)
# {
#   cat('GEV parameter loop no', i, '\n')
# 
#   # Defining GEV parameters:
#   shape_val <- parameterComb[i,1]
#   a_wght <- parameterComb[i,2]
# 
#   # Setup LKrig
#   LKinfo <- LKrigSetup(s,
#                        a.wght=a_wght,
#                        nlevel=1,
#                        nu=1,
#                        NC=16,  # Changed from 8
#                        NC.buffer=4,
#                        fixedFunction=NULL)
#   for(j in 1:R)
#   {
#     cat('Rep no', j, '\n')
#     # Simulate SAR coefficients
#     coefSARList <- LKrigSAREvd(LKinfo,
#                                loc=1,
#                                scale=shape_val,
#                                shape=shape_val,
#                                M=m,
#                                asList=FALSE)
#     coefSAR <- coefSARList$coefSAR
#     
#     # Simulate field
#     PHI <- LKrig.basis(s, LKinfo)
#     extremalFields[i,j,  ,] <- (PHI %*% coefSAR)
#   }
# }
# toc()
# dim(extremalFields)


M <- 4
sGrid <- list(x=1:M, y=1:M)
s <- make.surface.grid(sGrid)
plot(s, pch=19)

LKinfo <- LKrigSetup(s,
                     a.wght=4.01,
                     nlevel=1,
                     nu=1,
                     NC=4, 
                     NC.buffer=0,
                     fixedFunction=NULL)
B <- LKrigSAR(LKinfo, Level=1)
B <- spind2spam(B)
dim(B)



setwd("~/Desktop")
png("SAR-example.png",
    units="in", 
    width=7,
    height=3.6,
    res=200)
par(mfrow=c(1,2),
    mar=c(2,2,1.5,2.5),
    oma=c(2.5,2.5,1.5,2.5))
set.panel(1,2)
plot(s,
     pch=19,
     xlab = expression(s[1]),
     ylab = expression(s[2]), main='4 x 4 grid')
matB <- as.matrix(B)
image.plot(
  1:ncol(matB), 
  1:nrow(matB),
  matB,
  xlab = "Column", ylab = "Row",
  main = "SAR Matrix B",
  col = viridis::viridis(100) ) # or your palette of choice
dev.off()

## Working directory:
# save(parameterComb,
#      file="TruePar-NoNugget-Test-small.RData")

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
load("TruePar-NoNugget-Test-small.RData")
dim(parameterComb)

load("storeZRep30-NoNugget-Test-small.RData")
dim(extremalFields)

# ### -- Maximum Likelihood Estimation --: 
# R <- 60
# storeMLEs30 <- array(NA,
#                     dim=c(n_sims,R,2))# 100   2
# dim(storeMLEs30)
# 
# source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/LatestLikelihoodFunc-Oct7-2024.R")
# tic()
# for(i in 1:n_sims) # n_sims
# {
#   print(paste0('Loop', i))
#   
#   # Parameter transformation
#   para_start <- c(log(parameterComb[i,1]),
#                   Awght2Omega(parameterComb[i,2]))
#   
#   for(j in 1:R) {
#   print(paste0('Loop j', j))
#   # Calculate initial log-likelihood
#   initial_loglike <- tryCatch({
#     log_like(s,
#              NC=16,
#              NC.buff=4,
#              para_start,
#              Z=as.matrix(extremalFields[i,j,  ,1:30]))
#   }, error = function(e) {
#     cat("Error in log_like at iteration", i, ":", conditionMessage(e), "\n")
#     return(Inf)  # Assign a large value so that optimization avoids it
#   })
#   
#   cat('Likelihood at initial value: ',  initial_loglike, '\n')
#   
#   # Skip optimization if log-likelihood is not finite
#   if (!is.finite(initial_loglike)) {
#     cat("Skipping iteration", i, "due to invalid initial log-likelihood.\n")
#     next
#   }
#   
#   # Run the optimization
#   optim_result <- optim(par=para_start,
#                         fn=log_like,
#                         s=s,
#                         NC=16,
#                         NC.buff=4,
#                         Z=as.matrix(extremalFields[i,j,  ,1:30]),
#                         method= 'Nelder-Mead',
#                         # 'L-BFGS-B',
#                         #
#                         control=list(factr=1e7,
#                                      pgtol=1e-1,
#                                      maxit=15))
#   cat('MLE : ', optim_result$par, '\n')
#   cat('True : ', para_start, '\n')
#   
#   storeMLEs30[i, j,  ]  <- optim_result$par
#   }
# }
# toc()
# dim(storeMLEs30)

## -- Save the MLEs across Rep30 --:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# save(storeMLEs30,
#     file="storeMLE-for-Rep30-Repetition60-Sample-6x6grid-Latest.RData")
################################################################################
# save(storeMLEs1, file="storeMLE-for-Rep1-Sample-32x32grid-Latest.RData")
# save(storeMLEs10, file="storeMLE-for-Rep10-Sample-32x32grid-Latest.RData")
# save(storeMLEs30, file="storeMLE-for-Rep30-Sample-32x32grid-Latest.RData")

# load("storeMLE-for-Rep1-Sample-32x32grid-Latest.RData")
# load("storeMLE-for-Rep10-Sample-32x32grid-Latest.RData")
# load("storeMLE-for-Rep30-Sample-32x32grid-Latest.RData")
################################################################################
load('storeMLE-for-Rep30-Repetition60-Sample-6x6grid-Latest.RData')

### MSE: 
dim(storeMLEs30) # 36 60  2

xi_mle <- exp(storeMLEs30[, , 1]) # 36 60: exp( )
kappa_mle <- sqrt(exp(2*storeMLEs30[ , , 2]))  # 36 60:  (Awght2Omega)

xi_true <- (parameterComb[ ,1])
kappa_true <- sqrt(parameterComb[ ,2]-4)


## Load  CNN estimation:
# load("storeMLE-for-Rep1-Sample-32x32grid-Latest.RData")
file_path_30 <- path.expand("~/Downloads/cnnEst-small-testSample.npy")

# Import numpy and load the file
np <- import("numpy")
mat30 <- np$load(file_path_30)

# Load the NumPy file
cnnRep30 <- (mat30)
dim(cnnRep30)

xi_cnn <- (cnnRep30[,,1])
kappa_cnn <-  sqrt(exp(cnnRep30[,,2]))

#–– 3. Compute overall MSEs ––
xi_true_mat <- matrix(xi_true, 
                      nrow=36,
                      ncol=60, 
                      byrow=FALSE)
kappa_true_mat <- matrix(kappa_true,
                         nrow=36,
                         ncol=60, 
                         byrow=FALSE)
# 2. Compute squared errors
sq_err_xi_ml <- (xi_mle - xi_true_mat)^2   # still 36×60
sq_err_kappa_ml <- (kappa_mle - kappa_true_mat)^2   # still 36×60

sq_err_xi_cnn <- (xi_cnn - xi_true_mat)^2   # still 36×60
sq_err_kappa_cnn <- (kappa_cnn - kappa_true_mat)^2   # still 36×60

# 3. Compute MSE per parameter (mean over the 60 reps)
rmse_xi_ml <- (rowMeans(sq_err_xi_ml))^0.5   # length‑36
rmse_kappa_ml <- (rowMeans(sq_err_kappa_ml))^0.5   # length‑36

rmse_xi_cnn <- (rowMeans(sq_err_xi_cnn))^0.5 # length‑36
rmse_kappa_cnn <- (rowMeans(sq_err_kappa_cnn))^0.5   # length‑36

# 4.  Plotting:
setwd("~/Desktop")
png("MLvsCNN-MSE-Rep30.png",
    units="in", 
    width=7,
    height=3.6,
    res=200)
par(mfrow=c(1,2),
    mar=c(2,2,1.5,2.5),
    oma=c(2.5,2.5,1.5,2.5))
set.panel(1,2)

bplot.xy((parameterComb[,1]),
         log(rmse_xi_ml/rmse_xi_cnn))
         # ylim=c(0.1, 40))
abline(h=0, col='magenta', lty=2, lwd=1)
mtext(expression("True "),
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext(expression(log(RMSE[MLE]/RMSE[CNN])),
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression(xi),
      side=3,
      line=1)
bplot.xy(sqrt(parameterComb[,2]-4),
         log(rmse_kappa_ml/rmse_kappa_cnn))
         # ylim=c(0,4))
abline(h=0, col='magenta', lty=2, lwd=1)
mtext( expression("True "),
      side=1,
      cex=0.9,
      line=3,
      las=1)
mtext(expression(log(RMSE[MLE]/RMSE[CNN])),
      side=2,
      cex=0.9,
      las=3,
      line=2.5)
mtext(expression(kappa),
      side=3,
      line=1)
dev.off()


################################################################################
set.panel(1,3)
bplot.xy(parameterComb[,1],
         exp(storeMLEs30[,1, 1])-parameterComb[,1],
         col=alpha('grey', 0.1),
        # ylim=c(-0.2,4),
         #ylim=c(0.2,1.6),
         pch=19,
         xlab='True',
         ylab='MLEs',
         main='Rep 1')
#ylim=c(0.2,1),
#main=expression(xi))
abline(h=0, col='magenta')

bplot.xy(parameterComb[,1],
         exp(storeMLEs10[,1])-parameterComb[,1],
         col=alpha('grey', 0.1),
         ylim=c(-0.2,4),
         #ylim=c(0.2,1.6),
         pch=19,
         xlab='True',
         ylab='MLEs',
         main='Rep 10')
#ylim=c(0.2,1),
# main=expression(xi))
#abline(0,1, col='magenta')
abline(h=0, col='magenta')

bplot.xy(parameterComb[,1],
         exp(storeMLEs30[,1])-parameterComb[,1],
         col=alpha('grey', 0.1),
         ylim=c(-0.2,4),
         #ylim=c(0.2,1.6),
         main='Rep 30',
         pch=19,
         xlab='True',
         ylab='MLEs')
#ylim=c(0.2,1),
#main=expression(xi))
abline(h=0, col='magenta')

set.panel(1,3)
bplot.xy(Awght2Omega(parameterComb[,2]),
         storeMLEs1[,2]-Awght2Omega(parameterComb[,2]),
         col=alpha('grey', 0.1),
         ylim=c(-0.2,4),
         #ylim=c(0.2,1.6),
         pch=19,
         xlab='True',
         ylab='MLEs',
         main='Rep 1')
#ylim=c(0.2,1),
#main=expression(xi))
abline(h=0, col='magenta')

bplot.xy(Awght2Omega(parameterComb[,2]),
         storeMLEs10[,2]-Awght2Omega(parameterComb[,2]),
         col=alpha('grey', 0.1),
         ylim=c(-0.2,4),
         #ylim=c(0.2,1.6),
         pch=19,
         xlab='True',
         ylab='MLEs',
         main='Rep 10')
#ylim=c(0.2,1),
# main=expression(xi))
#abline(0,1, col='magenta')
abline(h=0, col='magenta')

bplot.xy(Awght2Omega(parameterComb[,2]),
         storeMLEs30[,2]-Awght2Omega(parameterComb[,2]),
         col=alpha('grey', 0.1),
         ylim=c(-0.2,4),
         #ylim=c(0.2,1.6),
         main='Rep 30',
         pch=19,
         xlab='True',
         ylab='MLEs')
#ylim=c(0.2,1),
#main=expression(xi))
abline(h=0, col='magenta')


abline(0,1, col='magenta')

# source Awght2Omega function:
bplot.xy(Awght2Omega(parameterComb[,2]),
         storeMLEs1[,2]-Awght2Omega(parameterComb[,2]),
         #ylim=c(-0.2,4),
         ylim=c(-0.4,6),
         pch=19,
         xlab='True',
         ylab='MLEs',
         col=alpha('grey', 0.6),
         main=expression('log'~ kappa))
abline(h=0, col='magenta')
bplot.xy(Awght2Omega(parameterComb[,2]),
         storeMLEs10[,2]-Awght2Omega(parameterComb[,2]),
         ylim=c(-0.4,6),
         #ylim=c(-0.2,4),
         pch=19,
         xlab='True',
         ylab='MLEs',
         col=alpha('grey', 0.6),
         main=expression('log'~ kappa))
abline(h=0, col='magenta')
abline(0,1, col='magenta')


  
  


