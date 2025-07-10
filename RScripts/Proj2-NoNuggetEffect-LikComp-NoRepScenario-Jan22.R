## To do:
# Let's see if first the mles are close to true value: seems like it works
# then train the CNN on no nugget case to estimate the parameter using CNN
# bias correction and percentile bootstrap

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

# # -- Generating training set --: 
n_grid <- 32
n_sims <- n_grid*n_grid
print(paste0('full size: ', n_sims))


# ### -- Parameter configuration --:
# # Test parameter configuration
AWghtGrid <- 4 + exp(seq(log(0.005), log(1.8), length.out=n_grid))
summary(AWghtGrid)
length(AWghtGrid)

shapeGrid <- seq(0.2, 0.8, length.out=n_grid)
summary(shapeGrid)
length(shapeGrid)

lambdaGrid <- rep(0.00001, length.out=n_sims)
length(lambdaGrid)

parameterComb <- expand.grid(shapeGrid, AWghtGrid)
dim(parameterComb)

parameterComb <- cbind(parameterComb, lambdaGrid)
dim(parameterComb)


# ## -- Generation of test sample, withtout replication --:
# m <- 30 # no of replication in the spatial field
# extremalFields <- array(NA,
#                         dim=c(n_sims, n, m))
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
#   
#   # Simulate SAR coefficients
#   coefSARList <- LKrigSAREvd(LKinfo,
#                              loc=1,
#                              scale=shape_val,
#                              shape=shape_val,
#                              M=m,
#                              asList=FALSE)
#     
#   coefSAR <- coefSARList$coefSAR
#     
#   # Simulate field
#   PHI <- LKrig.basis(s, LKinfo)
#   extremalFields[i, , ] <- (PHI %*% coefSAR)
# }
# toc()
# dim(extremalFields)
# image.plot(matrix(extremalFields[500, ,2], 16,16))
# 
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# dim(parameterComb) # parameter across grid: Xi and Kappa2
# save(parameterComb,
#      file="TruePar-NoNugget-LikComp-01-28-25-Test.RData")
# 
# 
# dim(extremalFields)
# save(extremalFields,
#      file="storeZRep30-NoNugget-LikComp-01-28-25-Test.RData")

## Load the test sample generated across the grid-32x32:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
load("TruePar-NoNugget-LikComp-01-28-25-Test.RData")
dim(parameterComb)

load("storeZRep30-NoNugget-LikComp-01-28-25-Test.RData")
dim(extremalFields)


# ## -- Maximum Likelihood Estimation --: 
storeMLEs30 <- matrix(NA,
                    nrow=n_sims,
                    ncol=2)# 100   2
dim(storeMLEs30)

storeLogLiks30  <- rep(NA, length.out=n_sims)
length(storeLogLiks30)

source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/LatestLikelihoodFunc-Oct7-2024.R")
#source("~/Desktop/Research/Tentative Results-Proj2/R script/LatestLikelihoodFunc-Oct7-2024.R")
tic()
for(i in 1:10) # n_sims
{
  print(paste0('Loop', i))

  # Parameter transformation
  para_start <- c(log(parameterComb[i,1]),
                  Awght2Omega(parameterComb[i,2]))

  # Calculate initial log-likelihood
  initial_loglike <- tryCatch({
    log_like(s,
             NC=16,
             NC.buff=4,
             para_start,
             Z=as.matrix(extremalFields[i, ,1:30]))
    }, error = function(e) {
    cat("Error in log_like at iteration", i, ":", conditionMessage(e), "\n")
    return(Inf)  # Assign a large value so that optimization avoids it
  })

  cat('Likelihood at initial value: ',  initial_loglike, '\n')

  # Skip optimization if log-likelihood is not finite
  if (!is.finite(initial_loglike)) {
    cat("Skipping iteration", i, "due to invalid initial log-likelihood.\n")
    next
  }

  # Run the optimization
  optim_result <- optim(par=para_start,
          fn=log_like,
          s=s,
          NC=16,
          NC.buff=4,
          Z=extremalFields[i, ,1:30],
          method= 'L-BFGS-B',
            # 'Nelder-Mead',
            # 'L-BFGS-B',
            #
          control=list(factr=1e7,
                       pgtol=1e-1,
                       maxit=15))

  cat('MLE : ', optim_result$par, '\n')
  cat('True : ', para_start, '\n')

  # Store the MLEs and log-likelihood
  storeLogLiks30[i] <- optim_result$value
  storeMLEs30[i, ]  <- optim_result$par
}
toc()
# # 8 hrs  40 mins elapsed time
# dim(storeMLEs30) # 1024x2
# 



# # Print the results
# print("Log-likelihoods:")
# print(length(storeLogLiks30))
# 
# print("MLEs:")
# print(dim(storeMLEs30))

## -- Save the MLEs across Rep30 --:
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# save(storeMLEs1, file="storeMLE-for-Rep1-Sample-32x32grid-Latest.RData")
# save(storeMLEs10, file="storeMLE-for-Rep10-Sample-32x32grid-Latest.RData")
# save(storeMLEs30, file="storeMLE-for-Rep30-Sample-32x32grid-Latest.RData")

load("storeMLE-for-Rep1-Sample-32x32grid-Latest.RData")
load("storeMLE-for-Rep10-Sample-32x32grid-Latest.RData")
load("storeMLE-for-Rep30-Sample-32x32grid-Latest.RData")


# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# save(storeLogLiks, file="storeLogLiks-for-Rep30-Sample-32x32grid.RData")
#################################

# # ## -- Save the MLEs across Rep10 --:
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# storeMLEsRep10 <- storeMLEs
# save(storeMLEsRep10, file="storeMLE-for-Rep10-Sample-32x32grid.RData")
# 
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# storeLogLiks10 <- storeLogLiks
# save(storeLogLiks10, file="storeLogLiks-for-Rep10-Sample-32x32grid.RData")
#################################

# ## -- Save the MLEs across Rep1 --:
# storeMLEsRep1 <- storeMLEs1
# 
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# save(storeMLEs1, file="storeMLE-for-Rep1-Sample-32x32grid.RData")
# 
# storeLogLiksRep1 <- storeLogLiks1
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# save(storeLogLiks1, file="storeLogLiks-for-Rep1-Sample-32x32grid.RData")
#################################

# load("storeMLE-for-Rep10-Sample-32x32grid.RData")

# colors <- rep(c("purple", "yellow", "green4"), length.out=300)  # Group colors
# 
# class(storeMLEs)
# 
# # Base Scatter Plot
# plot((storeMLEs[,1]),
#      (storeMLEs[,2]), 
#      col=colors,
#      pch=16, 
#      xlab=expression(sigma[epsilon]), 
#      ylab=expression(rho),
#      main="")
# 
# # Adding the red cross (mean point)
# mean_x <- mean(sigma_epsilon)
# mean_y <- mean(rho)
# points(mean_x, mean_y, col="red", pch=3, cex=2, lwd=3)  # Cross Marker

set.panel(1,3)
bplot.xy(parameterComb[,1],
         exp(storeMLEs1[,1])-parameterComb[,1],
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

#### --- Profile Likelihood ---:
parameterComb 
# Define a range (grid) for the parameter of interest, say theta
theta_grid <- seq(from = min(Awght2Omega(parameterComb[,2])), 
                  to   = max(Awght2Omega(parameterComb[,2])), 
                  length.out = 32)

# Prepare a place to store the profile likelihood results
profile_lik <- numeric(length(theta_grid))
profile_phi <- numeric(length(theta_grid))  # Store corresponding phi values

# Loop over your grid of theta values
for(k in seq_along(theta_grid)) {
  cat("Profiling for theta =", theta_grid[k], "\n")
  
  # For each fixed theta, define an objective function that only depends on phi
  profile_fn <- function(phi,
                         s, 
                         NC,
                         NC.buff,
                         Z,
                         fixed_theta) {
    # Combine the fixed theta with the current phi
    par_vec <- cbind(phi, theta_grid)
    
    log_like(s,
             NC = NC,
             NC.buff = NC.buff, 
             par_vec,
             Z = Z)
  }
  
  # Use an initial guess for phi (nuisance parameter). You might use the value from parameterComb or a constant.
  phi_start <- log(parameterComb[1,1])  # For example
  
  # Optimize over phi while holding theta fixed
  opt_result <- optim(par = phi_start,
                      fn = profile_fn,
                      s = s,
                      NC = 16,
                      NC.buff = 4,
                      Z = extremalFields[1, , 1:30],  # or loop over simulations as needed
                      fixed_theta = theta_grid[k],
                      method = 'Nelder-Mead')
                      # control = list(factr = 1e7,
                      #                pgtol = 1e-1,
                      #                maxit = 50))
  
  # Record the optimized value and the corresponding likelihood
  profile_phi[k] <- opt_result$par
  profile_lik[k] <- opt_result$value
  cat("Optimized phi =", opt_result$par, "with profile likelihood =", opt_result$value, "\n")
}

# Now you can plot the profile likelihood as a function of theta:
plot(theta_grid, profile_lik, type = "b", xlab = "theta (fixed)", ylab = "Profile Log-Likelihood")


