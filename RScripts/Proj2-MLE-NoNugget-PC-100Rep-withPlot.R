## To do:
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

# -- Generating training set --: 
n_sims <- 3
print(paste0('full size: ', n_sims))


# ### -- Parameter configuration --:
# Test parameter configuration
set.seed(111)
AWghtGrid <- 4 + exp(runif(n_sims, log(0.005), log(1.8)))
summary(AWghtGrid)
length(AWghtGrid)
set.seed(222)
shapeGrid <- runif(n_sims, 0.2, 0.8)
summary(shapeGrid)
length(shapeGrid)

lambdaGrid <- rep(0.00001, length.out=n_sims)
length(lambdaGrid)

parameterComb <- cbind(shapeGrid, AWghtGrid, lambdaGrid)
dim(parameterComb)

# ### -- Generate Test Set for comparison study with likelihood estimation --:
m <- 30 # no of replication in the spatial field
R <- 1
extremalFields <- array(NA,
                        dim=c(n_sims, R, n, m))
dim(extremalFields)
set.seed(123)
tic()
for (i in 1:n_sims)
{
  cat('GEV parameter loop no', i, '\n')

  # Defining GEV parameters:
  shape_val <- parameterComb[i,1]
  a_wght <- parameterComb[i,2]

  # Setup LKrig
  LKinfo <- LKrigSetup(s,
                       a.wght=a_wght,
                       nlevel=1,
                       nu=1,
                       NC=16,  # Changed from 8
                       NC.buffer=4,
                       fixedFunction=NULL)
  tic()
  for(j in 1:R)
  {
    cat('For rep', j, '\n')
    # Simulate SAR coefficients
    coefSARList <- LKrigSAREvd(LKinfo,
                               loc=1,
                               scale=shape_val,
                               shape=shape_val,
                               M=m,
                               asList=FALSE)

    coefSAR <- coefSARList$coefSAR

    # Simulate field
    PHI <- LKrig.basis(s, LKinfo)
    extremalFields[i, j, ,] <- (PHI %*% coefSAR)
  }
  toc()
}
toc()
dim(extremalFields)

### -- Latest -- 
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# dim(parameterComb)
# parameterCombML <- parameterComb
# save(parameterComb,
#      file="smallPC-NoNugget-LikelihoodComparison-Jan29-Test.RData")
# 
# extremalFieldsML <- extremalFields
# dim(extremalFieldsML)
# save(extremalFieldsML,
#      file="small-storeZRep30-NoNugget-LikComp-Jan29-Test.RData")
##################################################################

# ### Loading the test samples:
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# load("small-storeZRep30-NoNugget-LikComp-Jan29-Test.RData")
# dim(extremalFieldsML) # 30 Rep
# 
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# load("smallPC-NoNugget-LikelihoodComparison-Jan29-Test.RData")
# dim(parameterComb)


## Loading the CNN estimates:
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE/CNN-Estimates")
file_path_30 <- path.expand("est_para_testRep30_ML_smalls.npy")

# Import numpy and load the file
np <- import("numpy")
mat30 <- np$load(file_path_30)


# Load the NumPy file
cnnRep30 <- as.array(mat30)
dim(cnnRep30)


# ## -- Maximum Likelihood Estimation --: 
storeMLEs <- array(NA, c(n_sims, R, 2) )# 100   2
dim(storeMLEs)

storeLogLiks  <- matrix(NA, nrow=n_sims, ncol=R)
dim(storeLogLiks)

source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/LatestLikelihoodFunc-Oct7-2024.R")

tic()
for(i in 1:n_sims)
{
  print(paste0('Loop', i))

  # Calculate initial log-likelihood
  para_start <- c(log(parameterComb[i,1]), Awght2Omega(parameterComb[i,2])) # log Xi, log Kappa
  initial_loglike <- log_like(s,
                              NC=16,
                              NC.buff=4,
                              para_start,
                              Z=as.matrix(extremalFields[i, j, ,1:30]))
  cat('Likelihood at initial value: ',  initial_loglike, '\n')
  tic()
  for(j in 1:R)
  {
    print(paste0('Inner Loop:', j))
    ## Run the optimization using L-BFGS-B for faster convergence
    tic()
    # optim_result <- optim(par = para_start,
    #                       fn = log_like,
    #                       s = s,
    #                       NC = 16,
    #                       NC.buff = 4,
    #                       Z = (extremalFields[i, j, , 1:30]-median(extremalFields[i, j, , 1:30]))/sd(extremalFields[i, j, , 1:30]),
    #                       method = 'BFGS',  # Change method if gradients are available
    #                       control = list(reltol = 1e-8,   # Adjust relative tolerance for precision
    #                                      maxit = 100      # Increase maximum iterations if needed
    #                                      ) )
    optim_result <- optim(par=para_start,    # Starting values
                          fn=log_like,       # Log-likelihood function
                          s=s,
                          NC=16,
                          NC.buff=4,
                          Z = (extremalFields[i, j, , 1:30]/IQR(extremalFields[i, j, , 1:30])),
                          # Z=extremalFields[i,j, ,1:30],
                          method='Nelder-Mead',
                          control=list(factr=1e7,
                                       pgtol=1e-1,
                                       maxit=15))
    toc()

    cat('MLE : ', optim_result$par, '\n')
    cat('True : ', para_start, '\n')

    # Store the MLEs and log-likelihood
    storeLogLiks[i,j] <- optim_result$value
    storeMLEs[i,j, ]  <- optim_result$par

    # Check for finite initial log-likelihood
    if (!is.finite(initial_loglike))
    {
      stop(paste('Log-likelihood is not finite at loop',
                 i,
                 j,
                 'with parameters:',
                 para_start))
    }
  }
  toc()
}
toc()
# 3486.677 sec elapsed

# Print the results
print("Log-likelihoods:")
print(storeLogLiks)
print("MLEs:")

# print(dim(storeMLEs))
# 
# # setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
# # storeMLEsSmall <- storeMLEs
# # dim(storeMLEsSmall )
# # save(storeMLEsSmall, file="storeMLEs-smallTest-Jan29.RData")

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/MLE")
load("storeMLEs-smallTest-Jan29.RData")
dim(storeMLEsSmall)

colors <- rep(c("purple", "yellow", "green4"), length.out=300)  # Group colors

class(storeMLEs)

storeMLEsSmall[1, ,1]
plot(storeMLEsSmall[1, ,1])

setwd("~/Desktop")
png("ML-CNN-Est-Small.png",
    units="in", 
    width=27,
    height=8.8,
    res=200)
par(mfrow=c(1,3), 
    mar=c(5,4,4,4) + 0.3,
    oma=c(1.5,3.5,3.5,8))
set.panel(1,3)

## Plot 1
plot(storeMLEsSmall[1, ,2],
     exp(storeMLEsSmall[1, ,1]),
     xlim=c(min(storeMLEsSmall[1, ,2], cnnRep30[1, ,2]),
            max(storeMLEsSmall[1, ,2], cnnRep30[1, ,2])),
     ylim=c(min(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1]),
            max(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1])),
     pch=19,
     cex=3,
     #lwd=5.5,
     col=alpha("#63B", alpha=0.5),
     yaxt='n',
     xaxt='n',
     xlab='',
     ylab= '')
points(cnnRep30[1, ,2],
       cnnRep30[1, ,1],
       xlim=c(min(storeMLEsSmall[1, ,2], cnnRep30[1, ,2]),
              max(storeMLEsSmall[1, ,2], cnnRep30[1, ,2])),
       ylim=c(min(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1]),
         max(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1])),
       pch=19,
       cex=3,
       lwd=5.5,
       col=alpha("#8DD641", alpha=0.5),
       yaxt='n',
       xaxt='n',
       xlab='',
       ylab= '')
points(log(parameterComb[1,2]-4),
       parameterComb[1,1],
       xlim=c(min(storeMLEsSmall[1, ,2], cnnRep30[1, ,2]),
              max(storeMLEsSmall[1, ,2], cnnRep30[1, ,2])),
       ylim=c(min(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1]),
              max(exp(storeMLEsSmall[1, ,1]), cnnRep30[1, ,1])),
       col= "#E8336C" ,
       cex=3,
       lwd=9,
       pch=3,
       yaxt='n',
       xaxt='n',
       xlab='',
       ylab= '')
axis(2,
     las=1,
     cex.axis=2,
     lwd=2.2,
     padj=0.9)
axis(1,
     cex.axis=2,
     lwd=2.2,
     padj=0.9)
mtext(expression(xi),
      side=2.2,
      cex=3,
      line=5.4,
      las=2)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=2.2,
      line=5.8)

## Plot 2
plot(storeMLEsSmall[2, ,2],
     exp(storeMLEsSmall[2, ,1]),
     xlim=c(min(storeMLEsSmall[2, ,2], cnnRep30[2, ,2]),
            max(storeMLEsSmall[2, ,2], cnnRep30[2, ,2])),
     ylim=c(min(exp(storeMLEsSmall[2, ,1]), cnnRep30[2, ,1]),
            max(exp(storeMLEsSmall[2, ,1]), cnnRep30[2, ,1])),
     pch=19,
     cex=3,
     col=alpha("#63B", alpha=0.5),
     lwd=5.5,
     yaxt='n',
     xaxt='n',
     xlab='',
     ylab= '')
points(cnnRep30[2, ,2],
       cnnRep30[2, ,1],
       xlim=c(min(storeMLEsSmall[2, ,2], cnnRep30[2, ,2]),
              max(storeMLEsSmall[2, ,2], cnnRep30[2, ,2])),
       ylim=c(min(exp(storeMLEsSmall[2, ,1]), cnnRep30[2, ,1]),
              max(exp(storeMLEsSmall[2, ,1]), cnnRep30[2, ,1])),
       pch=19,
       cex=3,
       col=alpha("#8DD641", alpha=0.5),
       lwd=5.5,
       yaxt='n',
       xaxt='n',
       xlab='',
       ylab= '')
points(log(parameterComb[2,2]-4),
       parameterComb[2,1],
       xlim=c(min(storeMLEsSmall[2, ,2], cnnRep30[2, ,2]),
              max(storeMLEsSmall[2, ,2], cnnRep30[2, ,2])),
       ylim=c(min(exp(storeMLEsSmall[2, ,1]), cnnRep30[2, ,1]),
              max(exp(storeMLEsSmall[2, ,1]), cnnRep30[2, ,1])),
       col= "#E8336C" ,
       cex=3,
       lwd=9,
       pch=3,
       yaxt='n',
       xaxt='n',
       xlab='',
       ylab= '')
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2.2,
      cex=3,
      line=5.4,
      las=2)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=2.2,
      line=5.8)


## Plot 3
plot(storeMLEsSmall[3, ,2],
     exp(storeMLEsSmall[3, ,1]),
     xlim=c(min(storeMLEsSmall[3, ,2], cnnRep30[3, ,2]),
            max(storeMLEsSmall[3, ,2], cnnRep30[3, ,2])),
     ylim=c(min(exp(storeMLEsSmall[3, ,1]), cnnRep30[3, ,1]),
            max(exp(storeMLEsSmall[3, ,1]), cnnRep30[3, ,1])),
     pch=19,
     cex=3,
     col=alpha("#63B", alpha=0.5),
     lwd=5.5,
     yaxt='n',
     xaxt='n',
     xlab='',
     ylab= '')
points(cnnRep30[3, ,2],
       cnnRep30[3, ,1],
       xlim=c(min(storeMLEsSmall[3, ,2], cnnRep30[3, ,2]),
              max(storeMLEsSmall[3, ,2], cnnRep30[3, ,2])),
       ylim=c(min(exp(storeMLEsSmall[3, ,1]), cnnRep30[3, ,1]),
              max(exp(storeMLEsSmall[3, ,1]), cnnRep30[3, ,1])),
       pch=19,
       cex=3,
       col=alpha("#8DD641", alpha=0.5),
       lwd=5.5,
       yaxt='n',
       xaxt='n',
       xlab='',
       ylab= '')
points(log(parameterComb[3,2]-4),
       parameterComb[3,1],
       xlim=c(min(storeMLEsSmall[3, ,2], cnnRep30[3, ,2]),
              max(storeMLEsSmall[3, ,2], cnnRep30[3, ,2])),
       ylim=c(min(exp(storeMLEsSmall[3, ,1]), cnnRep30[3, ,1]),
              max(exp(storeMLEsSmall[3, ,1]), cnnRep30[3, ,1])),
       col= "#E8336C" ,
       cex=3,
       lwd=9,
       pch=3,
       yaxt='n',
       xaxt='n',
       xlab='',
       ylab= '')
axis(2,
     las=1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
axis(1,
     cex.axis=2.2,
     lwd=2,
     padj=0.9)
mtext(expression(xi),
      side=2.2,
      cex=3,
      line=5.4,
      las=2)
mtext(expression('log('~kappa~')'),
      side=1,
      cex=2.2,
      line=5.8)
dev.off()






