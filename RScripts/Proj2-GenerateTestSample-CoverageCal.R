## Simulate the LK fields and compare with matern gaussian field: 
rm(list=ls())

## Required script
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")

## Generating samples with replication for each parameter configuration
## Location
M <- 16 # changed from 15
sGrid <- list(x=1:M, y=1:M)
s <- make.surface.grid(sGrid)
plot(s, pch=19)

## Generating the training set    
## training grid size 50 and validation is 20, let's keep testing 10
n <- nrow(s)
print(paste0('Field size: ', n))

## Latest used
n_sims <- 7*7*7 # changed 1000
print(paste0('full size: ', n_sims))


# ## -- Parameter configuration --:
# # # Test parameter configuration
set.seed(111)
# AWghtGrid <- 4 + exp(runif(n=n_sims, log(0.005), log(1.8))) # (0.001, 2)
AWghtGrid <- 4 + exp(seq(log(0.005), log(1.8), length.out=7)) # (0.001, 2)
summary(AWghtGrid)
length(AWghtGrid)

# set.seed(222)
# lambdaGrid <- exp(runif(n=n_sims, log(0.00001), log(0.05))) # (0.00001, 0.099)
lambdaGrid <- exp(seq(log(0.00001), log(0.05), length.out=7)) # (0.00001, 0.099)
# lambdaGrid <- rep(0.001, length.out=40)
summary(lambdaGrid)
length(lambdaGrid)

set.seed(333)
# shapeGrid <- runif(n=n_sims, 0.2, 0.8) # (0.01, 0.9)
shapeGrid <- seq(0.2, 0.8, length.out=7) # (0.01, 0.9)
summary(shapeGrid)
length(shapeGrid)

# parameterComb <- cbind(shapeGrid, AWghtGrid, lambdaGrid)
# dim(parameterComb)
parameterComb <- expand.grid(shapeGrid, AWghtGrid, lambdaGrid)
dim(parameterComb)


# # set.seed(111)
# # AWghtGrid <- 4 + exp(seq(log(0.005), log(1.8), length.out=n_sims)) # (0.001, 2)
# # summary(AWghtGrid)
# # hist(AWghtGrid)
# # 
# # set.seed(222)
# # lambdaGrid <- exp(seq(log(0.00001), log(0.05), length.out=n_sims)) # (0.00001, 0.099)
# # summary(lambdaGrid)
# # hist(lambdaGrid)
# # 
# # set.seed(333)
# # shapeGrid <- seq(0.2, 0.8, length.out=n_sims) # (0.01, 0.9)
# # summary(shapeGrid)
# # hist(shapeGrid)
# # 
# # parameterComb <- expand.grid(shapeGrid, 
# #                              AWghtGrid, 
# #                              lambdaGrid)
# # dim(parameterComb)


plot(parameterComb[,1],
     parameterComb[,2],
     pch=19,
     xlab=expression(xi),
     ylab=expression('4 + '~kappa^2),
     col=alpha('grey', alpha=0.6))
plot(parameterComb[,1],
     parameterComb[,3],
     pch=19,
     xlab=expression(xi),
     ylab=expression(tau^2),
     col=alpha('grey', alpha=0.6))
plot(parameterComb[,2],
     parameterComb[,3],
     pch=19,
     xlab=expression('4 + '~kappa^2),
     ylab=expression(tau^2),
     col=alpha('grey', alpha=0.6))

# ### -- GENERATE TRAINING SET --:
m <- 30
R <- 500
extremalFieldsIV <- array(NA,
                        dim=c(43, R, n, m)) ## 121:172 done, next to 173:300 and 300:343
dim(extremalFieldsIV)

tic()  
set.seed(123)
for(i in 1:43) # next should be 101:216
{
  # Defining GEV parameters:
  shape_val <- parameterComb[(i+120+52+128), 1]
  cat('Xi: ', shape_val, '\n')

  a_wght <- parameterComb[(i+120+52+128), 2]
  cat('Kappa^2 : ', a_wght, '\n')

  lambda_val <- parameterComb[(i+120+52+128), 3]
  cat('lambda^2 : ', lambda_val, '\n')

  # Loop across the parameter configuration:
  cat('GEV parameter loop no', (i+120+52+128) , '\n')


  # Setup LKrig
  LKinfo <- LKrigSetup(s,
                       a.wght=a_wght,
                       nlevel=1,
                       nu=1,
                       NC=16,  # Changed from 8
                       NC.buffer=4,
                       fixedFunction=NULL)

  # Add measurement error
  scaleParameter <- log(1 + (lambda_val))  # Variance of the log-normal is tau2
  locParameter <- (-scaleParameter / 2)


  tic()
  for(j in 1:R)
  {

    cat('Loop across replication', j, '\n')

    # Simulate SAR coefficients
    coefSARList <- LKrigSAREvd(LKinfo,
                               loc=1,
                               scale=shape_val,
                               shape=shape_val,
                               M=m,
                               asList=FALSE)
#
    coefSAR <- coefSARList$coefSAR

    # Simulate field
    PHI <- LKrig.basis(s, LKinfo)
    ySim <- (PHI %*% coefSAR)

    # Log-normal matrix generation, unique for each iteration
    LogNormalMat <- exp(matrix(rnorm(n*m,
                                     mean=locParameter,
                                     sd=sqrt(scaleParameter)),
                               n, m))  # n x m matrix

    # Spatial Extremal Fields with Nugget effect
    Z <- ySim*LogNormalMat
    extremalFieldsIV[i, j,  , ] <- Z
  }
  toc()
}
toc()  # End timing for the entire process
dim(extremalFieldsIV)


# 10 hrs 
t <- 2
image.plot(matrix(extremalFieldsIV[t, 10, ,1], 16, 16))
parameterComb[120+52+128+t,]


# ## Shuffling the dataset:
dim(extremalFieldsIV)
# dim(parameterComb)

# plot(parameterComb[,1], parameterComb[,2], pch=19)
# plot(parameterComb[,2], parameterComb[,3], pch=19)
# plot(parameterComb[,1], parameterComb[,3], pch=19)

# parameterCombCvg  <- parameterComb
extremalFieldsCvgIV  <- extremalFieldsIV # size 52

## Defining test set: save(storeParameterRepTrain, file="ParameterGrid-Frechet-Train.RData")
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/Data/Coverage")
# save(parameterCombCvg,
#      file="PC-cvg-Test-Rep30-Feb22-2025.RData")
save(extremalFieldsCvgIV, file="storeZRep-cvg-Test-Rep30-Feb23-2025IV.RData")
load("storeZRep-cvg-Test-Rep30-Feb23-2025III.RData")
