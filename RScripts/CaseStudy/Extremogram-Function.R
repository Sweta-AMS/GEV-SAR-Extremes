rm(list=ls())
extremogram <- function(loc, 
                        y,
                        threshold,
                        id = NULL,
                        d = NULL,
                        lon.lat = FALSE,
                        dmax = NULL,
                        N = NULL, 
                        breaks = NULL) {
  
  y <- cbind(y)  # Ensure y is a matrix
  
  # Generate pairwise indices if not provided
  if (is.null(id)) {
    n <- nrow(loc)
    is = rep(1:n, n)
    js = rep(1:n, rep(n, n))
    ind <- is > js
    id <- cbind(is, js)[ind, ]
  }
  
  # Compute pairwise distances if not provided
  if (is.null(d)) {
    loc <- as.matrix(loc)
    if (lon.lat) {
      d <- rdist.earth.vec(loc[id[, 1], ], loc[id[, 2], 
      ])
    }
    else {
      d <- rdist.vec(loc[id[, 1], ], loc[id[, 2], ])
    }
  }
  
  colMeans <- apply(y, 
                    2, 
                    mean,
                    na.rm = TRUE)
  
  yCntr = sweep(y,
                2, 
                colMeans)
  
  y1Cntr = yCntr[id[, 1], ]
  y2Cntr = yCntr[id[, 2], ]
  
  # Compute madogram values (corrected formula: 0.5 × mean(|y1 - y2|))
  # vg <- as.integer(y1 > threshold & y2 >threshold), na.rm = TRUE)
  # vg <- (y1 > threshold & y2 > threshold, na.rm = TRUE)
  vg <- rowMeans(cbind((y1Cntr  > threshold & y2Cntr  > threshold)),
                 na.rm = TRUE)

  call <- match.call()
  
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  
  ind <- d <= dmax & !is.na(vg)
  out <- list(d = d[ind], 
              extgram = vg[ind],
              call = call)
  
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("extgram", class(out))
  out
}


# # Example
# threshold <- quantile(annual_precip_subset[, 1], 0.95)  # 95th percentile
# result <- extremogram(loc = coords_land_subset, 
#                       y = annual_precip_subset[, 1], 
#                       threshold = threshold,
#                       lon.lat = TRUE,
#                       N = 50)
# fields::bplot.xy(result$d, result$extremogram, 
#                  breaks = result$centers, 
#                  outline = FALSE)
# bplot.xy( result$d,
#           result$extremogram,
#           #breaks=lookBins$breaks,
#           #xlim=c(0,100),
#           outline=FALSE)
# points(result$centers, result$stats["mean", ], col = "red", pch = 16)
# 
# 

# t <- 3
# values <- (annual_precip_subset[,1] - median(annual_precip_subset[,1]))
# lookBins<- madogram (coords_subset, values, N=50)
# bplot.xy( lookBins$d,
#           lookBins$mgram,
#           breaks=lookBins$breaks,
#           #xlim=c(0,100),
#           outline=FALSE)
# vMean <- lookBins$stats["mean",]
# binCenter<- lookBins$centers
# points( binCenter,vMean, col="red",pch=16)


### Evaluation at the simulated data
## Big Xi=0.8, moderate Kappa2=1, small Tau2=0.0001


## Location
M <- 16
sGrid <- list(x=1:M, y=1:M)
s <- make.surface.grid(sGrid)
plot(s, pch=19)

n <- nrow(s)
print(paste0('Field size: ', n))


Kappa2 <- 0.01
Xi <- 0.9
Tau2 <- 0.00001


# ### -- Generate Test Set for comparison study with likelihood estimation --:
## Required script
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")
m <- 30 # no of replication in the spatial field
R <- 1
n_sims <-1
extremalFields <- array(NA,
                        dim=c(n_sims, R, n, m))
dim(extremalFields)
set.seed(123)


tic()
for (i in 1:n_sims)
{
  cat('GEV parameter loop no', i, '\n')
  
  # Defining GEV parameters:
  shape_val <- Xi
  a_wght <- 4+Kappa2
  
  # Setup LKrig
  LKinfo <- LKrigSetup(s,
                       a.wght=a_wght,
                       nlevel=1,
                       nu=1,
                       NC=40,  # Changed from 8
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
    ySim <- (PHI %*% coefSAR)
    
    # Add measurement error
    scaleParameter <- log(1 + (Tau2))  # Variance of the log-normal is tau2
    locParameter <- (-scaleParameter / 2)
    
    # Log-normal matrix generation, unique for each iteration
    LogNormalMat <- exp(matrix(rnorm(n*m,
                                     mean=locParameter,
                                     sd=sqrt(scaleParameter)),
                               n, m))  # n x m matrix
    
    # Spatial Extremal Fields with Nugget effect
    Z <- ySim*LogNormalMat
    extremalFields[i, j,  , ] <- Z
  }
  toc()
}
toc()

dim(extremalFields)



# Example
i<- 20
threshold <- quantile((extremalFields[1,1, , ]-median(extremalFields[1,1, , ]))/IQR(extremalFields[1,1, , ]), 0.9)  # 95th percentile
result <- extremogram(loc = s, 
                      y = (extremalFields[1,1, , ]-median(extremalFields[1,1, , ]))/IQR(extremalFields[1,1, , ]),
                        # coefSAR,
                        #  
                      
                      dmax = 1000,
                      threshold = threshold,
                      lon.lat = TRUE,
                      N = 50)
fields::bplot.xy(result$centers,
                 result$stats['mean',], 
                 #breaks = result$centers, 
                 outline = FALSE)
abline(h=0, col='magenta3')


# bplot.xy( result$d,
#           result$extremogram,
#           #breaks=lookBins$breaks,
#           xlim=c(0,500),
#           outline=FALSE)
points(result$centers, result$stats["mean", ], col = "red", pch = 16)


image.plot(matrix(extremalFields[1,1, ,20], 40, 40))

