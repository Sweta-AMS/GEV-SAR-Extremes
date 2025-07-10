madogram <- function (loc,
                      y, 
                      id = NULL,
                      d = NULL,
                      lon.lat = FALSE, 
                      dmax = NULL, 
                      N = NULL, 
                      breaks = NULL) 
{
  y <- cbind(y)
  if (is.null(id)) {
    n <- nrow(loc)
    is = rep(1:n, n)
    js = rep(1:n, rep(n, n))
    ind <- is > js
    id <- cbind(is, js)[ind, ]
  }
  
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
  
  colMeans <- apply(y, 2, mean, na.rm = TRUE)
  yCntr = sweep(y, 2, colMeans)
  y1Cntr = yCntr[id[, 1], ]
  y2Cntr = yCntr[id[, 2], ]
  
  vg <- 0.5 * rowMeans(cbind(abs(y1Cntr - y2Cntr)), na.rm = TRUE)
  
  
  call <- match.call()
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  ind <- d <= dmax & !is.na(vg)
  out <- list(d = d[ind], mgram = vg[ind], call = call)
  
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("mgram", class(out))
  out
}
  
  
  
#   
#   function(loc,
#                      y,
#                      id = NULL,
#                      d = NULL,
#                      lon.lat = FALSE,
#                      dmax = NULL,
#                      N = NULL,
#                      breaks = NULL) {
# 
#   y <- cbind(y)  # Ensure y is a matrix
# 
#   # Generate pairwise indices if not provided
#   if (is.null(id)) {
#     n <- nrow(loc)
#     id <- which(upper.tri(diag(n)), arr.ind = TRUE)
#     colnames(id) <- c("i", "j")  # Rename for clarity
#   }
# 
#   # Compute pairwise distances if not provided
#   if (is.null(d)) {
#     loc <- as.matrix(loc)
#     if (lon.lat) {
#       d <- rdist.earth.vec(loc[id[, "i"], , drop = FALSE], loc[id[, "j"], , drop = FALSE])
#     } else {
#       d <- rdist.vec(loc[id[, "i"], , drop = FALSE], loc[id[, "j"], , drop = FALSE])
#     }
#   }
# 
#   colMeans <- apply(y,
#                     2,
#                     median,
#                     na.rm = TRUE)
# 
#   yCntr <- sweep(y, 2, colMeans)
# 
#   # Extract paired values
#   y1 <- yCntr[id[, "i"], , drop = FALSE]
#   y2 <- yCntr[id[, "j"], , drop = FALSE]
# 
#   # Compute madogram values (corrected formula: 0.5 × mean(|y1 - y2|))
#   vg <- 0.5 * rowMeans(abs(y1 - y2), na.rm = TRUE)
# 
#   # Store function call
#   call <- match.call()
# 
#   # Determine maximum distance threshold if not provided
#   if (is.null(dmax)) {
#     dmax <- max(d, na.rm = TRUE)
#   }
# 
#   # Sort by distance
#   od <- order(d)
#   d <- d[od]
#   vg <- vg[od]
# 
#   # Filter values within the distance threshold
#   valid_idx <- d <= dmax & !is.na(vg)
#   out <- list(d = d[valid_idx], mgram = vg[valid_idx], call = call)
# 
#   # Compute binned statistics if needed
#   if (!is.null(breaks) | !is.null(N)) {
#     out <- c(out, stats.bin(d[valid_idx], vg[valid_idx], N = N, breaks = breaks))
#   }
# 
#   # Assign class and return result
#   class(out) <- c("mgram", class(out))
#   return(out)
# }





# 
# ### -- Generate Test Set for comparison study with likelihood estimation --:
# ## Location
# M <- 16
# sGrid <- list(x=1:M, y=1:M)
# s <- make.surface.grid(sGrid)
# plot(s, pch=19)
# 
# n <- nrow(s)
# print(paste0('Field size: ', n))
# 
# 
# Kappa2 <- 0.01
# Xi <- 0.9
# Tau2 <- 0.00001
# 
# ## Required script
# setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
# source("requiredPackages.R")
# source("LKrigSimSAR.R")
# source("LKrigSAREvd.R")
# m <- 30 # no of replication in the spatial field
# R <- 1
# n_sims <-1
# extremalFields <- array(NA,
#                         dim=c(1, R, n, m))
# dim(extremalFields)
# set.seed(123)
# 
# 
# tic()
# for (i in 1:n_sims)
# {
#   cat('GEV parameter loop no', i, '\n')
#   
#   # Defining GEV parameters:
#   shape_val <- Xi
#   a_wght <- 4+Kappa2
#   
#   # Setup LKrig
#   LKinfo <- LKrigSetup(s,
#                        a.wght=a_wght,
#                        nlevel=1,
#                        nu=1,
#                        NC=40,  # Changed from 8
#                        NC.buffer=4,
#                        fixedFunction=NULL)
#   tic()
#   for(j in 1:R)
#   {
#     cat('For rep', j, '\n')
#     # Simulate SAR coefficients
#     coefSARList <- LKrigSAREvd(LKinfo,
#                                loc=1,
#                                scale=shape_val,
#                                shape=shape_val,
#                                M=m,
#                                asList=FALSE)
#     
#     coefSAR <- coefSARList$coefSAR
#     
#     # Simulate field
#     PHI <- LKrig.basis(s, LKinfo)
#     ySim <- (PHI %*% coefSAR)
#     
#     # Add measurement error
#     scaleParameter <- log(1 + (Tau2))  # Variance of the log-normal is tau2
#     locParameter <- (-scaleParameter / 2)
#     
#     # Log-normal matrix generation, unique for each iteration
#     LogNormalMat <- exp(matrix(rnorm(n*m,
#                                      mean=locParameter,
#                                      sd=sqrt(scaleParameter)),
#                                n, m))  # n x m matrix
#     
#     # Spatial Extremal Fields with Nugget effect
#     Z <- ySim*LogNormalMat
#     extremalFields[i, j,  , ] <- Z
#   }
#   toc()
# }
# toc()
# dim(extremalFields)
# 
# 
# values = (extremalFields[1,1, , ]-median(extremalFields[1,1, , ]))/IQR(extremalFields[1,1, , ])
# #values <- (annual_precip_subset[,1] - median(annual_precip_subset[,1]))
# # coords_subset
# lookBins<- madogram(s, 
#                     values,
#                     breaks = seq(0,16, length.out=16))
# bplot.xy( lookBins$d,
#           lookBins$mgram,
#           breaks=lookBins$breaks,
#           #xlim=c(0,100),
#           outline=FALSE)
# vMean <- lookBins$stats["mean",]
# binCenter<- lookBins$centers
# points( binCenter,vMean, col="red",pch=16)
# 
# 
# # Define bins
# bins <- seq(0, 16, length.out = 16)  # Define distance bins
# 
# # Initialize madogram storage
# madogram_vals <- numeric(length(bins) - 1)  # Store results
# 
# # Compute pairwise distances
# library(fields)  # For rdist
# dist_matrix <- rdist(s)  # Compute pairwise distances between locations
# 
# # Compute madogram manually
# for (b in 1:(length(bins) - 1)) {
#   
#   # Find all location pairs within the bin range
#   bin_indices <- which(dist_matrix >= bins[b] & dist_matrix < bins[b + 1], arr.ind = TRUE)
#   
#   # Extract corresponding values
#   diffs <- abs(values[bin_indices[,1], ] - values[bin_indices[,2], ])
#   
#   # Compute madogram (mean of absolute differences)
#   madogram_vals[b] <- mean(diffs, na.rm = TRUE) / 2  
# }
# 
# # Print results
# madogram_vals
# points(bins[1:15], madogram_vals)
