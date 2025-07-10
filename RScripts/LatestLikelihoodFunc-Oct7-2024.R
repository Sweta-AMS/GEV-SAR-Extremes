source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Awght2Omega.R")
library(MASS)
log_like <- function(s,
                     NC=NC,
                     NC.buff=NC.buff,
                     theta,
                     Z)
{
  # Defining GEV parameters
  shape_val <- exp(as.numeric(theta[1]))
  a_wght <- omega2Awght(theta[2])

  # lambda_val <- theta[3]


  LKinfo <-  LKrigSetup(s,
                        a.wght=a_wght,
                        nlevel=1,
                        nu=1,
                        NC=NC,
                        NC.buffer=NC.buff,
                        fixedFunction = NULL)

  B <- LKrigSAR(LKinfo,
                Level=1) # no of basis function used
  B <- as.matrix(spind2spam(B))  # Convert to spam format

  # Simulate field
  PHI <- as.matrix(LKrig.basis(s, LKinfo))
  y <- Z

  # Compute residuals
  coef_vals <- ginv((t(PHI)%*%PHI) + 1e-8*diag(1, ncol(PHI))) %*% (t(PHI)%*%y)
  e_vals <- B%*%coef_vals

  # Compute GEV density
  density_vals <- devd(e_vals,
                       loc=1,
                       scale=shape_val,
                       shape=shape_val)

  density_vals[density_vals==0] <- 1e-10

  log_det_B <- determinant(B,
                           logarithm=TRUE)$modulus

  # Negative log-likelihood
  llike <- -sum(log(density_vals)) -log_det_B
  return(llike)
}

#
# 
# log_like <- function(s,
#                      NC=NC,
#                      NC.buff=NC.buff,
#                      theta,
#                      Z){
#   # Parameters
#   shape_val <- exp(as.numeric(theta[1]))
#   a_wght <- omega2Awght(theta[2])
#   lambda_val <- theta[3]
#   
#   LKinfo <-  LKrigSetup(s,
#                         a.wght=a_wght,
#                         nlevel=1,
#                         nu=1,
#                         NC=NC,
#                         NC.buffer=NC.buff,
#                         fixedFunction = NULL)
#   
#   B <- LKrigSAR(LKinfo,
#                   Level=1) # no of basis function used
#   B <- as.matrix(spind2spam(B))  # Convert to spam format
# 
#   
#   ## Basis matrix
#   PHI <- as.matrix(LKrig.basis(s, LKinfo))
#   
#   # Solve PHI %*% c_vals = Z using QR
#   svd_decomp <- svd(PHI)
#   d_inv <- diag(1 / svd_decomp$d)  # Invert nonzero singular values
#   c_vals <- svd_decomp$v %*% d_inv %*% t(svd_decomp$u) %*% Z
#   
#   
#   # c_vals <- as.matrix(ginv(PHI) %*% Z)
#   
#   # Residuals e = B * c
#   e_vals <- B %*% c_vals
#   
#   # Check for invalid residuals (Frechet case requires e > 0)
#   e_vals[e_vals <= 0] <- 1e-08
#   
#   # GEV log-density
#   log_density <- sum(devd(e_vals,
#                           loc=1, 
#                           scale=shape_val,
#                           shape=shape_val, 
#                           log=TRUE))
#   
#   # Log-determinant adjustment
#   log_det_B <- determinant(B,
#                            logarithm=TRUE)$modulus
#   
#   # Negative log-likelihood
#   negative_llike <- - (log_density + as.numeric(log_det_B))
#   
#   return(negative_llike)
# }


