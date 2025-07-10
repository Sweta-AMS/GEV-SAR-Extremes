# Function to generate spatial extremal fields using LKrig and GEV parameters
#
# Args:
#   s: Matrix of spatial locations (n x 2 for lon/lat)
#   m: Number of ensemble members (climate model realizations)
#   parameters: Matrix of GEV parameters (n_sims x 3), where:
#       - Column 1: Shape parameter (ξ)
#       - Column 2: Kappa^2 parameter
#       - Column 3: Lambda^2 parameter
#
# Returns:
# extremalFields: 3D array of simulated fields (n_sims x n x m)

generateSpatialExtFields <- function(s, m, parameters) {

  n <- nrow(s)  # Number of spatial locations
  n_sims <- nrow(parameters)  # Number of parameter sets
  
  # Pre-allocate the output array
  extremalFields <- array(NA, dim = c(n_sims, n, m))
  
  tic()  # Start timing
  # set.seed(123)  # Ensure reproducibility
  
  for (i in 1:n_sims) {
    # Extract parameters for the current simulation
    shape_val <- parameters[i, 1]
    a_wght <- parameters[i, 2]
    lambda_val <- parameters[i, 3]
    
    # Setup LKrig
    LKinfo <- LKrigSetup(s,
                         a.wght = a_wght,
                         nlevel = 1,
                         nu = 1,
                         NC = 16,
                         NC.buffer = 4,
                         fixedFunction = NULL)
    
    # Compute Log-Normal parameters
    scaleParameter <- log(1 + lambda_val)  
    locParameter <- -scaleParameter / 2  
    
    # Generate SAR coefficients
    coefSARList <- LKrigSAREvd(LKinfo,
                               loc = 1,
                               scale = shape_val,
                               shape = shape_val,
                               M = m,
                               asList = FALSE)
    
    coefSAR <- coefSARList$coefSAR
    
    # Compute basis matrix
    PHI <- LKrig.basis(s, LKinfo)
    ySim <- PHI %*% coefSAR  # Simulated spatial field
    
    # Generate Log-Normal noise
    LogNormalMat <- exp(matrix(rnorm(n * m,
                                     mean=locParameter,
                                     sd=sqrt(scaleParameter)),
                               n, m))
    
    # Store results
    extremalFields[i, , ] <- ySim * LogNormalMat
  }
  toc()  # End timing
  return(extremalFields)  # Explicitly return results
}

