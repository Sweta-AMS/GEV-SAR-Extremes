rm(list = ls())

library(fields)
library(Matrix)

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")

time_SAR_solve <- function(M, kappa2 = 0.01) {
  
  # --------------------------------------------------
  # Spatial grid
  # --------------------------------------------------
  sGrid <- list(x = 1:M, y = 1:M)
  s <- make.surface.grid(sGrid)
  n <- nrow(s)
  
  # --------------------------------------------------
  # LKrig setup
  # --------------------------------------------------
  LKinfo <- LKrigSetup(
    s,
    a.wght = 4 + kappa2,
    nlevel = 1,
    nu = 1,
    NC = M,
    NC.buffer = 0,
    fixedFunction = NULL
  )
  
  # --------------------------------------------------
  # Build SAR matrix
  # --------------------------------------------------
  t_build <- system.time({
    B_sparse <- LKrigSAR(LKinfo, Level = 1)
    B_sparse <- spind2spam(B_sparse)
  })["elapsed"]
  
  # --------------------------------------------------
  # Right-hand side
  # --------------------------------------------------
  e <- rnorm(n)
  
  # --------------------------------------------------
  # Sparse solve
  # --------------------------------------------------
  t_sparse <- system.time({
    c_sparse <- solve(B_sparse, e)
  })["elapsed"]
  
  # --------------------------------------------------
  # Dense solve (for comparison only)
  # --------------------------------------------------
  B_dense <- as.matrix(B_sparse)
  
  t_dense <- system.time({
    c_dense <- solve(B_dense, e)
  })["elapsed"]
  
  # --------------------------------------------------
  # Output
  # --------------------------------------------------
  data.frame(
    Grid = paste0(M, "x", M),
    n = n,
    BuildTime = round(t_build, 3),
    SparseSolve = round(t_sparse, 3),
    DenseSolve = round(t_dense, 3)
  )
}

# --------------------------------------------------
# Run experiment
# --------------------------------------------------
grid_sizes <- c(4, 16, 24, 32)

timing_results <- do.call(
  rbind,
  lapply(grid_sizes, time_SAR_solve)
)

print(timing_results)


## Simulation time:

# ---------------------------------------------------------
# Function: build SAR system (ONE-TIME COST)
# ---------------------------------------------------------
build_SAR_system <- function(M, kappa2 = 0.01) {
  
  # Spatial grid
  sGrid <- list(x = 1:M, y = 1:M)
  s <- make.surface.grid(sGrid)
  n <- nrow(s)
  
  # Setup LKrig
  t_build <- system.time({
    LKinfo <- LKrigSetup(
      s,
      a.wght = 4 + kappa2,
      nlevel = 1,
      nu = 1,
      NC = M,
      NC.buffer = 0,
      fixedFunction = NULL
    )
    
    # SAR matrix (sparse)
    B <- LKrigSAR(LKinfo, Level = 1)
    B <- spind2spam(B)
    
    # Basis matrix
    PHI <- LKrig.basis(s, LKinfo)
  })["elapsed"]
  
  list(
    s = s,
    n = n,
    B = B,
    PHI = PHI,
    build_time = t_build
  )
}

# ---------------------------------------------------------
# Function: simulate ONE spatial field
# ---------------------------------------------------------
simulate_field <- function(B, PHI, xi) {
  
  n <- nrow(B)
  
  # Step 1: innovations
  e <- rgev(n, loc = 1, scale = xi, shape = xi)
  
  # Step 2: SAR solve
  c <- solve(B, e)
  
  # Step 3: map to spatial field
  Z <- as.vector(PHI %*% c)
  
  return(Z)
}

# ---------------------------------------------------------
# Benchmark simulation time
# ---------------------------------------------------------
benchmark_simulation <- function(M, xi = 0.5, R = 100) {
  
  SAR <- build_SAR_system(M)
  
  t_sim <- system.time({
    for (r in 1:R) {
      Z <- simulate_field(SAR$B, SAR$PHI, xi)
    }
  })["elapsed"]
  
  data.frame(
    Grid = paste0(M, "x", M),
    n = M^2,
    BuildTime_sec = SAR$build_time,
    AvgSimTime_sec = t_sim / R
  )
}

# ---------------------------------------------------------
# Run benchmark
# ---------------------------------------------------------
grid_sizes <- c(8, 16, 24, 32)

results <- do.call(
  rbind,
  lapply(grid_sizes, benchmark_simulation)
)

print(results)

