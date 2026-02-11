## ================================================================
## Hill plots for Gaussian-SAR and Pareto-SAR
## Varying a.wght (kappa^2) ONLY
## ================================================================
rm(list = ls())
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")

source("requiredPackages.R")
source("LKrigSimSAR.R")

library(fields)
library(evir)

## ------------------------------------------------------------------
## 1. Grid setup
## ------------------------------------------------------------------
M <- 16
sGrid <- list(x = 1:M, y = 1:M)
s <- make.surface.grid(sGrid)
n <- nrow(s)

sx <- s[,1]; sy <- s[,2]
center_loc <- which(
  sx %in% c(M/2, M/2 + 1) &
    sy %in% c(M/2, M/2 + 1)
)[1]

cat("Center location:", center_loc, "\n")

## ------------------------------------------------------------------
## 2. Parameters (ONLY a.wght varies)
## ------------------------------------------------------------------
a_wght_vec <- c(4.001, 4.01, 5)   # kappa^2 = a.wght - 4
m <- 5000
set.seed(123)

## Lognormal scaling (same as GEV-SAR)
lambda_val     <- 0.01
scaleParameter <- log(1 + lambda_val)
locParameter   <- -scaleParameter / 2

## ------------------------------------------------------------------
## 3. SAR generators (nlevel = 1)
## ------------------------------------------------------------------

LKrigSARGaussian1 <- function(LKinfo, M) {
  B <- LKrigSAR(LKinfo, Level = 1)
  B <- spind2spam(B)
  mLevel <- LKinfo$latticeInfo$mLevel[1]
  E <- matrix(rnorm(mLevel * M), mLevel, M)
  as.matrix(solve(B, E))
}

LKrigSARPareto1 <- function(LKinfo, M) {
  B <- LKrigSAR(LKinfo, Level = 1)
  B <- spind2spam(B)
  mLevel <- LKinfo$latticeInfo$mLevel[1]
  U <- runif(mLevel * M)
  E <- matrix((1 - U)^(-1), mLevel, M)  # Pareto(alpha=1)
  as.matrix(solve(B, E))
}

## ------------------------------------------------------------------
## 4. Storage
## ------------------------------------------------------------------
results_gaussian <- list()
results_pareto   <- list()

## ------------------------------------------------------------------
## 5. Simulation loop (ONLY a.wght)
## ------------------------------------------------------------------
for (a_val in a_wght_vec) {
  
  cat("\n==============================\n")
  cat("a.wght =", a_val, "\n")
  cat("==============================\n")
  
  LKinfo <- LKrigSetup(
    s,
    a.wght = a_val,
    nlevel = 1,
    nu = 1,
    NC = 16,
    NC.buffer = 4
  )
  
  PHI <- LKrig.basis(s, LKinfo)
  
  ## ---- Gaussian SAR ----
  coefG <- LKrigSARGaussian1(LKinfo, M = m)
  yG <- PHI %*% coefG
  
  ## ---- Pareto SAR ----
  coefP <- LKrigSARPareto1(LKinfo, M = m)
  yP <- PHI %*% coefP
  
  ## ---- SAME lognormal scaling ----
  LogNormalMat <- exp(matrix(
    rnorm(n * m, mean = locParameter, sd = sqrt(scaleParameter)),
    n, m
  ))
  
  yG <- yG * LogNormalMat
  yP <- yP * LogNormalMat
  
  ## ---- Center location ----
  margG <- as.numeric(yG[center_loc, ])
  margP <- as.numeric(yP[center_loc, ])
  
  ## ---- Hill estimators ----
  hG <- hill(margG)
  hP <- hill(margP)
  
  ## ---- Hill estimate (upper 10%, SAME as GEV-SAR) ----
  idxG <- floor(0.9 * length(hG$y)) : length(hG$y)
  idxP <- floor(0.9 * length(hP$y)) : length(hP$y)
  
  alpha_hat_G <- median(hG$y[idxG], na.rm = TRUE)
  alpha_hat_P <- median(hP$y[idxP], na.rm = TRUE)
  
  results_gaussian[[as.character(a_val)]] <- list(
    hill_obj  = hG,
    alpha_hat = alpha_hat_G
  )
  
  results_pareto[[as.character(a_val)]] <- list(
    hill_obj  = hP,
    alpha_hat = alpha_hat_P
  )
}

## ------------------------------------------------------------------
## 6. PLOT: Gaussian-SAR Hill plots (WITH blue dotted estimate)
setwd("~/Desktop")
pdf("Hill_Gaussian_SAR_kappa_panels.pdf", width = 10, height = 4)
par(
  mfrow = c(1, length(a_wght_vec)),
  mar  = c(2.2, 3.2, 2.2, 1),
  oma  = c(3.5, 4.5, 3.5, 1)
)
axis_cex  <- 1.8
title_cex <- 1.8
for (a_val in a_wght_vec) {
  
  res <- results_gaussian[[as.character(a_val)]]
  h   <- res$hill_obj
  alpha_hat <- res$alpha_hat
  
  x_perc <- (1:length(h$y)) / length(h$y) * 100
  
  plot(
    x_perc, h$y,
    type = "l",
    col  = "black",
    lwd  = 1.8,
    xlab = "",
    ylab = "",
    main = bquote(kappa^2 == .(a_val - 4)),
    cex.axis = axis_cex,
    cex.main = title_cex
  )
  
  ## Hill estimate
  abline(h = alpha_hat, col = "blue", lwd = 3, lty = 3)
}

mtext("Upper order statistics (in %)",
      side = 1, outer = TRUE, line = 1, cex = 1.2)

mtext(expression(alpha[k]),
      side = 2, outer = TRUE, line = 1, cex = 1.2)

mtext("Hill plots for Gaussian-SAR under varying spatial dependence",
      side = 3, outer = TRUE, line = 1, cex = 1.5, font = 2)

dev.off()


## ------------------------------------------------------------------
## 7. PLOT: Pareto-SAR Hill plots (WITH true + estimate)
## ------------------------------------------------------------------
setwd("~/Desktop")
pdf("Hill_Pareto_SAR_kappa_panels.pdf", width = 10, height = 4)
par(
  mfrow = c(1, length(a_wght_vec)),
  mar  = c(2.2, 3.2, 2.2, 1),
  oma  = c(3.5, 4.5, 3.5, 1)
)

for (a_val in a_wght_vec) {
  
  res <- results_pareto[[as.character(a_val)]]
  h   <- res$hill_obj
  alpha_hat <- res$alpha_hat
  
  x_perc <- (1:length(h$y)) / length(h$y) * 100
  
  plot(
    x_perc, h$y,
    type = "l",
    col  = "black",
    lwd  = 1.8,
    xlab = "",
    ylab = "",
    main = bquote(kappa^2 == .(a_val - 4)),
    cex.axis = axis_cex,
    cex.main = title_cex
  )
  
  ## True alpha = 1
  abline(h = 1, col = "red", lwd = 3, lty = 2)
  
  ## Hill estimate
  abline(h = alpha_hat, col = "blue", lwd = 3, lty = 3)
}

mtext("Upper order statistics (in %)",
      side = 1, outer = TRUE, line = 1, cex = 1.2)

mtext(expression(alpha[k]),
      side = 2, outer = TRUE, line = 1, cex = 1.2)

mtext("Hill plots for Pareto-SAR under varying spatial dependence",
      side = 3, outer = TRUE, line = 1, cex = 1.5, font = 2)

dev.off()
