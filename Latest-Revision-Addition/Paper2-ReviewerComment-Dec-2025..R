## ================================================================
rm(list = ls())
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")

source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")

library(fields)
library(evir)

library(evd)

## ------------------------------------------------------------------
## 1. Grid setup
## ------------------------------------------------------------------
M <- 16
sGrid <- list(x = 1:M, y = 1:M)
s <- make.surface.grid(sGrid)
n <- nrow(s)
set.panel(1,1)
plot(s)

sx <- s[,1]
sy <- s[,2]

center_idx <- which(
  sx %in% c(M/2, M/2 + 1) &
    sy %in% c(M/2, M/2 + 1)
)
center_loc <- center_idx[1]
cat("Center location:", center_loc, "\n")

## ------------------------------------------------------------------
## 2. Parameter sets
## ------------------------------------------------------------------
a_wght_vec <- c(4.001, 4.01, 5)

xi_vec     <- c(0.2, 0.5, 0.8)

lambda_val     <- 0.01
scaleParameter <- log(1 + lambda_val)
locParameter   <- -scaleParameter / 2

m <- 5000

## ------------------------------------------------------------------
## 3. Storage object
## ------------------------------------------------------------------
results <- list()
set.seed(123)

## ------------------------------------------------------------------
## 4. Loop over spatial dependence (a.wght = κ²)
## ------------------------------------------------------------------
for (a_val in a_wght_vec) {
  
  cat("\n==============================\n")
  cat("Simulating for a.wght =", a_val, "\n")
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
  
  results[[as.character(a_val)]] <- list(
    tail_summary = data.frame(
      xi_true    = xi_vec,
      hill_alpha = NA_real_,
      hill_xi    = NA_real_,
      fgev_xi    = NA_real_
    ),
    marginals = list(),
    q_scaled  = list(),
    hill_obj  = list(),
    fields    = list(),
    chi       = list()
  )
  
  for (xi_val in xi_vec) {
    
    cat("  xi =", xi_val, "\n")
    
    coefSARList <- LKrigSAREvd(
      LKinfo,
      loc   = 1,
      scale = xi_val,
      shape = xi_val,
      M     = m,
      asList = FALSE
    )
    coefSAR <- coefSARList$coefSAR
    ySim <- PHI %*% coefSAR
    
    LogNormalMat <- exp(matrix(
      rnorm(n*m, mean = locParameter, sd = sqrt(scaleParameter)),
      n, m
    ))
    
    fields <- ySim * LogNormalMat
    results[[as.character(a_val)]]$fields[[as.character(xi_val)]] <- fields
    
    marginal <- fields[center_loc,]
    results[[as.character(a_val)]]$marginals[[as.character(xi_val)]] <- marginal
    
    h <- hill(marginal)
    results[[as.character(a_val)]]$hill_obj[[as.character(xi_val)]] <- h
    
    idx <- floor(0.9 * length(h$y)) : length(h$y)
    alpha_hat <- median(h$y[idx])
    xi_hat_hill <- 1 / alpha_hat
    
    fit <- fgev(marginal)
    xi_hat_gev <- fit$estimate["shape"]
    
    results[[as.character(a_val)]]$tail_summary[
      results[[as.character(a_val)]]$tail_summary$xi_true == xi_val,
      c("hill_alpha","hill_xi","fgev_xi")
    ] <- c(alpha_hat, xi_hat_hill, xi_hat_gev)
    
    pvals <- ppoints(length(marginal))
    qfrechet <- function(p, alpha) (-log(p))^(-1/alpha)
    alpha_true <- 1 / xi_val
    
    q_theory <- qfrechet(pvals, alpha_true)
    scale_fac <- quantile(marginal,0.95) / quantile(q_theory,0.95)
    
    q_scaled <- scale_fac * q_theory
    results[[as.character(a_val)]]$q_scaled[[as.character(xi_val)]] <- q_scaled
  }
}

# ## ================================================================
# ## 6. PLOTTING TO PDF
# ## ================================================================
# 
# while (dev.cur() > 1) dev.off()
## ================================================================
## ================================================================
setwd("~/Desktop")
pdf("Hill_panels_kappa_xi.pdf", width = 10, height = 10)
par(
  mfrow = c(length(a_wght_vec), length(xi_vec)),
  mar  = c(2.2, 3.2, 2.2, 1),
  oma  = c(3.5, 4.5, 4, 1)
)
axis_cex  <- 1.8
title_cex <- 1.8

# ---- Compute ymax per column ----
col_ymax <- numeric(length(xi_vec))
for (j in seq_along(xi_vec)) {
  xi_val <- xi_vec[j]
  all_alpha <- c()
  
  for (a_val in a_wght_vec) {
    h_obj <- results[[as.character(a_val)]]$hill_obj[[as.character(xi_val)]]
    all_alpha <- c(all_alpha, h_obj$y)
  }
  
  col_ymax[j] <- max(all_alpha, na.rm = TRUE)
}

# ---- Panel plotting ----
for (i in seq_along(a_wght_vec)) {
  
  a_val <- a_wght_vec[i]
  
  for (j in seq_along(xi_vec)) {
    
    xi_val <- xi_vec[j]
    
    h_obj <- results[[as.character(a_val)]]$hill_obj[[as.character(xi_val)]]
    est_alpha <- results[[as.character(a_val)]]$tail_summary$hill_alpha[
      results[[as.character(a_val)]]$tail_summary$xi_true == xi_val
    ]
    true_alpha <- 1 / xi_val
    
    x_perc <- (1:length(h_obj$y)) / length(h_obj$y) * 100
    is_bottom_row <- (i == length(a_wght_vec))
    
    plot(
      x_perc, h_obj$y,
      type = "l",
      lwd  = 1.8,
      col  = "black",
      main = bquote(kappa^2 == .(a_val-4) ~ "," ~ xi == .(xi_val)),
      xlab = "", ylab = "",
      xaxt = ifelse(is_bottom_row, "s", "n"),
      
      ## >>>>>>> FIXED LOWER Y LIMIT = 1.5 <<<<<<<<
      ylim = c(1.2, col_ymax[j]),
      
      cex.main = title_cex,
      cex.axis = axis_cex
    )
    
    abline(h = true_alpha, col="red",  lwd=3, lty=2)
    abline(h = est_alpha,  col="blue", lwd=3, lty=3)
    
    legend(
      "topleft",
      legend = c(
        bquote(alpha[true] == .(round(true_alpha, 3))),
        bquote(alpha[est]  == .(round(est_alpha, 3)))
      ),
      col = c("red","blue"),
      lty = c(2,3), lwd = 2,
      bty = "n", cex = 1.4
    )
    
    if (j == 1) {
      mtext(expression(alpha[k]), side = 2, line = 3.5, cex = 1.2)
    }
  }
}

mtext("Upper order statistics (in %)",
      side = 1, outer = TRUE, line = 1.8, cex = 1.2)

mtext(
  expression("Hill Plots for GEV-SAR for varying " * kappa^2 * ", " * xi),
  side = 3, outer = TRUE, line = 1, cex = 1.8, font = 2
)
dev.off()


#############################################################
### 6B. QQ PLOTS (MARGINAL DISTRIBUTION CHECK)
#############################################################
setwd("~/Downloads")
pdf("QQ_panels_kappa_xi.pdf", width = 12, height = 12)
par(
  mfrow = c(length(a_wght_vec), length(xi_vec)),
  mar  = c(3, 3, 3, 1),
  oma  = c(5, 5, 4, 1)
)
axis_cex  <- 1.8   # larger numeric ticks
title_cex <- 2.0
label_cex <- 2.0
for (i in seq_along(a_wght_vec)) {
  a_val <- a_wght_vec[i]
  
  for (j in seq_along(xi_vec)) {
    xi_val <- xi_vec[j]
    
    # empirical marginals
    marg <- results[[as.character(a_val)]]$marginals[[as.character(xi_val)]]
    marg_sorted <- sort(marg)
    
    # theoretical Frechet quantiles (already scaled by you)
    q_theory <- results[[as.character(a_val)]]$q_scaled[[as.character(xi_val)]]
    q_sorted <- sort(q_theory)
    
    # Compute robust fit line slope + intercept (25–75% quantiles)
    emp_q <- quantile(marg_sorted, probs = c(0.25, 0.75))
    th_q  <- quantile(q_sorted,   probs = c(0.25, 0.75))
    slope <- diff(emp_q) / diff(th_q)
    intercept <- emp_q[1] - slope * th_q[1]
    
    # ---- QQ Plot ----
    plot(
      q_sorted, marg_sorted,
      pch = 19,
      col = rgb(0, 0, 0, 0.35),
      cex = 0.8,
      main = bquote(kappa^2 == .(a_val-4) ~ "," ~ xi == .(xi_val)),
      xlab = "", ylab = "",
      cex.main = title_cex,
      cex.axis = axis_cex
    )
    
    # reference robust line
    abline(intercept, slope, col = "red", lwd = 2.5, lty = 2)
    
    # perfect 1:1 line (baseline)
    abline(0, 1, col = "blue", lwd = 2, lty = 3)
  }
}

# GLOBAL AXIS LABELS
mtext("Theoretical Scaled Frechét Quantiles", 
      side = 1, outer = TRUE, line = 2, cex = label_cex)
mtext("Empirical Marginal Quantiles", side = 2,
      outer = TRUE, line = 2, cex = label_cex)

# GLOBAL TITLE
mtext(
  expression("QQ Plots for Center-Point Marginals under varying " * kappa^2 * ", " * xi),
  side = 3, outer = TRUE, line = 0.5, cex = 2, font = 2
)
dev.off()



################################################################################
## ================================================================
## ================================================================
##  Chi(h) panels — clean layout, no overlap, journal-ready
## ================================================================
compute_chi_fixed_center <- function(fields,
                                     coords, 
                                     center_loc,
                                     nbins=10,
                                     q=0.95){
  
  # Distance matrix
  D <- rdist(coords)
  d <- D[center_loc, ]
  d[center_loc] <- NA
  
  # Center exceedances
  Yc <- fields[center_loc, ]
  u_c <- quantile(Yc, q)
  Ic  <- (Yc > u_c)
  
  p_c <- mean(Ic)
  if (p_c == 0)
    stop("No exceedances at center — lower threshold q.")
  
  chi_j <- rep(NA_real_, length(d))
  
  for (j in seq_along(d)) {
    if (j == center_loc) next
    
    Yj <- fields[j, ]
    u_j <- quantile(Yj, q)
    Ij  <- (Yj > u_j)
    
    # χ_j = P(Y_j > u_j | Y_c > u_c)
    chi_j[j] <- mean(Ij & Ic)/p_c
  }
  
  # Bin by distance
  edges <- seq(0, max(d, na.rm = TRUE), length.out = nbins + 1)
  mids  <- (edges[-1] + edges[-length(edges)]) / 2
  
  chi_h <- numeric(nbins)
  for (b in seq_len(nbins)) {
    idx <- which(d > edges[b] & d <= edges[b + 1])
    chi_h[b] <- if (length(idx) > 0)
      mean(chi_j[idx], na.rm = TRUE) else NA
  }
  list(distance = mids, chi = chi_h)
}
###################################################################
# Add storage lists to results
for (a_val in a_wght_vec) {
  results[[as.character(a_val)]]$chi_95 <- list()
  results[[as.character(a_val)]]$chi_97 <- list()
  results[[as.character(a_val)]]$chi_99 <- list()
}

# Compute χ for each (κ², ξ, q)
for (a_val in a_wght_vec) {
  for (xi_val in xi_vec) {

    fields_dep <- results[[as.character(a_val)]]$fields[[as.character(xi_val)]]

    # ---- q = 0.95 -----------------------------------------
    results[[as.character(a_val)]]$chi_95[[as.character(xi_val)]] <-
      compute_chi_fixed_center(fields_dep, 
                               s, 
                               center_loc,
                               nbins=10,
                               q=0.95)

    # ---- q = 0.97 -----------------------------------------
    results[[as.character(a_val)]]$chi_97[[as.character(xi_val)]] <-
      compute_chi_fixed_center(fields_dep, 
                               s, 
                               center_loc,
                               nbins=10,
                               q=0.97)

    # ---- q = 0.99 -----------------------------------------
    results[[as.character(a_val)]]$chi_99[[as.character(xi_val)]] <-
      compute_chi_fixed_center(fields_dep,
                               s, 
                               center_loc,
                               nbins=10, 
                               q=0.99)
  }
}
cat("\n Finished computing χ(h) for q = 0.95, 0.97, 0.99\n")
###############################################################################




# 0) Close any open devices safely
while (dev.cur() > 1) dev.off()

# 1) Output file (robust PDF device)
out_file <- file.path(path.expand("~/Desktop"),
                      "Chi_panels_kappa_xi_multiple_u_latest_median.pdf")

grDevices::cairo_pdf(filename = out_file, width = 11, height = 13)


# 2) Panel layout — IMPORTANT FIXES HERE
par(
  mfrow = c(length(a_wght_vec), length(xi_vec)),
  mar  = c(2.4, 3.0, 3.0, 0.8),   # ↑ more space ABOVE each panel
  oma  = c(3.6, 4.4, 5.2, 0.6)    # ↑ more space for global labels + title
)

# 3) Font sizes
axis_cex  <- 1.8   # numeric tick labels
title_cex <- 1.6   # panel titles (smaller to avoid overlap)
label_cex <- 1.8
line_wd   <- 2.4

cols <- c("black", "blue", "red")

# 4) Loop over panels
for (i in seq_along(a_wght_vec)) {
  a_val <- a_wght_vec[i]
  
  for (j in seq_along(xi_vec)) {
    xi_val <- xi_vec[j]
    
    chi95 <- results[[as.character(a_val)]]$chi_95[[as.character(xi_val)]]
    chi97 <- results[[as.character(a_val)]]$chi_97[[as.character(xi_val)]]
    chi99 <- results[[as.character(a_val)]]$chi_99[[as.character(xi_val)]]
    
    ymax <- max(c(chi95$chi, chi97$chi, chi99$chi), na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    is_bottom <- (i == length(a_wght_vec))
    
    plot(
      chi95$distance, chi95$chi,
      type = "b",
      pch  = 16,
      col  = cols[1],
      lwd  = line_wd,
      ylim = c(0, 1.05 * ymax),
      xlab = "",
      ylab = "",
      xaxt = ifelse(is_bottom, "s", "n"),
      main = bquote(kappa^2 == .(a_val-4) ~ "," ~ xi == .(xi_val)),
      cex.axis = axis_cex,
      cex.main = 2
    )
    
    lines(chi97$distance, chi97$chi,
          type = "b", pch = 16, col = cols[2], lwd = line_wd)
    
    lines(chi99$distance, chi99$chi,
          type = "b", pch = 16, col = cols[3], lwd = line_wd)
    
    grid(lty = 2, col = "grey80")
    
    # y-axis label only on first column
    if (j == 1) {
      mtext(expression(chi(h)), side = 2, line = 2.4, cex = label_cex)
    }
  }
}

# 5) ONE common x-axis label
mtext("Distance from center point",
      side = 1, outer = TRUE, line = 1.9, cex = label_cex)

# 6) ONE clean global title — moved UP so it never overlaps
mtext(
  expression("Asymptotic Dependence " * chi(h) *
               " for varying " * kappa^2 * ", " * xi * ", and thresholds"),
  side = 3, outer = TRUE, line = 1.6, cex = 2.0
)

# 7) Close device
dev.off()

cat("✓ Saved figure to:\n", out_file, "\n")


################################################################################
while (dev.cur() > 1) dev.off()
# Output
out_file <- "~/Desktop/Chi_panels_kappa_xi_multiple_u_latest_median.pdf"
grDevices::cairo_pdf(out_file, width = 11, height = 13)

# Layout
par(
  mfrow = c(length(a_wght_vec), length(xi_vec)),
  mar  = c(2.2, 2.2, 3.0, 0.6),   # inner margins
  oma  = c(4.0, 4.6, 5.0, 0.6)    # outer margins for common labels
)

axis_cex  <- 1.7
title_cex <- 1.4
label_cex <- 1.7
line_wd   <- 2.4

cols <- c("black", "blue", "red")

for (i in seq_along(a_wght_vec)) {
  a_val <- a_wght_vec[i]
  
  for (j in seq_along(xi_vec)) {
    xi_val <- xi_vec[j]
    
    chi95 <- results[[as.character(a_val)]]$chi_95[[as.character(xi_val)]]
    chi97 <- results[[as.character(a_val)]]$chi_97[[as.character(xi_val)]]
    chi99 <- results[[as.character(a_val)]]$chi_99[[as.character(xi_val)]]
    
    ymax <- max(c(chi95$chi, chi97$chi, chi99$chi), na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    
    is_bottom <- (i == length(a_wght_vec))
    
    # plot(
    #   chi95$distance, chi95$chi,
    #   type = "b",
    #   pch  = 16,
    #   col  = cols[1],
    #   lwd  = line_wd,
    #   ylim = c(0, 1.05 * ymax),
    #   xlab = "",
    #   ylab = "",
    #   xaxt = ifelse(is_bottom, "s", "n"),
    #   main = bquote(kappa^2 == .(a_val) ~ "," ~ xi == .(xi_val)),
    #   cex.axis = axis_cex,
    #   cex.main = 2.3
    # )
    
    plot(
      chi95$distance, chi95$chi,
      type = "b",
      pch  = 16,
      col  = cols[1],
      lwd  = line_wd,
      ylim = c(0, 1.05 * ymax),
      xlab = "",
      ylab = "",
      xaxt = ifelse(is_bottom, "s", "n"),
      cex.axis = axis_cex
    )
    
    ## ---- ADD TITLE SEPARATELY ----
    title(
      main = bquote(kappa^2 == .(a_val-4) ~ "," ~ xi == .(xi_val)),
      cex.main = 2.3,   # title size
      line = 1.2        # <<< controls vertical spacing (key!)
    )
    
    lines(chi97$distance, chi97$chi, type = "b",
          pch = 16, col = cols[2], lwd = line_wd)
    
    lines(chi99$distance, chi99$chi, type = "b",
          pch = 16, col = cols[3], lwd = line_wd)
    
    grid(lty = 2, col = "grey80")
  }
}

## ---- ONE COMMON AXIS LABELS ----
mtext("Distance from center point",
      side = 1, outer = TRUE, line = 1.8, cex = label_cex)

mtext(expression(chi(h)),
      side = 2, outer = TRUE, line = 1.6, cex = label_cex)

## ---- ONE COMMON TITLE ----
mtext(
  expression("Asymptotic Dependence " * chi(h) *
               " for varying " * kappa^2 * ", " * xi * ", and thresholds"),
  side = 3, outer = TRUE, line = 1.6, cex = 2.0
)

dev.off()

cat("✓ Saved with ONE common x- and y-axis label:\n", out_file, "\n")



setwd("~/Downloads")
pdf("QQ_panels_kappa_xi.pdf", width = 12, height = 12)
par(mfrow=c(length(a_wght_vec), length(xi_vec)),
   mar=c(3, 3, 3, 1),
   oma=c(5, 5, 4, 1))
axis_cex  <- 1.8
title_cex <- 2.0
label_cex <- 2.0
for (i in seq_along(a_wght_vec)) {
  a_val <- a_wght_vec[i]
  
  for (j in seq_along(xi_vec)) {
    xi_val <- xi_vec[j]
    
    ## ---- Empirical marginals ----
    marg <- results[[as.character(a_val)]]$marginals[[as.character(xi_val)]]
    marg_sorted <- sort(marg)
    
    ## ---- Theoretical Frechét quantiles ----
    q_theory <- results[[as.character(a_val)]]$q_scaled[[as.character(xi_val)]]
    q_sorted <- sort(q_theory)
    
    ## ---- Trim extreme tails (1%–99%) ----
    n <- length(marg_sorted)
    idx_keep <- seq(ceiling(0.0001 * n), floor(0.9995 * n))
    
    marg_trim <- marg_sorted[idx_keep]
    q_trim    <- q_sorted[idx_keep]
    
    ## ---- Robust fit (25–75%) ----
    emp_q <- quantile(marg_trim, probs = c(0.05, 0.95))
    th_q  <- quantile(q_trim,   probs = c(0.05, 0.95))
    slope <- diff(emp_q) / diff(th_q)
    intercept <- emp_q[1] - slope * th_q[1]
    
    ## ---- QQ Plot ----
    plot(
      q_trim, marg_trim,
      pch = 19,
      col = rgb(0, 0, 0, 0.35),
      cex = 0.8,
      main = bquote(kappa^2 == .(a_val-4) ~ "," ~ xi == .(xi_val)),
      xlab = "", ylab = "",
      cex.main = title_cex,
      cex.axis = axis_cex
    )
    
    ## ---- Reference lines ----
    abline(intercept, slope, col = "magenta3",  lwd = 2.5, lty = 2)
    abline(0, 1, col = "blue", lwd = 2,   lty = 3)
  }
}

## ---- Global labels ----
mtext("Theoretical Scaled Fréchet Quantiles", 
      side = 1, outer = TRUE, line = 2, cex = label_cex)

mtext("Empirical Marginal Quantiles", 
      side = 2, outer = TRUE, line = 2, cex = label_cex)

mtext(
  expression("QQ Plots for Center-Point Marginals under varying " * kappa^2 * ", " * xi),
  side = 3, outer = TRUE, line = 0.5, cex = 2, font = 2
)
dev.off()

