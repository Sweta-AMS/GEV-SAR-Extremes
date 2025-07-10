rm(list=ls())

## -- Required Package library -- :
source("~/Desktop/Tentative Results-Proj2/R scripts/requiredPackages.R")
library(reticulate)
library(MASS)
library(survival)
library(fitdistrplus)
library(RcppCNPy)

## -- Load true paramater values --:
## Bias-correction on the training set
load("~/Downloads/ParameterConfiguration-Jan14-Train-Rep30-Latest.RData")
dim(storeParameterRep) # 50,000x3

# Expand the file path: Training set
file_path <- path.expand("~/Downloads/trainEstimates-n-30-Latest-Jan15-2025.npy") # Estimates over the train set

####  -- Test Set -- ####
# file_path_test <- path.expand("~/Downloads/test_estimates_NGP-n-30-Latest-Jan15-2025.npy") # Estimates over the test set
# load("~/ParameterConfiguration-Jan14-Test-Rep30-Latest.RData") # test set
####  -----  ####

# Import numpy and load the file
np <- import("numpy")
mat <- np$load(file_path)

# Load the NumPy file
cnnEst <- mat
dim(cnnEst)  # 100x100x3

# Bias correction for shape
shape_pred <- cnnEst[, 1]
kappa2_pred <- exp(cnnEst[, 2])
tau2_pred <- exp(cnnEst[, 3])

## check:
length(shape_pred) # should be 50000

shape_true <- storeParameterRep[ ,1]
kappa2_true <- (storeParameterRep[ ,2]-4)
tau2_true <- (storeParameterRep[ ,3])

## check:
length(shape_true)
plot(c(shape_pred),
     c(shape_true),
     pch=19,
     col=alpha(alpha=0.2,'grey'))
abline(0,1, col='magenta')

library(splines)
fit_shape <- lm(shape_true~ns(shape_pred, df=5))
shape_pred_bc <- predict(fit_shape,
                         newdata=data.frame(shape_pred=shape_pred)) 
plot(shape_true, 
     shape_pred_bc,
     pch=19,
     col=alpha(alpha=0.2,'grey'))
abline(0,1, col='magenta')

plot(kappa2_true, 
     kappa2_pred,
     pch=19,
     col=alpha(alpha=0.2,'grey'))
abline(0,1, col='magenta')

# Fit a natural spline
fit_kappa2_bc <- lm(kappa2_true~ns(kappa2_pred, df=5))
kappa2_pred_bc <- predict(fit_kappa2_bc,
                          data.frame(kappa2_pred=kappa2_pred))
plot(kappa2_true, 
     kappa2_pred_bc,
     col=alpha(alpha=0.2,"gray"),
     pch=19,
     main="Spline Fit")
abline(0, 1, col='magenta')


plot(tau2_true, 
     tau2_pred,
     col=alpha(alpha=0.2,"gray"),
     pch=19,
     main="Spline Fit")
abline(0, 1, col='magenta')

model_tau2_bc <- lm(tau2_true~ns(tau2_pred, df=5))
tau2_pred_bc <- predict(model_tau2_bc,
                        data.frame(tau2_pred=tau2_pred))
plot(tau2_true, 
     tau2_pred_bc,
     col=alpha(alpha=0.2,"gray"),
     pch=19,
     main="Spline Fit")
abline(0, 1, col='magenta')



# # Visualize the spline fit
# plot(kappa2_true, kappa2_pred, col = "gray", pch = 19, main = "Spline Fit")
# lines(sort(kappa2_true),
#       predict(model_spline, data.frame(kappa2_true=sort(kappa2_true))), col = "red")



fit_LM_kappa2 <- lm((kappa2_true)~(shape_pred+ kappa2_pred + kappa2_pred^2 +tau2_pred))
kappa2_pred_bc <- predict(fit_LM_kappa2, newdata = data.frame(shape_pred=shape_pred, 
                                                                    kappa2_pred=kappa2_pred,
                                                                    tau2_pred=tau2_pred))

plot((kappa2_true), 
     kappa2_pred_bc,
     #ylim=c(0,2),
     pch=19,
     col=alpha(alpha=0.2,'grey'))
abline(0,1, col='magenta')


length(shape_true)
# # Fit the regression model
# model <- glm(shape_true ~ shape_pred+logkappa2_pred+logtau2_pred)


## Bias correction:
library('splines')
BCshape <- lm(shape_true~bs(shape_pred, knots = c(0.2,0.6,1)))
bc_shape_pred <- BCshape$fitted.values


BCshape <- lm(shape_true~poly(shape_pred, 2))
bc_shape_pred <- BCshape$fitted.values


parameter_vals <- cbind(shape_pred, logkappa2_pred)
fit_tps__shape <- Tps(parameter_vals,
                      shape_true)


summary(kappa2_true)
summary(tau2_true)

plot(c(bc_shape_pred), c(shape_true))
abline(0,1, col='magenta')



BClogkappa2 <- lm(log(kappa2_true)~logkappa2_pred)
bc_logkappa2_pred <- BClogkappa2$fitted.values

plot(c(bc_logkappa2_pred), c(log(kappa2_true)))
abline(0,1, col='magenta')


BClogtau2 <- lm(log(tau2_true)~logtau2_pred)
bc_logtau2_pred <- BClogtau2$fitted.values

plot(c(logtau2_pred), c(log(tau2_true)))
abline(0,1, col='magenta')

##  Fitting Quantile Regression for CI of the prediction:
## Fit Quantile Regression: 
library(quantreg) 
dataShape <- data.frame(true=storeParameterRep[,1],
                        estimated=shape_pred)

datalnKappa2 <- data.frame(true=log(storeParameterRep[,2]-4),
                           estimated=logkappa2_pred)

datalnTau2 <- data.frame(true=log(storeParameterRep[,3]),
                         estimated=logtau2_pred)
# Fit quantile regression
quantiles <- seq(0.025, 0.975, length.out=2)
# c(0.2, 0.25, 0.5, 0.75, 0.9)
modelsShape <- lapply(quantiles, function(tau) {
  rq(true ~ estimated, data = dataShape, tau = tau)
})

modelslnKappa2 <- lapply(quantiles, function(tau) {
  rq(true ~ estimated, data = datalnKappa2, tau = tau)
})

modelslnTau2 <- lapply(quantiles, function(tau) {
  rq(true ~ estimated, data = datalnTau2, tau = tau)
})

# Summarize the models
fittedShape <- lapply(modelsShape, fitted)
QRShapeEst <- matrix(NA, nrow=2, ncol=10000)
for(i in 1:2){
  QRShapeEst[i, ] <- fittedShape[[i]]
}
dim(QRShapeEst)

# Summarize the models
fittedlnKappa2 <- lapply(modelslnKappa2, fitted)
QRlnKappa2Est <- matrix(NA, nrow=2, ncol=10000)
for(i in 1:2){
  QRlnKappa2Est[i, ] <- fittedlnKappa2[[i]]
}
dim(QRlnKappa2Est)

# Summarize the models
fittedlnTau2 <- lapply(modelslnTau2, fitted)
QRlnTau2Est <- matrix(NA, nrow=2, ncol=10000)
for(i in 1:2){
  QRlnTau2Est[i, ] <- fittedlnTau2[[i]]
}
dim(QRlnTau2Est)


## Compute lower and upper bound:
XiCI_lower_bound <- rep(NA, 10000)
XiCI_upper_bound <- rep(NA, 10000)
for(i in 1:10000){
  lower_bd <- QRShapeEst[1,i]
  # quantile(QRShapeEst[,i], 0.025)
  XiCI_lower_bound[i] <- lower_bd
  
  # quantile(QRShapeEst[,i], 0.975)
  upper_bd <- QRShapeEst[2,i]
  XiCI_upper_bound[i] <- upper_bd
}

## Compute lower and upper bound:
lnKappa2CI_lower_bound <- rep(NA, 10000)
lnKappa2CI_upper_bound <- rep(NA, 10000)
for(i in 1:10000){
  # lower_bd <- quantile(QRlnKappa2Est[,i], 0.025)
  lower_bd <- QRlnKappa2Est[1,i]
  lnKappa2CI_lower_bound[i] <- lower_bd
  
  upper_bd <- QRlnKappa2Est[2,i]
  lnKappa2CI_upper_bound[i] <- upper_bd
}


## Compute lower and upper bound:
lnTau2CI_lower_bound <- rep(NA, 10000)
lnTau2CI_upper_bound <- rep(NA, 10000)
for(i in 1:10000){
  lower_bd <- QRlnTau2Est[1,i]
  # quantile(QRlnTau2Est[,i], 0.025)
  lnTau2CI_lower_bound[i] <- lower_bd
  
  upper_bd <- QRlnTau2Est[2,i]
  lnTau2CI_upper_bound[i] <- upper_bd
}

# Combine data into a data frame
envelope_data_shape <- data.frame(
  shape_pred = shape_pred,
  Response = storeParameterRep[, 1],
  CI_lower =  XiCI_lower_bound,
  CI_upper =  XiCI_upper_bound
)

# ggplot2 envelope plot with mathematical expression
library(ggplot2)


plot_xi <- ggplot(envelope_data_shape, aes(x = shape_pred, y = Response)) +
  geom_point(color = "blue", size = 2) +  # Observed responses
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), 
              fill = "gray", alpha = 0.3) +  # Confidence envelope
  geom_line(aes(y = CI_lower), color = "red", linetype = "dashed") +  # Lower CI
  geom_line(aes(y = CI_upper), color = "red", linetype = "dashed") +  # Upper CI
  theme_minimal() +
  labs(
    title = expression(xi),    # Add xi as title
    x = "Estimate",            # x-axis label
    y = "True"                 # y-axis label
  ) +
  theme(
    plot.title = element_text(hjust = 0.45),  # Center the title
    plot.margin = unit(c(2, 2, 2, 2), "cm")  # Top, right, bottom, left margins
  )
# # # Save the plot as a PNG file
# ggsave(filename = "CI-shape.png",   # Name of the file
#   plot = plot_xi,             # Plot object to save
#   width = 8,                        # Width in inches
#   height = 6,                       # Height in inches
#   dpi = 300                         # Resolution in dots per inch
# )


## Kappa 

# Combine data into a data frame
envelope_data_lnKappa2 <- data.frame(
  lnKappa2_pred=logkappa2_pred,
  Response=log(storeParameterRep[, 2]-4),
  CI_lower=lnKappa2CI_lower_bound,
  CI_upper=lnKappa2CI_upper_bound
)

# ggplot2 envelope plot with mathematical expression
library(ggplot2)

plot_lnKappa2 <- ggplot(envelope_data_lnKappa2, aes(x=logkappa2_pred, y = Response)) +
  geom_point(color = "blue", size = 2) +  # Observed responses
  geom_ribbon(aes(ymin = lnKappa2CI_lower_bound, ymax = lnKappa2CI_upper_bound), 
              fill = "gray", alpha = 0.3) +  # Confidence envelope
  geom_line(aes(y = CI_lower), color = "red", linetype = "dashed") +  # Lower CI
  geom_line(aes(y = CI_upper), color = "red", linetype = "dashed") +  # Upper CI
  theme_minimal() +
  labs(
    title = expression('ln'~kappa^2),    # Add xi as title
    x = "Estimate",            # x-axis label
    y = "True"                 # y-axis label
  ) +
  theme(
    plot.title = element_text(hjust = 0.45),  # Center the title
    plot.margin = unit(c(2, 2, 2, 2), "cm")  # Top, right, bottom, left margins
  )


# Combine data into a data frame
envelope_data_lnTau2 <- data.frame(
  lntau2_pred=logtau2_pred,
  Response=log(storeParameterRep[, 3]),
  CI_lower=lnTau2CI_lower_bound,
  CI_upper=lnTau2CI_upper_bound
)

# ggplot2 envelope plot with mathematical expression
library(ggplot2)

plot_lnTau2 <- ggplot(envelope_data_lnTau2, aes(x=logtau2_pred, y = Response)) +
  geom_point(color = "blue", size = 2) +  # Observed responses
  geom_ribbon(aes(ymin = lnTau2CI_lower_bound, ymax = lnTau2CI_upper_bound), 
              fill = "gray", alpha = 0.3) +  # Confidence envelope
  geom_line(aes(y = CI_lower), color = "red", linetype = "dashed") +  # Lower CI
  geom_line(aes(y = CI_upper), color = "red", linetype = "dashed") +  # Upper CI
  theme_minimal() +
  labs(
    title = expression('ln'~tau^2),    # Add xi as title
    x = "Estimate",            # x-axis label
    y = "True"                 # y-axis label
  ) +
  theme(
    plot.title = element_text(hjust = 0.45),  # Center the title
    plot.margin = unit(c(2, 2, 2, 2), "cm")  # Top, right, bottom, left margins
  )

library(gridExtra)  # For arranging multiple plots
library(patchwork)
# Arrange all three plots in a row
combined_plot <- grid.arrange(plot_xi, plot_lnKappa2, plot_lnTau2, ncol = 3)

# Save the combined plot
ggsave(filename = "Combined_Plot_Panel.png",  # File name
       plot = combined_plot,                 # Combined plot
       width = 15,                           # Total width for 3 plots
       height = 5,                           # Height of the plots
       dpi = 300                             # Resolution
)

library(grid)

# Arrange all three plots in a row with reduced spacing
combined_plot <- grid.arrange(
  arrangeGrob(
    plot_xi,
    plot_lnKappa2,
    plot_lnTau2,
    ncol = 3,
    plot_layout(widths = c(4, -1, 4.5)) &
      theme(legend.position = "top")))

# Save the combined plot
ggsave(
  filename = "~/Desktop/Combined_Plot_Panel.png",  # File name
  plot = combined_plot,                 # Combined plot
  width = 18,                           # Total width for 3 plots
  height = 6,                           # Height of the plots
  dpi = 300                             # Resolution
)


## panel plot:
# Plot for 'xi'
plot(
  envelope_data_shape$shape_pred,
  envelope_data_shape$Response,
  pch = 19,                 # Solid points
  col = alpha(alpha=0.4, "blue"),             # Point color
  xlab = "Estimate",        # x-axis label
  ylab = "True",            # y-axis label
  main = expression(xi),    # Title with mathematical expression
  ylim = range(c(envelope_data_shape$CI_lower, envelope_data_shape$CI_upper)) # Set y limits
)

# Add confidence ribbon
polygon(
  x = c(envelope_data_shape$shape_pred, rev(envelope_data_shape$shape_pred)),
  y = c(envelope_data_shape$CI_lower, rev(envelope_data_shape$CI_upper)),
  col = rgb(0.5, 0.5, 0.5, 0.3),  # Semi-transparent gray
  border = NA                     # No border
)

# Add lower and upper confidence interval lines
lines(envelope_data_shape$shape_pred, envelope_data_shape$CI_lower, col = alpha(alpha=0.5, "red"), lty = 3)
lines(envelope_data_shape$shape_pred, envelope_data_shape$CI_upper, col = alpha(alpha=0.5, "red"), lty = 3)

