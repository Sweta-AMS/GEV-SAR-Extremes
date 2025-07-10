## Confidence bound: log Kappa
## Required script
rm(list=ls())

setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")


## Loading the training parameter and test true parameters:
shape_true_train <- trainParameters[ , 1]
kappa2_true_train <- (trainParameters[ , 2]-4)
tau2_true_train <- (trainParameters[ , 3]) # 50000

shape_true_test <- testParameters[ , 1]
kappa2_true_test <- (testParameters[ , 2]-4)
tau2_true_test <- (testParameters[ , 3]) # 10000

plot(0.5*log(kappa2_pred_test),
     0.5*log(kappa2_true_test),
     pch=19,
     col=alpha(alpha=0.4,'grey'),
     main=expression('Before update log('~kappa~')'))
abline(0,1, col='magenta')

## Using all the estimate parameters to fit the CI curve:
X_kappa2_train <- data_frame(shape=shape_pred_train,
                             kappa=0.5*log(kappa2_pred_train),
                             tau=0.5*log(tau2_pred_train))
  
  #data_frame(kappa=0.5*log(kappa2_pred_train),
  #                           kappa2=(0.5*log(kappa2_pred_train))^2,
  #                          tau=0.5*log(tau2_pred_train))
dim(X_kappa2_train)


# Quantile levels
quantiles <- c(0.025, 0.975) 

# Fit models for each quantile
modelsKappa2 <- lapply(quantiles, function(tau){
  rq(0.5*log(kappa2_true_train) ~., data=X_kappa2_train, tau=tau)
})

summary(modelsKappa2[[1]])

## Test: Summarize the models
X_kappa2_test <- data_frame(shape=shape_pred_test,
                            kappa=0.5*log(kappa2_pred_test),
                            tau=0.5*log(tau2_pred_test))
  
  #data_frame(kappa=0.5*log(kappa2_pred_test),
  #                          kappa2=(0.5*log(kappa2_pred_test))^2,
  #                          tau=0.5*log(tau2_pred_test))
dim(X_kappa2_test)


fittedKappa2Test <- lapply(modelsKappa2, function(model){
  predict(model, newdata=X_kappa2_test)})

QRKappa2EstTest <- matrix(NA,
                          nrow=length(quantiles),
                          ncol=nrow(X_kappa2_test))
for(i in 1:length(quantiles))
{
  QRKappa2EstTest[i, ] <- fittedKappa2Test[[i]]
}
dim(QRKappa2EstTest)

## Compute lower and upper bound:
Kappa2CI_lower_bound_test <- rep(NA, ncol(QRKappa2EstTest))
Kappa2CI_upper_bound_test <- rep(NA, ncol(QRKappa2EstTest))
for(i in 1:ncol(QRKappa2EstTest))
{
  Kappa2CI_lower_bound_test[i] <- QRKappa2EstTest[1, i]
  Kappa2CI_upper_bound_test[i] <- QRKappa2EstTest[2,i]
}


##  Confidence interval for  Training set
envelope_data_Kappa2_test <- data.frame(Kappa2_pred=0.5*log(kappa2_pred_test),
                                        Response=0.5*log(kappa2_true_test),
                                        CI_lower=Kappa2CI_lower_bound_test,
                                        CI_upper=Kappa2CI_upper_bound_test)

write.csv(envelope_data_Kappa2_test,
          "~/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/envelope_data_Kappa2_test.csv",
          row.names=FALSE)

plot_Kappa2 <- ggplot(envelope_data_Kappa2_test, aes(x=Kappa2_pred, y=Response)) +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper), 
              fill = "lightblue", alpha=0.4) + 
  geom_point(color="navy", size=0.8, alpha=0.4) +
  geom_point(aes(y=CI_lower), color="lightblue", size=0.8) +
  geom_point(aes(y=CI_upper), color="lightblue", size=0.8) +
  geom_abline(intercept=0, slope=1, color="magenta", linetype="dashed", size=1) +
  theme_minimal() +
  labs(
    title=expression("95% Confidence Interval of log("~kappa~")"),
    x="Estimates",
    y="True Values") +
  theme(
    plot.title=element_text(hjust=0.5, face="bold"),
    legend.position="top"
  ) 
print(plot_Kappa2)

# length(kappa2_pred_test)
# plot(0.5*log(kappa2_pred_test),
#      0.5*log(kappa2_true_test),
#      pch=19,
#      col=alpha(alpha=0.4,'grey'),
#      main=expression('Before update log('~kappa~')'))
# abline(0,1, col='magenta')
# 
# bplot.xy(0.5*log(kappa2_true_test),
#          0.5*log(kappa2_pred_test)- 0.5*log(kappa2_true_test),
#          #Kappa2PredBCTest- 0.5*log(kappa2_true_test),
#          pch=19,
#          xlab= 'True',
#          ylab= 'Residuals', 
#          col=alpha(alpha=0.4,'grey'),
#          main=expression('After update log('~kappa~')'))
# abline(h=0, col='magenta')
# 
# ### 
# # Load necessary libraries
# library(splines)   # For smoothing splines
# 
# # Compute residuals
# residual <- 0.5 * log(kappa2_pred_test) -  0.5 * log(kappa2_true_test)
# 
# # Compute transformed predictor
# log_kappa_true <- 0.5 * log(kappa2_true_test)
# 
# # Fit a smoothing spline to residuals
# fit_bC_kappa <- smooth.spline(x=log_kappa_true,
#                               y=residual,
#                               spar = 0.8)  # Adjust 'spar' for smoothness
# 
# # Extract fitted residuals from the spline model
# residualF <- predict(fit_bC_kappa,
#                      x=log_kappa_true)$y  # Smoothed residuals
# 
# 
# 
# # Plot residuals after smoothing
# bplot.xy(log_kappa_true,
#          0.5 *log(kappa2_pred_test), 
#          pch = 19,
#          col = alpha('grey', 0.4),
#          xlab = expression('True log('~kappa~')'), 
#          ylab = "Estimated",
#          main = expression('Before bias correciton log('~kappa~')')
# )
# abline(0, 1, col = 'magenta', lwd = 2)  # Reference line at 0
# 
# # Plot residuals after smoothing
# bplot.xy(log_kappa_true,
#          0.5 *log(kappa2_pred_test) - residualF, 
#          pch = 19,
#          col = alpha('grey', 0.4),
#          xlab = expression('True log('~kappa~')'), 
#          ylab = "Estimated",
#          main = expression('After bias correction log('~kappa~')')
# )
# abline(0, 1, col = 'magenta', lwd = 2)  # Reference line at 0
# 
# # Plot residuals after smoothing
# bplot.xy(log_kappa_true,
#          0.5 *log(kappa2_pred_test) - log_kappa_true,
#          pch = 19,
#          col = alpha('grey', 0.4),
#          xlab = expression('True '~xi), 
#          ylab = "Residuals",
#          main = expression('Before bias correction '~xi)
# )
# abline(h=0,
#        col='magenta',
#        lwd=2)  # Reference line at 0
# lines(fit_bC_kappa,
#       col="blue", 
#       lwd=2)  # Add smoothed spline line
# 
# 
# # Plot residuals after smoothing
# bplot.xy(log_kappa_true,
#          0.5 *log(kappa2_pred_test) - residualF - log_kappa_true,
#          pch = 19,
#          col = alpha('grey', 0.4),
#          xlab = expression('True '~xi), 
#          ylab = "Residuals",
#          main = expression('Before bias correction '~xi)
# )
# abline(h=0,
#        col='magenta',
#        lwd=2)  # Reference line at 0

