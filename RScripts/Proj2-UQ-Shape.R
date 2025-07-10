## Confidence bound: Xi
## Required script
rm(list=ls())
setwd("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")
source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")


shape_true_train <- trainParameters[ , 1]
kappa2_true_train <- (trainParameters[ , 2]-4)
tau2_true_train <- (trainParameters[ , 3])

shape_true_test <- testParameters[ , 1]
kappa2_true_test <- (testParameters[ , 2]-4)
tau2_true_test <- (testParameters[ , 3])

## Combine predictors into a data frame
X_shape_train <- data.frame(shape=shape_pred_train,
                            kappa=0.5*log(kappa2_pred_train),
                            tau=0.5*log(tau2_pred_train))


## -----------------------------------------------
# # Spline implementation
# X_shape_train<- ns(shape_pred_train, df=15)
# Combine predictors
# x1_train <- cbind(x1=x1_spline_train,
#                   x2=x2_spline_train,
#                   x3=x3_spline_train)
## -------------------------------------------------

## Quantile levels:
quantiles <- c(0.025, 0.975)

## Fit models for each quantile:
library(quantreg)
modelsShape <- lapply(quantiles, function(tau) {
 rq(shape_true_train ~ ., data = X_shape_train, tau=tau) # when dataframe under normal scenario
})
summary(modelsShape)

### Finding the Confidence bound on the Test
X_shape_test <- data.frame(shape=shape_pred_test,
                           kappa=0.5*log(kappa2_pred_test),
                           tau=0.5*log(tau2_pred_test))

fittedShapeTest <- lapply(modelsShape, function(model) {
  predict(model, newdata=X_shape_test)
})

QRShapeEstTest <- matrix(NA,
                         nrow=length(quantiles),
                         ncol=length(shape_pred_test))

for(i in 1:length(quantiles)) {
  QRShapeEstTest[i, ] <- fittedShapeTest[[i]]
}

# Checking dimensions correctly
dim(QRShapeEstTest)  


## Compute lower and upper bound:
XiCI_lower_bound_test <- rep(NA, ncol(QRShapeEstTest))
XiCI_upper_bound_test <- rep(NA, ncol(QRShapeEstTest))
for(i in 1:ncol(QRShapeEstTest))
{
  XiCI_lower_bound_test[i] <- QRShapeEstTest[1, i]
  XiCI_upper_bound_test[i] <- QRShapeEstTest[2,i]
}


##  Confidence interval for  Training set
envelope_data_Xi_test <- data.frame(shape_pred=shape_pred_test,
                                    Response=shape_true_test,
                                    CI_lower=XiCI_lower_bound_test,
                                    CI_upper=XiCI_upper_bound_test)

write.csv(envelope_data_Xi_test,
          "~/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/envelope_data_Xi_test.csv",
          row.names=FALSE)

plot_xi <- ggplot(envelope_data_Xi_test, aes(x=shape_pred, y=Response)) +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper), 
              fill = "lightblue", alpha=0.4) + 
  geom_point(color="navy", size=0.8, alpha=0.4) +
  geom_point(aes(y=CI_lower), color="lightblue", size=0.8) +
  geom_point(aes(y=CI_upper), color="lightblue", size=0.8) +
  geom_abline(intercept=0, slope=1, color="magenta", linetype="dashed", size=1) +
  theme_minimal() +
  labs(
    title=expression("95% Confidence Interval of "~xi),
    x="Estimates",
    y="True Values") +
  theme(
    plot.title=element_text(hjust=0.5, face="bold"),
    legend.position="top"
  )
print(plot_xi)