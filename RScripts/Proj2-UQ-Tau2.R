## Confidence bound: log Tau
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


plot(0.5*log(tau2_pred_test),
     0.5*log(tau2_true_test),
     pch=19,
     col=alpha(alpha=0.4,'grey'),
     main=expression('Before update log('~tau~')'))
abline(0,1, col='magenta')


X_tau2_train <- data_frame(shape=shape_pred_train,
                           kappa=0.5*log(kappa2_pred_train),
                           tau=0.5*log(tau2_pred_train))

dim(X_tau2_train)

# Quantile levels
quantiles <- c(0.025, 0.975)

# Fit models for each quantile
modelsTau2 <- lapply(quantiles, function(tau){
  rq(0.5*log(tau2_true_train)~., data=X_tau2_train, tau = tau)
})



## Test: Summarize the models
X_tau2_test <- data_frame(shape=shape_pred_test,
                          kappa=0.5*log(kappa2_pred_test),
                          tau=0.5*log(tau2_pred_test))
dim(X_tau2_test)

fittedTau2Test <- lapply(modelsTau2, function(model){
  predict(model, newdata=X_tau2_test)})

QRTau2EstTest <- matrix(NA,
                        nrow=length(quantiles),
                        ncol=nrow(X_tau2_test))
for(i in 1:length(quantiles))
{
  QRTau2EstTest[i, ] <- fittedTau2Test[[i]]
}
dim(QRTau2EstTest)

## Compute lower and upper bound:
Tau2CI_lower_bound_test <- rep(NA, ncol(QRTau2EstTest))
Tau2CI_upper_bound_test <- rep(NA, ncol(QRTau2EstTest))
for(i in 1:ncol(QRTau2EstTest))
{
  Tau2CI_lower_bound_test[i] <-  QRTau2EstTest[1, i]
  Tau2CI_upper_bound_test[i] <-  QRTau2EstTest[2,i]
}


##  Confidence interval for  Training set
envelope_data_Tau2_test <- data.frame(Tau2_pred=0.5*log(tau2_pred_test),
                                      Response=0.5*log(tau2_true_test),
                                      CI_lower=Tau2CI_lower_bound_test,
                                      CI_upper=Tau2CI_upper_bound_test)

write.csv(envelope_data_Tau2_test,
          "~/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/envelope_data_Tau2_test.csv",
          row.names=FALSE)


plot_Tau2 <- ggplot(envelope_data_Tau2_test, aes(x=Tau2_pred, y=Response)) +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper), 
              fill = "lightblue", alpha=0.4) + 
  geom_point(color="navy", size=0.8, alpha=0.4) +
  geom_point(aes(y=CI_lower), color="lightblue", size=0.8) +
  geom_point(aes(y=CI_upper), color="lightblue", size=0.8) +
  geom_abline(intercept=0, slope=1, color="magenta", linetype="dashed", size=1) +
  theme_minimal() +
  labs(
    title=expression("95% Confidence Interval of log("~tau~")"),
    x="Estimates",
    y="True Values") +
  theme(
    plot.title=element_text(hjust=0.5, face="bold"),
    legend.position="top"
  )
print(plot_Tau2)

