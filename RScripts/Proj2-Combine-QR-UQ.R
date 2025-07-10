rm(list=ls())

## Required script
setwd("/Users/srai@mines.edu/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts")

source("requiredPackages.R")
source("LKrigSimSAR.R")
source("LKrigSAREvd.R")
source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Proj2-Data2Evaluate.R")

# source("~/Desktop/Research/Proj2-SAR/RScripts/Proj2-UQ-Shape.R") # Shape
# source("~/Desktop/Research/Proj2-SAR/RScripts/Proj2-UQ-Kappa2.R") # Kappa2
# source("~/Desktop/Research/Proj2-SAR/RScripts/Proj2-UQ-Tau2.R")  # Tau2

library(gridExtra)  
library(patchwork)

## Xi:
# Read the CSV file
shapeDF <- read.csv("~/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/envelope_data_Xi_test.csv", 
                    header=TRUE)

# View the first few rows
head(shapeDF)
dim(shapeDF)

envelope_data_Xi_test <- data.frame(shape_pred=shapeDF$shape_pred, #  XiPredBCTest
                                    Response=shapeDF$Response,     # shape_true_test
                                    CI_lower= shapeDF$CI_lower,    # XiCI_lower_bound_test,
                                    CI_upper=shapeDF$CI_upper)     # XiCI_upper_bound_test
plot_xi <- ggplot(envelope_data_Xi_test, aes(x=shape_pred, y=Response)) +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper), 
              fill = "lightblue", alpha=0.2) + 
  geom_point(color="navy", size=0.8, alpha=0.2) +
  geom_point(aes(y=CI_lower), color="lightblue", size=1) +
  geom_point(aes(y=CI_upper), color="lightblue", size=1) +
  geom_abline(intercept=0, slope=1, color="magenta", linetype="dashed", size=1) +
  theme_minimal() +
  labs(
    title=expression("95% Confidence Interval of "~xi),
    x="Estimates",
    y="True Values") +
  theme(
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    panel.border = element_rect(color="grey", fill=NA, size=1),  # Adds a black border around the panel
    plot.title=element_text(size=15, hjust=0.5, face="bold"),
    legend.position="top",
    axis.title.x = element_text(size=13),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=13), 
    axis.text.y = element_text(size=13)
  )
# theme(
#   text = element_text(size=20),        # Base text size
#   axis.title.x = element_text(size=24),  # X-axis label size
#   axis.title.y = element_text(size=24),  # Y-axis label size
#   axis.text.x = element_text(size=18),   # X-axis tick labels
#   axis.text.y = element_text(size=18),   # Y-axis tick labels
#   plot.title = element_text(size=26, face="bold", hjust=0.5)  # Title size, bold, centered
# )
print(plot_xi)

## Kappa 
kappaDF <- read.csv("~/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/envelope_data_Kappa2_test.csv", 
                    header=TRUE)

# View the first few rows
head(kappaDF)
dim(kappaDF)

##  Confidence interval for  Training set
envelope_data_Kappa2_test <- data.frame(Kappa2_pred=kappaDF$Kappa2_pred,  # Kappa2PredBCTest
                                        Response=kappaDF$Response,        # 0.5*log(kappa2_true_test),
                                        CI_lower=kappaDF$CI_lower,        # Kappa2CI_lower_bound_test,
                                        CI_upper=kappaDF$CI_upper)        #Kappa2CI_upper_bound_test)

plot_Kappa2 <- ggplot(envelope_data_Kappa2_test, aes(x=Kappa2_pred, y=Response)) +
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper), 
              fill = "lightblue", alpha = 0.4) +  # Transparent shading
  geom_point(color="navy", size=0.8, alpha=0.2) +
  geom_point(aes(y=CI_lower), color="lightblue", size=1) +
  geom_point(aes(y=CI_upper), color="lightblue", size=1) +
  geom_abline(intercept=0, slope=1, color="magenta", linetype="dashed", size=1) +
  theme_minimal() +
  labs(
    title=expression("95% Confidence Interval of log("~kappa~")"),
    x="Estimates",
    y="True Values") +
  theme(
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    panel.border = element_rect(color="grey", fill=NA, size=1),  # Adds a black border around the panel
    plot.title=element_text(size=15, hjust=0.5, face="bold"),
    legend.position="top",
    axis.title.x = element_text(size=13),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=13), 
    axis.text.y = element_text(size=13)
  )
print(plot_Kappa2)


## Tau
tauDF <- read.csv("~/Desktop/Research/Proj2-SAR/GEV-SAR/Result_UQ_BPs/envelope_data_Tau2_test.csv", 
                  header=TRUE)

# View the first few rows
head(tauDF)
dim(tauDF)

envelope_data_Tau2_test <- data.frame(Tau2_pred=tauDF$Tau2_pred,   #Tau2PredBCTest
                                      Response=tauDF$Response,     # 0.5*log(tau2_true_test),
                                      CI_lower=tauDF$CI_lower,     # Tau2CI_lower_bound_test,
                                      CI_upper=tauDF$CI_upper)     # Tau2CI_upper_bound_test)

plot_Tau2 <- ggplot(envelope_data_Tau2_test, aes(x=Tau2_pred, y=Response)) +
  # Shaded confidence interval
  geom_ribbon(aes(ymin=CI_lower, ymax=CI_upper), 
              fill = "lightblue", alpha=0.2) +  # Transparent shading
  
  # Scatter points for observed responses (True vs Predicted)
  geom_point(color="navy", size=0.8, alpha=0.2) +
  # Lower and upper confidence bounds
  geom_point(aes(y=CI_lower), color="lightblue", size=1) +
  geom_point(aes(y=CI_upper), color="lightblue", size=1) +
  # Perfect agreement line (y = x)
  geom_abline(intercept=0, slope=1, color="magenta", linetype="dashed", size=1) +
  # Customize theme
  theme_minimal() +
  labs(
    title=expression("95% Confidence Interval of log("~tau~")"),
    x="Estimates",
    y="True Values") +
  theme(
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),  # Removes minor grid lines
    panel.border = element_rect(color="grey", fill=NA, size=1),  # Adds a black border around the panel
    plot.title=element_text(size=15, hjust=0.5, face="bold"),
    legend.position="top",
    axis.title.x = element_text(size=13),
    axis.title.y = element_text(size=13),
    axis.text.x = element_text(size=13), 
    axis.text.y = element_text(size=13)
  )
print(plot_Tau2)

# Arrange all three plots in a row
combined_plot <- grid.arrange(plot_xi,
                              plot_Kappa2,
                              plot_Tau2,
                              ncol=3)

# Save the combined plot
ggsave(filename = "Combined_Plot_Panel.png",  # File name
  plot=combined_plot,                 # Combined plot
  width=14,                           # Total width for 3 plots
  height=4.5,                           # Height of the plots
  dpi=300                             # Resolution
)

# library(grid)
# 
# # Arrange all three plots in a row with reduced spacing
# combined_plot <- grid.arrange(
#   arrangeGrob(
#     plot_xi,
#     plot_Kappa2,
#     plot_Tau2,
#     ncol = 3,
#     plot_layout(widths = c(4, -1, 4.5)) &
#     theme(legend.position = "top")))
# 
# # Save the combined plot
# ggsave(
#   filename = "~/Desktop/Combined_Plot_Panel.png",  # File name
#   plot = combined_plot,                 # Combined plot
#   width = 18,                           # Total width for 3 plots
#   height = 6,                           # Height of the plots
#   dpi = 300                             # Resolution
# )
# 

