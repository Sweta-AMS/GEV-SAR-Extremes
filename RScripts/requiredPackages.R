### -- Package Required -- ###
suppressMessages({
  library("ggplot2")
  library("parallel")
  library("quantreg")
  library("scales")
  library("fields")
  library("mgcv")
  library("maps")
  library("dplyr")
  library("plyr")
  library("ggpubr")
  library("viridis") # viridis colour scale
  library("reshape2")
  library("rlist")
  library("colorspace")
  library("grDevices")
  library("LatticeKrig")
  
  ## Load ismev library for marginal GPD fits:
  library("ismev")
  library("evd")
  library("MASS")
  library("fitdistrplus")
  library("qmap")
  library("survival")
  library("tictoc")
  library("extRemes")
  library("splines")
  library("scales") # plotting against time
  
  library("reticulate")
  library("survival")
  library("RcppCNPy")
  
  ## Simulate r-Pareto Process: 
  library('mev')
  library('LatticeKrig')
  library('Matrix')
  library('lubridate')
})