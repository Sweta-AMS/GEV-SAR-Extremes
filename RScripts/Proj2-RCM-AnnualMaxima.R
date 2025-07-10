source("~/Desktop/Research/Proj2-SAR/GEV-SAR/RScripts/Read-RCM-Data.R")
setwd("~/Desktop/Tentative-Reslt")

png("RCM.png",
    units="in",
    width=21,
    height=7,
    res=200)
par(mfrow=c(1,3),
    mar=c(4,4,4,4) + 0.3,
    oma=c(0.4,2.3,2.8,4))
set.seed(123)
t <- round(runif(n=3, 1, 31))
for(i in 1:length(t))
{
  ind <- t[i]
  data <- (annual_maxima[ , ,ind] - median(annual_maxima[,,ind]))/sd(annual_maxima[ , ,ind])
  image.plot(lon_data,
             lat_data,
             data,
             zlim=c(-1.4, quantile(data, 0.999)),
             yaxt='n',
             xaxt='n',
             xlab='',
             ylab= '',
             horizontal=TRUE,
             legend.width=1.2,
             cex.main=7,
             axis.args=list(cex.axis=3,padj=0.9, lwd=2))
  axis(2,
       las=1,
       cex.axis=3,
       lwd=2,
       padj=0.9)
  axis(1,
       cex.axis=3,
       lwd=2,
       padj=0.9)
}
dev.off()