LKrigSARGaussian<- function(LKinfo,M,
                       asList=FALSE ){
  coefSAR<- NULL
  nLevel<- LKinfo$nlevel
  for( K in 1:nLevel){
    B <- LKrigSAR(LKinfo,Level=K)
    B<- spind2spam(B) # convert to spam format 
    mLevel<- LKinfo$latticeInfo$mLevel[K]
    # here using Gaussian but can change to other distribution 
    U<- runif(mLevel*M )
    E<-  qnorm(U)
    E<- matrix(E, mLevel, M)
    #
    y<- as.matrix(solve(B, E)) # uses sparse methods
    if(!asList){
    coefSAR<- rbind(coefSAR, y)
    }
    else{ coefSAR<- c( coefSAR, list(y))}
  }
  return( coefSAR)
}