LKrigSimSAR<- function(LKinfo, M, rep=m){
  coefSAR<- NULL
  nLevel<- LKinfo$nlevel
  for( K in 1: nLevel){
    B <- LKrigSAR(LKinfo,Level=K)
    B<- spind2spam(B) # convert to spam format 
    mLevel<- LKinfo$latticeInfo$mLevel[K]
    # here using Gaussian but can change to other distribution 
    E<- matrix(
      rnorm(mLevel*M), mLevel, M)
    #
    coefTmp<- solve(B, E) # uses sparse methods
    coefSAR<- rbind(coefSAR, coefTmp)
  }
  return( coefSAR)
}