LKrigSAREvd<- function(LKinfo, M,
                       asList=FALSE,
                       loc=0, scale=1, shape=0){
  coefSAR<- NULL
  nLevel<- LKinfo$nlevel
  B_list <- list() 
  for(K in 1:nLevel){
    B <- LKrigSAR(LKinfo, Level=K)
    B <- spind2spam(B) # convert to spam format 
    mLevel <- LKinfo$latticeInfo$mLevel[K] # from LatticeKrig package
    
    # here using Gaussian but can change to other distribution 
    U <- runif(mLevel*M )
    E <-  qevd(U,
              loc = loc, scale = scale,
              shape=shape, type="GEV")
    E <- matrix(E, mLevel, M)

    c <- as.matrix(solve(B, E)) # uses sparse methods
    
    if(!asList){
    coefSAR <- rbind(coefSAR, c)
    }
    else{ coefSAR <- c(coefSAR, list(c))}
    
    B_list[[K]] <- B # Store B matrix at each level
  }
  return(list('coefSAR'= coefSAR, 'SAR matrix'= B_list))
}
