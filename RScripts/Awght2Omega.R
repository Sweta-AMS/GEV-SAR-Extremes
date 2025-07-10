omega2Awght<- function( omega)
{
  # for a rectangle this should be:
  Awght <- 4 +  exp(omega)^2
  return( Awght )
}

Awght2Omega<- function( Awght)
{
  # for a rectangle this should be:
  #  Awght <- 4 +  exp(omega)^2
  omega  <- log(Awght-4)/2
}
