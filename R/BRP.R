BRP <- function( fit, pella){

  r <- fit$params['r']
  K <- as.numeric( fit$params['K'])
  p <- ifelse( pella, fit$params['p'], 1)

  Bmsy <- as.numeric(K*(1/(p+1))^(1/p))
  MSY<-as.numeric(r*K*(1/(p+1))^(1+1/p))
  x<-as.numeric(MSY/K)
  Fmsy <- as.numeric((r/p)*(1-(1/(p+1))))

  BRPs <- c( K=K, B_MSY=Bmsy, F_MSY=Fmsy, MSY=MSY, MSYoverK=x)
  return(BRPs)
}
