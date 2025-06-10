predict_model <- function( fit, x, pella, model, env = NULL){

  x <- x
  r <- fit$params['r']
  K <- fit$params['K']
  p <- ifelse( pella, fit$params['p'], 1)

  if( model != 'Base') cs <- fit$params[ifelse(pella,4,3):length(fit$params)]

  if ( model == "Mult") { y <- as.numeric(exp(as.matrix(env) %*% cs)*((r/p)*x*(1-(x/K)^(p))))
  } else if ( model == "Add") { y <- as.numeric((r/p)*x*(1-(x/K)^(p)) + (as.matrix(env) %*% cs)*x)
  } else { y <- (r/p)*x*(1-(x/K)^(p))}

  return(y)

}
