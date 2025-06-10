#' @title Print a knobi object
#'
#' @name print.knobi
#'
#' @description The default print method for \code{\link{knobi_fit}}, \code{\link{knobi_env}} and \code{\link{knobi_proj}} object
#'
#' @param x,... Fitted model objects of class \code{knobi} produced by \code{knobi_fit()}, \code{knobi_env()} or \code{knobi_proj()}.
#'
#' @details Prints out the formula and the parameters estimates of the base KBPM fit or the environmental KBPM fit, furthermore it also reports the KBPM projections.
#'
#' @author
#' \itemize{
#' \item{Anxo Paz}
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{M. Grazia Pennino}
#' }
#'
#' @seealso
#' \code{\link{knobi_fit}}, \code{\link{knobi_env}}, \code{\link{knobi_proj}}
#'

utils::globalVariables(c("Sc"))

#' @export

print.knobi <- function(x, ...){

  if(!is.null(x$plots) & is.null(x$plots3D)){
    cat("\n Projections: \n \n")
    print(subset(x$df, Sc!='input'), row.names = F)
    cat("\n \n")

  } else if(is.null(x$mult)){

    if(x$control$pella==TRUE){
      cat("\n Formula:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p) \n \n")
    } else {
      cat("\n Formula:\n","SP_t = r*B_t*(1-B_t/K) \n \n")}

    cat("Parameter estimates:\n")

    cat( "r ",x$params['r'],"\n")
    cat( "K ",x$params['K'],"\n")

    if(x$control$pella==TRUE){
      cat("p ",x$params['p'],"\n \n")
    } else {cat("\n")}

  } else {

    if( x$input$basecontrol$pella){

      cat("\n Multiplicative model:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p)*exp(c*X_t) \n \n")
    } else {
      cat("\n Multiplicative model:\n","SP_t = r*B_t*(1-B_t/K)*exp(c*X_t) \n \n")
    }

    cat("Parameter estimates:\n")

    npms<-length(x$mult$params)

    for(i in 1:npms){
      cat(names(x$mult$params)[i]," ",x$mult$params[i],"\n")
    }

    cat("\n")

    if( x$input$basecontrol$pella){

      cat("\n Additive model:\n","SP_t = (r/p)*B_t*(1-(B_t/K)^p)+c*X_t*B_t \n \n")
    } else {
      cat("\n Additive model:\n","SP_t = r*B_t*(1-B_t/K)+c*X_t*B_t \n \n")
    }

    cat("Parameter estimates:\n")

    for(i in 1:npms){
      cat(names(x$add$params)[i]," ",x$add$params[i],"\n")
    }

    cat("\n")

  }

}
