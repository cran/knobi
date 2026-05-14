#' @title knobi_fit results
#'
#' @description
#' 
#' Results of the \code{\link{knobi_fit}} first example, to be used in the other functions' exemples
#' 
#' @details The results can be obtained using the following R code:
#' 
#' \preformatted{
#' data(knobi_dataset)
#' hake_n <- knobi_dataset$hake_n
#' data<-list()
#' data$SSB<-hake_n$SSB
#' data$Catch<-hake_n$catches
#' data$F_input<-hake_n$F
#' data$RP<-list(F_MSY=0.259, B_MSY=207398, MSY=75052, K=NA)
#' data$years<-hake_n$Year
#' control<-list(pella = "TRUE") 
#' knobi_results<-knobi_fit(data,control)
#' }
#'
#' @usage data(knobi_results)

"knobi_results"