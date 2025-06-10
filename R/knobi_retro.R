#' @title KBPM retrospective analysis
#'
#' @name knobi_retro
#'
#' @description This function performs a retrospective analysis that evaluates the robustness of the KBPM fit to the systematic deletion of recent data.
#'
#' @param knobi_results The object resulting from the \code{\link{knobi_fit}} function, which are the results of the base KBPM fit.
#' @param nR Number of retrospective peels to fit the model peeling off from 1 to nR years of data. 5 by default. See details.
#' @param yR Optional. A vector representing the final years of the catch time series for each of the retrospective models. See details.
#' @param yR0 Optional. A vector representing the starting years of the catch time series for each of the retrospective models. This vector must be the same length as the yR vector. By default, the catch time series is assumed to start in the same year as the original fit.
#' @param env_results Optional. The object resulting from the \code{\link{knobi_env}} function. A list containing the results of the environmental KBPM fit. If this argument is provided, the retrospective analysis is carried out for the base and the environmental KBPM models.
#' @param plot_out Logical. If TRUE, a file containing the plot of the retrospective fits is created. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_dir Optional. Directory to create the folder and save the plots. Required when 'plot_out=TRUE'. The default value is taken from the input of the \code{\link{knobi_fit}} function.
#' @param plot_filename Optional. Name of the folder that will contain the plots. Required when 'plot_out=TRUE'. The default value is taken from the input of the \code{\link{knobi_fit}} function.
#'
#' @details There are different options for defining retrospective fits:  
#' 
#' (1) Usage of \code{nR} argument. This argument specifies the number of retrospective peels. By using this argument, it is implied that the retrospective fits will consist of systematically deleting the last year of data, up to the number of years specified by \code{nR}.  
#' 
#' (2) Usage of \code{yR} argument. This argument specifies the final years of the catch time series for each of the retrospective models, providing greater flexibility in choosing the years from which to delete information. The number of retrospective fits will correspond to the length of the \code{yR} vector. Additionally, different starting years can be set using the \code{yR0} argument. 
#' 
#' If both arguments, \code{nR} and \code{yR}, are provided, the package will prioritize the use of \code{yR}.
#' 
#' As described in the \code{\link{knobi_env}} function details, in the case of the environmental models, both the estimated biological reference points and the plotted production curve correspond to a value of the scaled environmental variable equal to the mean of the time series, i.e. \eqn{X_{t}=0}, which cancels out the environmental effect in the equations defining both models. For more details, such as the calculation of BRPs as a function of the environmental variable, see vignettes.
#'
#' @return A list containing the results of the retrospective analysis, including parameter estimates and reference points for each model. The estimated surplus production curves from the retrospective analysis are also plotted, with a panel where each graph represents the curves for each model in case of considering environmental models. The plot is displayed in the plot window and saved (if plot_out=TRUE) in the specified directory or in the current directory.
#'
#' @author
#' \itemize{
#' \item{Anxo Paz}
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{M. Grazia Pennino}
#' }
#'
#' @examples
#'
#' \donttest{
#'
#' library(knobi)
#'
#' # See knobi_fit example to obtain the knobi_results object
#' knobi_retrospectives<-knobi_retro(knobi_results,plot_out=T)  # default nR=5
#' knobi_retrospectives
#'
#' knobi_retro(knobi_results,nR=3)
#' knobi_retro(knobi_results,yR=c(2010,2015))
#' knobi_retro(knobi_results,yR=c(2010,2015),yR0=c(1995,2000))
#'
#' # See knobi_env example to obtain the env_results object
#' knobi_retro(knobi_results,env_results=knobi_environmental,yR=c(2010,2015),yR0=c(1995,2000))
#'
#' }
#'
#' @export


utils::globalVariables(c("y","B","SP"))

knobi_retro <- function( knobi_results, env_results=NULL, nR=5, yR=NULL, yR0=NULL, plot_out=F, plot_filename=NULL, plot_dir=NULL){
  
  df <- knobi_results$df
  
  Year <- df$Year
  lastyear <- max(Year)
  
  x <- df$x
  
  pars <- knobi_results$params
  
  pella <- knobi_results$control$pella
  
  
  if(plot_out==T){
    
    old_dir <- getwd()
    on.exit( setwd(old_dir))
    
    if (is.null(plot_dir)) plot_dir <- knobi_results$control$plot_settings$plot_dir
    setwd(plot_dir)
    
    if (is.null(plot_filename)) plot_filename <- knobi_results$control$plot_settings$plot_filename
    if ( !plot_filename %in% list.dirs( full.names = FALSE)) dir.create( plot_filename)
    
    setwd(paste0(plot_dir,"/",plot_filename))
    
  }
  
  r <- pars[['r']]
  K <- pars[['K']]
  if ( pella){ p <- pars[['p']]}
  
  av <- seq( 0, K, length.out = 3*length(x))
  bv <- predict_model( knobi_results, av, pella, 'Base')
  
  df_aux <- data.frame( av, bv)
  
  
  if(is.null(yR)){
    
    nR <- nR
    yR <- lastyear-1:nR
    yR0 <- rep(Year[1],length(yR))
    
  } else {
    
    yR <- yR
    nR <- length(yR)
    if( is.null(yR0)) yR0<-rep( Year[1], length(yR))
    if( length(yR)!=length(yR0)) stop("yR and yR0 must have the same length")
    
  }
  
  modelretro <- list()
  brps <- NULL
  params <- NULL
  names_retro <- NULL
  
  df_plot <- data.frame( B=df_aux$av, SP=df_aux$bv, factor=paste(min(df$Year), "-", lastyear))
  
  for (i in 1:nR){
    
    newdf <- subset( df, Year<=yR[i] & Year>=yR0[i])
    iname <- paste( yR0[i], "-", yR[i])
    names_retro <- c( names_retro, iname)
    
    Data <- list(data=newdf, start_r=0.5, start_K=max(newdf$x), start_p=1)
    
    ipars <- kbpm_fit( Data, pella = pella, model = "Base")$par
    attr( ipars,'status') <- NULL
    params <- rbind( params, ipars)
    fit <- list( params = ipars)
    
    ibrps <- BRP( fit, pella)
    fit$brps <- ibrps
    
    brps <- rbind( brps, ibrps)
    rownames(brps)[i] <- rownames(params)[i] <- iname
    
    modelretro[[iname]] <- fit
    
    iK <- ipars['K']
    aretro <- seq(0, iK, length.out = 3*length(x))
    bretro <- predict_model( fit, aretro, pella, 'Base')
    
    i_df_plot <- data.frame( B=aretro, SP=bretro, factor=iname)
    df_plot <- rbind( df_plot, i_df_plot)
    
  }
  
  
  brps <- rbind( knobi_results$BRPs, brps)
  params <- rbind( knobi_results$params, params)
  
  rownames(brps)[1] <- rownames(params)[1] <- paste(min(df$Year), "-", lastyear)
  
  Retros <- list( BRPs = brps, params = params)
  
  
  if( is.null(env_results)){
    
    df_plot$factor <- as.factor( df_plot$factor)
    
    if( knobi_results$control$method == "SSB"){ baxis<-"SSB"} else { baxis<-"Biomass"}
    
    
    max_y <- max( df_plot$SP, df$y)
    min_y <- min( df_plot$SP, df$y)
    
    retro_plot <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::ylim(min_y,max_y) +
      ggplot2::geom_point( data=df[c(1,nrow(df)),], ggplot2::aes( x=x, y=y), color='darkblue', size=3, show.legend=FALSE) +
      ggplot2::geom_text( data=df[c(1,nrow(df)),], ggplot2::aes( x=x, y=y, label=Year), color='darkblue', vjust=-1, size=4, show.legend = FALSE) +
      ggplot2::geom_point( data=df, ggplot2::aes( x=x, y=y), color='darkblue', show.legend=FALSE) +
      ggplot2::geom_path( data=df, ggplot2::aes( x=x, y=y), color='darkblue') +
      ggplot2::labs( title="Retrospective Surplus Production Curves", subtitle=knobi_results$input$Stock, x=baxis, y="SP") +
      ggplot2::geom_line( data=df_plot, ggplot2::aes( x=B, y=SP, color=factor), linewidth=1) +
      ggplot2::scale_color_manual( values=c(2:(nR+1),1)) +
      ggplot2::guides( size="none",col=ggplot2::guide_legend(title="")) +
      ggplot2::theme( legend.position = c(.89,0.75), legend.background = ggplot2::element_rect(fill = "transparent"),
                      plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5),
                      axis.line=ggplot2::element_line())
    
    print(retro_plot)
    
    if(plot_out==TRUE){
      
      p <- grDevices::recordPlot()
      grDevices::jpeg("fits_retro.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(p)
      grDevices::dev.off()
      
      setwd(old_dir)
      
    }
    
  } else {
    
    Retros <- list( base = list( BRPs = brps, params = params))
    
    df <- env_results$df
    evar <- env_results$selected_var
    
    df_plot$Model='Base'
    
    addres <- env_results$add
    multres <- env_results$mult
    
    add_brps <- mult_brps <- NULL
    add_params <- mult_params <- NULL
    
    for (i in 1:nR){
      
      newdf <- subset( df, Year<=yR[i] & Year>=yR0[i])
      iname <- names_retro[i]
      
      ienv <- as.matrix( newdf[,evar])
      
      Data <- list( data=newdf, start_r=0.5, start_K=max(newdf$x), start_p=1, start_c=1)
      
      iadd_pars <- kbpm_fit( Data, ienv, pella, 'Add')
      imult_pars <- kbpm_fit( Data, ienv, pella, 'Mult')
      
      add_params <- rbind( add_params, iadd_pars)
      mult_params <- rbind( mult_params, imult_pars)
      
      iadd_fit <- list( params = iadd_pars)
      imult_fit <- list( params = imult_pars)
      
      iadd_brps <- BRP( iadd_fit, pella)
      imult_brps <- BRP( imult_fit, pella)
      
      add_brps <- rbind( add_brps, iadd_brps)
      mult_brps <- rbind( mult_brps, imult_brps)
      
      rownames(add_brps)[i] <- rownames(add_params)[i] <-
        rownames(mult_brps)[i] <- rownames(mult_params)[i] <- iname
      
      addiK <- iadd_fit$params['K']; multiK <- imult_fit$params['K']
      
      add_aretro <- seq(0, addiK, length.out = 3*length(x))
      add_bretro <- predict_model( iadd_fit, add_aretro, pella, 'Base')
      
      mult_aretro <- seq(0, multiK, length.out = 3*length(x))
      mult_bretro <- predict_model( imult_fit, mult_aretro, pella, 'Base')
      
      iadd_df_plot <- data.frame( B=add_aretro, SP=add_bretro, factor=iname, Model='Add')
      imult_df_plot <- data.frame( B=mult_aretro, SP=mult_bretro, factor=iname, Model='Mult')
      
      df_plot <- rbind( df_plot, iadd_df_plot, imult_df_plot)
      
    }
    
    
    add_brps <- rbind( addres$BRPs, add_brps)
    mult_brps <- rbind( multres$BRPs, mult_brps)
    
    add_params <- rbind( addres$params, add_params)
    mult_params <- rbind( multres$params, mult_params)
    
    rownames(add_brps)[1] <- rownames(add_params)[1] <-
      rownames(mult_brps)[1] <- rownames(mult_params)[1] <- paste(min(df$Year), "-", lastyear)
    
    Retros$add <- list( BRPs = add_brps, params = add_params)
    Retros$mult <- list( BRPs = mult_brps, params = mult_params)
    
    
    df_plot$factor <- as.factor( df_plot$factor)
    df_plot$Model <- as.factor( df_plot$Model)
    
    if( knobi_results$control$method == "SSB"){ baxis<-"SSB"} else { baxis<-"Biomass"}
    
    max_y <- max( df_plot$SP, df$y)
    min_y <- min( df_plot$SP, df$y)
    
    retro_plot <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::ylim(min_y,max_y) +
      ggplot2::geom_point( data=df[c(1,nrow(df)),], ggplot2::aes( x=x, y=y), color='darkblue', size=3, show.legend=FALSE) +
      ggplot2::geom_text( data=df[c(1,nrow(df)),], ggplot2::aes( x=x, y=y, label=Year), color='darkblue', vjust=-1, size=3, show.legend = FALSE) +
      ggplot2::geom_point( data=df, ggplot2::aes( x=x, y=y), color='darkblue', show.legend=FALSE) +
      ggplot2::geom_path( data=df, ggplot2::aes( x=x, y=y), color='darkblue') +
      ggplot2::labs( title="Retrospective Surplus Production Curves by Model", subtitle=knobi_results$input$Stock, x=baxis, y="SP") +
      ggplot2::geom_line( data=df_plot, ggplot2::aes( x=B, y=SP, color=factor), linewidth=1) +
      ggplot2::scale_color_manual( values=c(2:(nR+1),1)) +
      ggplot2::facet_grid(rows=dplyr::vars(Model), scales='free_y') +
      ggplot2::guides( size="none",col=ggplot2::guide_legend(title="")) +
      ggplot2::theme( legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                      plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5),
                      axis.line=ggplot2::element_line())
    
    print(retro_plot)
    
    if(plot_out==TRUE){
      
      p <- grDevices::recordPlot()
      grDevices::jpeg("fits_retro.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(p)
      grDevices::dev.off()
      
    }
    
  }
  
  return( Retros)
  
}

