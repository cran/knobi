#' @title KBPM projections
#'
#' @name knobi_proj
#'
#' @description This function forecasts population and fishery dynamics under a defined set of management targets. It projects the time series of biomass (or SSB) and the surplus production for a set of future catch or fishing mortality scenarios.
#'
#' @param knobi_results The output object of \code{\link{knobi_fit}} function (main package function).
#' @param env_results Optional. The output object of \code{\link{knobi_env}} function.
#' @param c Optional. A vector, data frame, or matrix specifying the catch values for each projected year. Multiple catch scenarios can be defined, with each column in the data frame or matrix representing a different scenario. For a single scenario, the length of the vector must match the number of projected years; for multiple scenarios, the number of rows must match the projected years. Projections are based on either future catch values or fishing mortality values, so only one of the two arguments, 'c' or 'f', is required.
#' @param f Optional. A vector, data frame, or matrix specifying the fishing mortality values for each projected year. Multiple fishing mortality scenarios can be defined, with each column in the data frame or matrix representing a different scenario. For a single scenario, the length of the vector must match the number of projected years; for multiple scenarios, the number of rows must match the projected years. Projections are based on either future catch values or fishing mortality values, so only one of the two arguments, 'c' or 'f', is required.
#' @param env Optional. Environmental variable(s) projections required if the environmental fit is provided to forecast the population and fishery dynamics. This fit considers the variable(s) selected in the \code{\link{knobi_env}} function. If the 'multicovar' argument of \code{\link{knobi_env}} is FALSE, this argument is a vector, data frame or matrix containing the values of the environmental covariates (unstandardized) (cols) for the projection years (rows). On the other hand, if the 'multicovar' argument of \code{\link{knobi_env}} is TRUE, the current argument must be a list, and each entry must be a data frame or matrix corresponding to each environmental scenario containing the values of the considered covariates for each scenario.
#' @param plot_out Logical. If TRUE, a file containing the plot of the retrospective fits is created. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_dir Optional. Directory to create the folder and save the plots. Required when 'plot_out=TRUE'. The default value is taken from the input of the \code{\link{knobi_fit}} function.
#' @param plot_filename Optional. Name of the folder that will contain the plots. Required when 'plot_out=TRUE'. The default value is taken from the input of the \code{\link{knobi_fit}} function.
#'
#' @return A list containing the projection results. \itemize{
#' \item df: data frame containing the projected time series. The rows correspond to the projection years, while the columns represent stock quantities: biomass or SSB, surplus production (SP), catch (C), and fishing mortality (F). Three additional columns specify the catch or fishing mortality scenario (Sc), the model used, and the environmental scenario (EnvSc).
#' \item plots: list containing the plots with the projections. Each plot is a panel plot that includes four subplots representing biomass or SSB, surplus production, fishing mortality, and catches for each catch or fishing mortality scenario and for each environmental scenario (if they is provided).}
#' The resulting plots are displayed in the plot window and are also saved (if plot_out="TRUE") in the  provided directory or in the same directory as \code{link{knobi_fit}}.
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
#' data(knobi_results) # loading results of the knobi_fit example
#' data(knobi_environmental) # loading results of the knobi_env example
#' data(knobi_dataset)
#'
#' ### Projecting through catch with no environmental information
#'
#' # Then, create the data frame containing the selected catch for the projected
#' # years. In this illustration, within each scenario, the catch values are
#' # constant through the projected years. Three scenarios are considered:
#' # (i) catch value equal to the last historical catch multiplied by 1,
#' # (ii) last historical catch multiplied by 1.2 and
#' # (iii) last historical catch multiplied by 0.8.
#'
#' catch<-rep(knobi_results$df$C[length(knobi_results$df$C)],5)
#'
#' C<-data.frame(catch=catch,
#'               catch08=0.8*catch,
#'               catch12=1.2*catch)
#'
#' # Then, knobi_proj function can be applied
#'
#' knobi_proj(knobi_results, c=C)
#'
#'
#' ### With environmental information
#'
#' # In this case, in addition to the previous example, the 'knobi_env' example
#' # has to be run at first, where AMO variable was selected in the fit
#'
#' # We include the future values of the environmental variable(s) in a data
#' # frame containing the environmental covariable values for the projected
#' # years. Three scenarios are considered:
#' # (i) Constant AMO equal to last year's AMO,
#' # (ii) Constant AMO equal to last year's AMO with a 50% increment
#' # (iii) Constant AMO equal to last year's AMO with a 50% decrease
#'
#' Env <- knobi_dataset$Env
#' last_AMO <- Env$AMO[length(Env$AMO)]
#' env <- data.frame( AMOi=rep(last_AMO,5),
#'                    AMOii=rep(last_AMO*1.5,5),
#'                    AMOiii=rep(last_AMO*0.5,5))
#'
#' # Based on the previous objects we can apply the projection function.
#'
#' knobi_proj(knobi_results, knobi_environmental, c=C, env=env)
#'
#'
#' ### Through fishing mortality without environmental information
#'
#' # Alternatively, projections can be based on fishing mortality.
#' # The scenarios presented below have been created from the estimated F_msy of
#' # knobi_fit analysis.
#'
#' fmsy<-knobi_results$BRPs['F_MSY']
#' ff<-rep(fmsy,8)
#' f<-data.frame(f=ff,f12=ff*1.2,f08=ff*0.8)
#'
#' knobi_proj(knobi_results, f=f)
#'
#'
#' ### Through fishing mortality with environmental information
#'
#' knobi_proj(knobi_results, f=f[1:5,], env_results=knobi_environmental, env=env)
#'
#'
#' # In case of multicovar<-TRUE in knobi_env, a list is required in which
#' # each item is a data frame for each environmental scenario
#'
#' # env<-list(climate_1=data.frame(AMO=c(0.2,0.2,0.3,0.3,0.4),
#' #                                NAO=c(0.2,0.2,0.3,0.3,0.4)),
#' #           climate_2=data.frame(AMO=c(0.2,0.3,0.4,0.5,0.6),
#' #                                NAO=c(0.2,0.2,0.3,0.3,0.4)))
#' # 
#' # knobi_proj(knobi_results, knobi_environmental2, c=C, env=env)
#' }
#'
#' @export


utils::globalVariables(c("B", "Sc", "SP","C","FM","Value","Variable","EnvSc"))

knobi_proj <- function( knobi_results, env_results=NULL, c=NULL, f=NULL, env=NULL,
                         plot_out=FALSE, plot_filename=NULL, plot_dir=NULL){


  years <- knobi_results$df$Year
  lastyear <- years[length(years)]

  pella <- knobi_results$control$pella

  if( is.null(c) & is.null(f)) stop('You must provide catch or f time series')

  # models ---------------

  if( is.null(f)){ byc <- T
    model <- function(Bt1,Bt,Xt,K,r,p) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))-Xt-Bt1
    model_a <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))+Et%*%c*Bt-Xt-Bt1
    model_m <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))*exp(Et%*%c)-Xt-Bt1
  } else { byc <- F
    model <- function(Bt1,Bt,Xt,K,r,p) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))-Xt*((Bt1+Bt)/2)-Bt1
    model_a <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))+Et%*%c*Bt-Xt*((Bt1+Bt)/2)-Bt1
    model_m <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))*exp(Et%*%c)-Xt*((Bt1+Bt)/2)-Bt1
  }

  
  if(plot_out==T){
    
    old_dir <- getwd()
    on.exit( setwd(old_dir))
    
    if (is.null(plot_dir)) plot_dir <- knobi_results$control$plot_settings$plot_dir
    setwd(plot_dir)
    
    if (is.null(plot_filename)) plot_filename <- knobi_results$control$plot_settings$plot_filename
    if ( !plot_filename %in% list.dirs( full.names = FALSE)) dir.create( plot_filename)
    
    setwd(paste0(plot_dir,"/",plot_filename))
    
  }


  # scenarios & years ------------------

  if( byc) Xt <- as.matrix(c) else Xt <- as.matrix(f)

  n_esc <- ncol(Xt)
  if(is.null(colnames(Xt))) colnames(Xt) <- if( n_esc==1){ 'Projection'} else { paste0("Projection_",1:n_esc)}

  sc_names <-colnames(Xt)

  ny <- nrow(Xt)
  nyears <- c( lastyear:c(lastyear+ny))
  pyears <- nyears[-1]
  fyears <- c(years,pyears)
  ly <- length(nyears)


  # init ---------------------

  Bt_ini <- if( knobi_results$control$method == 'SSB') knobi_results$input$SSB else knobi_results$input$Biomass
  B0 <- Bt_ini[length(Bt_ini)]

  params<-knobi_results$params
  r<-params['r']; K<-params['K']; p<-ifelse(pella,params['p'],1)

  base_Bt <- array( c( B0, rep(0,ly-1)), c( ly, n_esc), dimnames = list( nyears+1, sc_names))

  base_Baver <- base_SP <- base_C <- base_F <-
    array(rep(0,ly-1),c(ly-1,n_esc), dimnames = list( pyears,sc_names))

  df <- pdf0 <- subset( knobi_results$df, Year >= lastyear-5)
  df$y <- df$base
  df$Sc <- 'input'
  df$Model <- 'base'
  df$base <- NULL
  
  names(df)[which(names(df)=='F')] <- 'FM'

  
  # base loop ------------------

  for(j in sc_names){

    for(i in 1:ny){

      Bi <- base_Bt[i,j]

      if( Bi == 0){ Bt1 <- 0} else {

        Xi <- Xt[i,j]

        v <- NULL

        try( v <- stats::uniroot( model, c(0,2*K), Bt=Bi, Xt=Xi, K=K, r=r, p=p)$root, silent=TRUE)
        Bt1 <- ifelse( is.null(v), 0, v)

        if( byc){ base_C[i,j] <- Xi; base_F[i,j] <- Xi/((Bi+Bt1)/2)} else {
          base_F[i,j] <- Xi; base_C[i,j] <- Xi*((Bi+Bt1)/2)}}

      if(Bt1<=1e-10){

        base_Bt[c(i+1),j] <- 0
        base_Baver[i,j] <- Bi/2
        base_SP[i,j] <- (r/p)*((Bi)/2)*(1-((Bi)^p)/(K^p*2^p))
        base_C[i,j] <- Bi+base_SP[i,j]
        base_F[i,j] <- base_C[i,j]/(ifelse(Bi==0,1,0.5*Bi))

      } else {

        base_Bt[c(i+1),j]<-Bt1
        base_SP[i,j]<-as.numeric(Bt1-Bi+base_C[i,j])
        base_Baver[i,j]<-(Bt1+Bi)/2

      }
    }

    df <- rbind( df, data.frame( x = base_Baver[,j], y = base_SP[,j], Year = pyears, C = base_C[,j],
                       FM = base_F[,j], Sc = j, Model = 'base'))

  }

  if(any(base_Bt==0)) for(i in which(base_Bt[nrow(base_Bt),]==0))
      warning(paste0('Introduced catch or F in "',sc_names[i],'" scenario lead to stock collapse'))


  # environmental -------------

  if(is.null(env_results)==FALSE & is.null(env)==TRUE) {stop('Environmental data is required')}
  if(is.null(env_results)==TRUE & is.null(env)==FALSE) {stop('Environmental fit results are required')}

  if(is.null(env_results)==FALSE){

    df$EnvSc <- NA

    edf <- subset( env_results$df, Year >= lastyear-5)
    
    names(edf)[which(names(edf)=='F')] <- 'FM'

    for(mi in c('mult','add'))
      df <- rbind( df, data.frame( x = edf$x, y = edf[[mi]], Year = edf$Year, C = edf$C,
                                      FM = edf$FM, Sc = 'input', Model = mi, EnvSc = NA))


    env <- env

    selvar <- env_results$selected_var
    basevar <- env_results$scaled_var

    nvar <- length(selvar)
    csn <- if(nvar>1) paste0('c',1:nvar) else 'c'

    add <- env_results$add$params
    r_a <- add['r']; K_a <- add['K']; p_a <- ifelse(pella,add['p'],1)
    c_a <- add[csn]

    mult <- env_results$mult$params
    r_m <- mult['r']; K_m <- mult['K']; p_m <- ifelse(pella,mult['p'],1)
    c_m <- mult[csn]

    multicovar <- ifelse( length(env_results$selected_var)>1, T, F)

    if( multicovar){

      nsc <- length(env)
      scnames <- if( is.null(names(env))) paste0( 'Environmental_',1:nsc) else names(env)
      varnames <- colnames( env[[1]]); nvars <- ncol(env[[1]])
      if( any(selvar != varnames)) stop( 'Different environmental variables are provided')
      Et <-  array( NA, dim=c(dim(env[[1]]),nsc), dimnames=(list(NULL,varnames,scnames)))
      for(i in 1:nsc) for(j in 1:nvars)
        Et[,j,i] <- (env[[i]][,j] - attr(basevar,"scaled:center")[j])/attr(basevar,"scaled:scale")[j]

    } else {

      env <- as.matrix(env)
      nsc <- ncol(env)
      scnames <- if( is.null(colnames(env))) paste0( 'Environmental_',1:nsc) else colnames(env)
      Et <-  array( NA, dim=c( nrow(env),1,nsc), dimnames=(list(NULL,NULL,scnames)))
      for(i in 1:nsc) Et[,,i] <- (env[,i] - attr(basevar,"scaled:center"))/attr(basevar,"scaled:scale")

    }

    add_Bt <- mult_Bt <- array( c( B0, rep(0,ly-1)), c( ly, n_esc, nsc), dimnames = list( nyears+1, sc_names, scnames))

    add_Baver <- add_SP <- add_C <- add_F <- mult_Baver <- mult_SP <- mult_C <- mult_F <-
      array( 0, c( ly-1, n_esc, nsc), dimnames = list( nyears[-1], sc_names, scnames))


    for(j in sc_names){

      for(n in scnames){

        for(i in 1:ny){

          Bi_a <- add_Bt[i,j,n]
          Bi_m <- mult_Bt[i,j,n]
          Xi <- Xt[i,j]
          Ei <- Et[i,,n]

          if( Bi_a == 0){ Bt1 <- 0} else {

            va <- NULL

            try( va <- stats::uniroot( model_a, c(0,2*(K_a+Ei%*%c_a)),
                              Bt=Bi_a, Xt=Xi, K=K_a, r=r_a, p=p_a, c=c_a, Et=Ei)$root, silent=TRUE)

            Bt1 <- if(is.null(va)) 0 else va

            if( byc){ add_C[i,j,n] <- Xi; add_F[i,j,n] <- Xi/((Bi_a+Bt1)/2)} else {
              add_F[i,j,n] <- Xi; add_C[i,j,n] <- Xi*((Bi_a+Bt1)/2)}
            }

          if( Bt1 <= 1e-10){

            add_Bt[c(i+1),j,n] <- 0
            add_SP[i,j,n] <- (r_a/p_a)*((Bi_a)/2)*(1-((Bi_a)^p_a)/(K_a^p_a*2^p_a))
            add_Baver[i,j,n] <- Bi_a/2
            add_C[i,j,n] <- Bi_a+add_SP[i,j,n]+Ei%*%c_a*Bi_a
            add_F[i,j,n] <- add_C[i,j,n]/(ifelse(Bi_a==0,1,0.5*Bi_a))

          } else {

            add_Bt[c(i+1),j,n] <- Bt1
            add_SP[i,j,n] <- as.numeric(Bt1-Bi_a+add_C[i,j,n])
            add_Baver[i,j,n] <- (Bt1+Bi_a)/2

          }


          if( Bi_m == 0){ Bt1 <- 0} else {

            vm <- NULL

            try( vm <- stats::uniroot( model_m, c(0,2*(K_a+Ei%*%c_a)),
                              Bt=Bi_m, Xt=Xi, K=K_m, r=r_m, p=p_m, c=c_m, Et=Ei)$root, silent=TRUE)

            Bt1 <- if(is.null(vm)) 0 else vm

            if( byc){ mult_C[i,j,n] <- Xi; mult_F[i,j,n] <- Xi/((Bi_m+Bt1)/2)} else {
              mult_F[i,j,n] <- Xi; mult_C[i,j,n] <- Xi*((Bi_m+Bt1)/2)}
            }

          if( Bt1 <= 1e-10){

            mult_Bt[c(i+1),j,n] <- 0
            mult_SP[i,j,n] <- (r_m/p_m)*((Bi_m)/2)*(1-((Bi_m)^p_m)/(K_m^p_m*2^p_m))
            mult_Baver[i,j,n] <- Bi_m/2
            mult_C[i,j,n] <- Bi_m+mult_SP[i,j,n]*exp(Ei%*%c_m)
            mult_F[i,j,n] <- mult_C[i,j,n]/(ifelse(Bi_m==0,1,0.5*Bi_m))

          } else {

            mult_Bt[c(i+1),j,n] <- Bt1
            mult_SP[i,j,n] <- as.numeric(Bt1-Bi_m+mult_C[i,j,n])
            mult_Baver[i,j,n] <- (Bt1+Bi_m)/2

          }



        }

        df <- rbind( df, data.frame( x = add_Baver[,j,n], y = add_SP[,j,n], Year = pyears, C = add_C[,j,n],
                           FM = add_F[,j,n], Sc = j, Model = 'add', EnvSc = n),
                     data.frame( x = mult_Baver[,j,n], y = mult_SP[,j,n], Year = pyears, C = mult_C[,j,n],
                           FM = mult_F[,j,n], Sc = j, Model = 'mult', EnvSc = n))

      }

  }

  }

  ndf <- df

  colnames(df)[which(colnames(df) == 'x')] <- ifelse( knobi_results$control$method=='SSB','SSB','Biomass')
  colnames(ndf)[which(colnames(ndf) == 'x')] <- 'B'
  colnames(df)[which(colnames(df) == 'y')] <- colnames(ndf)[which(colnames(ndf) == 'y')] <- 'SP'

  dummyb <- subset( ndf, Year == nyears[1] & Model =='base')

  if(is.null(env_results)==FALSE){
    dummya <- subset( ndf, Year == nyears[1] & Model =='add')
    dummym <- subset( ndf, Year == nyears[1] & Model =='mult')}

  for(i in sc_names){
    idummyb <- dummyb
    idummyb$Sc <- i
    ndf <- rbind( ndf, idummyb)
  }

  if(is.null(env_results)==FALSE){ for(i in sc_names){ for(j in scnames){

      idummya <- dummya; idummym <- dummym
      idummya$Sc <- idummym$Sc <-i; idummya$EnvSc <- idummym$EnvSc <-j
      ndf <- rbind( ndf, idummya, idummym)

    }}}

  ndf <- subset( ndf, Year >= lastyear-5)

  baseplots <- list()

  bdf <- subset( ndf, Model != 'add' & Model != 'mult')

  bdfl <- tidyr::pivot_longer( bdf, cols = c( B, SP, C, FM),
                               names_to = "Variable",
                               values_to = "Value")

  bdfl$Variable[which(bdfl$Variable=='B')] <- if(knobi_results$control$method=='SSB') 'SSB (t)' else 'B (t)'
  bdfl$Variable[which(bdfl$Variable=='C')] <- 'Catch (t)'
  bdfl$Variable[which(bdfl$Variable=='SP')] <- 'SP (t)'
  bdfl$Variable[which(bdfl$Variable=='FM')] <- 'F'

  baseplots <- ggplot2::ggplot( bdfl, ggplot2::aes( x = Year, y = Value, linetype = Model, color = Sc)) +
    ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
    ggplot2::facet_grid(rows=dplyr::vars(Variable), scales='free_y') +
    ggplot2::labs(title = "Base KBPM projections", y='', color = "Scenario") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  print(baseplots)

  if (plot_out==TRUE){
    p  <-  grDevices::recordPlot()
    grDevices::jpeg(paste0('base_proj.jpeg'),width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  forecast <- list( df=df, plots = list( base = baseplots))

  rquants <- c('B_aver','Catch','FM','SP')

  if(is.null(env_results)==FALSE){

    plotsbyenv <- list()

    for( j in scnames){

      jdf <- subset( ndf, EnvSc == j | Sc == 'input' & Model != 'input')

      jdfl <- tidyr::pivot_longer( jdf, cols = c( B, SP, C, FM),
                                   names_to = "Variable",
                                   values_to = "Value")

      jdfl$Variable[which(jdfl$Variable=='B')] <- if(knobi_results$control$method=='SSB') 'SSB (t)' else 'B (t)'
      jdfl$Variable[which(jdfl$Variable=='C')] <- 'Catch (t)'
      jdfl$Variable[which(jdfl$Variable=='SP')] <- 'SP (t)'
      jdfl$Variable[which(jdfl$Variable=='FM')] <- 'F'

      plotsbyenv[[j]] <- ggplot2::ggplot( jdfl, ggplot2::aes( x = Year, y = Value, linetype = Model, color = Sc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = paste0( '"', j, '" scenario projections'), x = "Year", y = "", color = "Scenario") +
        ggplot2::facet_grid(rows=dplyr::vars(Variable), scales='free_y') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      print(plotsbyenv[[j]])

      if (plot_out==TRUE){
        p  <-  grDevices::recordPlot()
        grDevices::jpeg(paste0(j,'_proj.jpeg'),width=2500, height=2500,res=300)
        grDevices::replayPlot(p)
        grDevices::dev.off()}

    }


    plotsbyC <- list()

    for( j in sc_names){

      jdf <- subset( ndf, Sc == j | Sc == 'input' & Model != 'input')

      jdf$EnvSc[which(is.na(jdf$EnvSc))] <- 'input'

      jdfl <- tidyr::pivot_longer( jdf, cols = c( B, SP, C, FM),
                                   names_to = "Variable",
                                   values_to = "Value")

      jdfl$Variable[which(jdfl$Variable=='B')] <- if(knobi_results$control$method=='SSB') 'SSB (t)' else 'B (t)'
      jdfl$Variable[which(jdfl$Variable=='C')] <- 'Catch (t)'
      jdfl$Variable[which(jdfl$Variable=='SP')] <- 'SP (t)'
      jdfl$Variable[which(jdfl$Variable=='FM')] <- 'F'

      plotsbyenv[[j]] <- ggplot2::ggplot( jdfl, ggplot2::aes( x = Year, y = Value, linetype = Model, color = EnvSc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = paste0( '"', j, '" scenario projections'), x = "Year", y = "", color = "Scenario") +
        ggplot2::facet_grid(rows=dplyr::vars(Variable), scales='free_y') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5), 
                       legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       axis.line=ggplot2::element_line())

      print(plotsbyenv[[j]])

      if (plot_out==TRUE){
        p  <-  grDevices::recordPlot()
        grDevices::jpeg(paste0(j,'_proj.jpeg'),width=2500, height=2500,res=300)
        grDevices::replayPlot(p)
        grDevices::dev.off()}

    }

    forecast$plots[['env']] <- list( byC= plotsbyC, byEnv=plotsbyenv)

  }

  names(forecast$df)[which(names(forecast$df)=='FM')] <- 'F'
  
  class(forecast) <- 'knobi'

  return(forecast)

}

