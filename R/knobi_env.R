#' @title KBPM environmental analysis
#'
#' @name knobi_env
#'
#' @description Analyze and model the relationships between surplus production (SP) and environmental covariable(s) to test whether productivity changes in response to environmental fluctuations. The analysis is conducted in three steps: 
#' 
#' (1) Correlation Analysis: Assess the correlation between the standardized environmental variable(s), at different delays (lags), and the KBPM residuals using Pearson's correlation or autoregressive models. See details.
#' 
#' (2) Variable Selection: Determine which lagged environmental variable(s) will be included in the environmental KBPM models.
#' 
#' (3) Environmental Fitting: Fit the KBPM model incorporating environmental effects as both additive and multiplicative effects (see details).
#'
#' @param knobi_results The output object of \code{\link{knobi_fit}} function (main package function).
#' @param data A list containing the following data: \itemize{
#'   \item \code{env}: containing the values of each environmental variable, with each column representing a different variable and each row representing a year.
#'   \item \code{years}: years in which the environmental variable(s) are reported.}
#' @param control Optional. List containing the following settings: \itemize{
#'   \item \code{nlag}: this argument specifies the maximum lag of the environmental variable to test in the correlation analysis, meaning that lags less than or equal to \code{nlag} (a natural number) are evaluated. The correlation between \eqn{residuals_t} and \eqn{X_{t - lag}} is computed, where \eqn{X} is the environmental variable and \eqn{lag} takes values from the sequence \eqn{\{0, 1, \ldots, nlag\}}. The lagged environmental variable corresponding to the highest correlation with the KBPM residuals is included in the environmental model. By default, \code{nlag = 3}. See details.
#'   \item \code{lag}: an optional numerical vector specifying the lag value(s) to consider in the relationship between the KBPM surplus production and the environmental variable(s). The length of this vector must match the number of environmental variables included. This argument applies only if the \code{nlag} argument is not provided.
#'   \item \code{start_c}: optional. A numerical vector specifying the starting values for the environmental parameter \eqn{c} in the additive and multiplicative models, respectively. By default, \code{start_c = c(1, 1)}. See details.
#'   \item \code{ar_cor}: optional. Logical. By default, this argument is \code{FALSE}, meaning the correlation between the KBPM residuals and the environmental variable(s) is analyzed using the Pearson correlation measure. If set to \code{TRUE}, the relationship is instead analyzed by fitting autoregressive (AR) models, with each lagged environmental variable included as an explanatory covariate for KBPM residuals. The environmental variable associated with the model that has the lowest Akaike Information Criterion (AIC) is selected for inclusion in the environmental KBPM fit. See details.
#'   \item \code{plot3d}: optional. Logical. If set to \code{TRUE}, 3D plots are generated, displaying the surplus production curve across a range of values for the environmental variable(s). The default is \code{FALSE}.
#'   \item \code{multicovar}: optional. Logical. If \code{TRUE}, the environmental model incorporates all input environmental covariates simultaneously. By default, this argument is \code{FALSE}, meaning that only the environmental variable with the highest correlation (after lagging, if applicable) is included in the model.
#' }
#' @param plot_out Logical. If set to \code{TRUE}, a file containing the plot of the environmental fits is created. The default value is taken from the input in the \code{\link{knobi_fit}} function.
#' @param plot_dir Optional. Directory where the folder for saving the plots will be created. Required when \code{plot_out = TRUE}. The default value is taken from the input in the \code{\link{knobi_fit}} function.
#' @param plot_filename Optional. Name of the folder that will contain the plots. Required when \code{plot_out = TRUE}. The default value is taken from the input in the \code{\link{knobi_fit}} function.
#'
#' @details It is important to mention that the environmental variable(s), in a first step, are standardized, in order to make their scale and magnitude comparable. To do this, each variable is subtracted from its mean and divided by its standard deviation.
#' 
#' The additive environmental model adds the following term on the right-hand side of Eq. (1) or Eq. (2) described in \code{\link{knobi_fit}}: \eqn{c X_{t - lag} B_t}, being \eqn{X_{t - lag}} the environmental variable at time \eqn{t - lag} and \eqn{B_t} the biomass or SSB at time \eqn{t}.
#' 
#' The multiplicative environmental model multiplies the right-hand side of Eq. (1) or Eq. (2) by \eqn{\exp(c X_{t - lag})}.
#' 
#' In the case of these models, the estimated biological reference points (BRPs) correspond to a value of the scaled environmental variable equal to the mean of the time series, i.e., \eqn{X_t = 0}, which cancels out the effect of the parameter \eqn{c}. The estimates of the remaining parameters included in Eq. (1) or Eq. (2), and therefore for the BRPs as well, will differ from the base model due to the inclusion of environmental effects. For more details, such as the calculation of BRPs as a function of the environmental variable, see the vignettes.
#'
#' If \code{ar_cor = TRUE}, the correlation analysis between the \code{\link{knobi_fit}} residuals and the environmental variable(s) is conducted as follows:
#'
#' First, an AR model is fitted to the KBPM base residuals:
#' \deqn{r_t = \sum_{i = 1}^{\rho} \beta_i r_{t - i} + \epsilon_t}
#' where \eqn{r_t} is the KBPM base residual for year \eqn{t}, and \eqn{\rho} is the AR model order, estimated as the maximum lag at which the absolute value of the residuals' partial autocorrelation exceeds \eqn{qnorm(0.975) / \sqrt{N_r}}, with \eqn{N_r} being the length of the residuals series.
#'
#' AR models are then fitted to the residuals incorporating each lagged environmental variable \eqn{X_{t - lag}} as an explanatory covariate:
#' \deqn{r_t = \sum_{i = 1}^{\rho} \beta_i r_{t - i} + X_{t - lag} + \epsilon_t}
#' for \eqn{lag = 0, 1, \ldots, nlag}. Then, we have one autoregressive model for each lag.
#'
#' The lagged environmental variable whose AR model has the lowest AIC is selected for inclusion in the environmental KBPM fit. If none improve the AIC compared to the base AR model, no environmental covariate is added.
#' 
#' It is important to highlight that the results include the analysis of the AR model only with the base model residuals in order to determine the need for coupling environmental information, considering that it would not be necessary if this model shows a lower AIC, even reducing the number of parameters to fit.
#'
#' @return A list containing the results of the three-step environmental analysis: \itemize{
#'   \item \code{add}: estimates of the additive model parameters.
#'   \item \code{mult}: estimates of the multiplicative model parameters.
#'   \item \code{BRPs}: reference points (BRPs) estimates for each model (see details).
#'   \item \code{df}: data frame with the information used in the fit.
#'   \item \code{selected_var}: environmental variable(s) used in the fit.
#'   \item \code{selected_lag}: data frame providing the time lag of the environmental variable(s) in the KBPM fit.
#'   \item \code{lag_cor}: correlation between the environmental variable(s) and the KBPM residuals for each time lag.
#'   \item \code{env_aic}: if \code{ar_cor = TRUE}, AIC values of each autoregressive model (see details).
#'   \item \code{scaled_var}: standardized environmental variable(s) used in the fit.
#'   \item \code{plots3D}: list with 3D plot objects (if \code{plot3d = TRUE}).
#'   \item \code{residuals}: Pearson residuals from each model fit (base, additive, and multiplicative).
#'   \item \code{performance_metrics}: array of performance and accuracy measures for each model, including those from the \code{error_table} output of \code{\link{knobi_fit}}, plus an F-test comparing environmental models to the base model.
#' }
#'
#' The function also returns plots in the graphics window and saves them (if \code{plot_out = TRUE}) to the specified directory. The first plot shows the correlation analysis between the environmental variable(s) and the base residuals. The second shows fitted SP values from all models. If \code{multicovar = FALSE} and \code{plot3d = TRUE}, 3D plots are also returned. If \code{multicovar = TRUE}, a plot with Pearson correlations among environmental variables is included.
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
#' # First, run the example of knobi_fit function
#'
#' # Then, provide environmental data series
#'
#' Env <- knobi_dataset$Env
#'
#' # The environmental data series must start in the first year of the KBPM fit data
#' # minus the provided nlag or lag arguments
#' years <- knobi_results$df$Year # See knobi_fit example to obtain the knobi_results object
#' ind <- which(Env[,1]==years[1])
#' ind1 <- which(Env[,1]==years[length(years)])
#' nlag <- 5
#' Env <- Env[(ind-nlag):ind1,]
#'
#' # Now we create the environmental list
#' data <- list(env=data.frame(AMO=Env$AMO,NAO=Env$NAO),
#'            years=Env$years)
#' control <- list(nlag=nlag)
#'
#' knobi_environmental <- knobi_env(knobi_results,data,control)
#' knobi_environmental
#' knobi_environmental$BRPs  # use the '$' to access to all the fit information
#' }
#'
#' @export


utils::globalVariables(c("correlation","AIC","SP"))

knobi_env <- function( knobi_results, data, control=NULL, plot_out=FALSE, plot_filename=NULL, plot_dir=NULL){

  if( is.null( control)) { control <- list()}

  if( is.null( control$ar_cor)) { control$ar_cor <- FALSE}

  input <- knobi_results$input
  bdata <- knobi_results$df
  pella <- knobi_results$control$pella

  env0 <- as.data.frame(data$env)
  env_names <- names(env0)
  y_env <- data$years
  res <- knobi_results$residuals

  if( is.null(control$start_c)){ start_c <- c(1,1)} else { start_c <- control$start_c}

  if(plot_out==T){

    old_dir <- getwd()
    on.exit( setwd(old_dir))

    if (is.null(plot_dir)) plot_dir <- knobi_results$control$plot_settings$plot_dir
    setwd(plot_dir)
    
    if (is.null(plot_filename)) plot_filename <- knobi_results$control$plot_settings$plot_filename
    if ( !plot_filename %in% list.dirs( full.names = FALSE)) dir.create( plot_filename)

    setwd(paste0(plot_dir,"/",plot_filename))

  }

  df <- cbind( res, bdata)

  f_year <- df$Year[1]; l_year <- df$Year[length(df$Year)]

  env <- list()
  res_env <- list()

  df_env <- df

  res_env$selected_lag <- array( NA, dim=c(length(env_names),2))
  colnames(res_env$selected_lag) <- c("lag","correlation")
  rownames(res_env$selected_lag) <- env_names

  if ( is.null(control$lag)){ lag <- ifelse( is.null(control$nlag), 3, control$nlag)
  } else { lag <- max(control$lag); res_env$selected_lag[,1] <- control$lag}

  # Env data ----------------------------

  data_env <- list()

  res_env$lag_cor <- array( NA, dim=c(length(env_names), lag+1))

  vec_env <- "lag_0"

  for (i in 1:lag) vec_env <- c(vec_env,paste0("lag_",i))

  colnames(res_env$lag_cor) <- vec_env
  rownames(res_env$lag_cor) <- env_names


  for(j in env_names){

    data_env[[j]] <- df

    ind <- which( y_env == f_year)
    ind1 <- which( y_env == l_year)

    if( length(ind) > 0){ data_env[[j]]$env0 <- env0[ind:ind1,j]} else {
      warning('The length of the environmental variable is not enough to use the number of input lags')}

    for ( i in 1:lag){
      ind <- which( y_env == (f_year-i))
      ind1 <- which( y_env == (l_year-i))
      if(length(ind) > 0){ data_env[[j]][,5+i] <- env0[[j]][ind:ind1]} else {
          warning('The length of the environmental variable is not enough to use the number of input lags')}
    }

    bname <- ifelse( knobi_results$control$method=="Biomass", 'B', 'SSB')

    colnames(data_env[[j]]) <- c("res","SP",bname,"years",vec_env)

    if(!is.na(input$Recruitment[1])) data_env[[j]]$R <- input$Recruitment

    env[[j]] <- data_env[[j]][,c("res",vec_env)]

    cor <- round(cor(as.matrix(env[[j]]),use="na.or.complete"),4)
    cor <- cor[,1]; cor <- cor[-1]

    if (is.null(control$lag)){

      res_env$selected_lag[j,2] <- cor[[which(max(abs(cor))==abs(cor))[1]]]
      res_env$selected_lag[j,1] <- which(max(abs(cor))==abs(cor))[1]-1

    } else {

      res_env$selected_lag[j,2] <- cor[res_env$selected_lag[j,1]+1]

    }

    df_env[,j] <- scale(env[[j]][,res_env$selected_lag[j,1]+2])

    res_env$lag_cor[j,] <- cor

  }


  # AR cor ------------------

  if( control$ar_cor==TRUE){

    pacf_res <- stats::pacf( res, plot=FALSE)$acf[,1,1]

    ref <- stats::qnorm(0.975)/sqrt(length(res))

    auto <- max(which(abs(pacf_res)>=ref),0)
    if(auto==0){ warning("KBPM base residuals are not autocorrelated. An AR(0) model is fitted for SP residuals.")}

    fit_base  <-  stats::arima0( res, order = c(auto,0,0))

    env_aic <- array( NA, dim=c(length(env_names),lag+1))
    colnames(env_aic) <- vec_env; rownames(env_aic) <- env_names

    for(j in env_names) for(i in vec_env) env_aic[j,i] <- stats::arima0( res, order = c(auto,0,0), xreg = env[[j]][,i])$aic

    env_aic <- cbind(base=rep(fit_base$aic,length(env_names)),env_aic)
    res_env$env_aic <- env_aic

    env_aic_c <- env_aic-fit_base$aic+2
    min_aic <- env_aic_c[which.min(env_aic_c)]
    if(min_aic>=0) warning("AR models considering environmental variable(s) do not really improve AR model considering only the residuals")

    colnames(res_env$selected_lag)[2] <- "aic"

    if (is.null(control$lag)){

      for(j in env_names){
        res_env$selected_lag[j,1] <- which(env_aic == min(env_aic[j,-1]), arr.ind=TRUE)[2]-2
        res_env$selected_lag[j,2] <- min(env_aic[j,-1])
      }

    } else {

      for(j in env_names){
        res_env$selected_lag[j,2] <- env_aic[j,(res_env$selected_lag[j,1]+2)]
      }
    }

    for(j in env_names){
      df_env[,j] <- scale(env[[j]][,res_env$selected_lag[j,1]+2])
    }

  }


  # corr plot ----------------

  lagf <- NULL
  corlist <- NULL

  if( control$ar_cor == FALSE){

    for(i in vec_env){
      lagf <- c(lagf,rep(i,length(env_names)))
      corlist <- c(corlist,res_env$lag_cor[,i])
    }

    envcorplot_df <- data.frame(correlation=corlist,lag=lagf,factor=rep(env_names,length(vec_env)))

    envcorplot <- ggplot2::ggplot(data=envcorplot_df,ggplot2::aes(x=lag,y=correlation,group=factor,color=factor)) +
      ggplot2::theme_bw() + ggplot2::geom_point() + ggplot2::geom_line(linetype = "dashed") + ggplot2::ylim(-1,1) +
      ggplot2::labs(title="Correlation between environmental variables and base KBPM Residuals", subtitle=knobi_results$input$Stock,
                    y="Correlation",x="") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

    print(envcorplot)

  } else {

    for(i in vec_env){
      lagf <- c(lagf,rep(i,length(env_names)))
      corlist <- c(corlist,res_env$env_aic[,i])
    }

    maximo <- max(env_aic)+2
    minimo <- min(env_aic)-2

    envcorplot_df <- data.frame(AIC=corlist,lag=lagf,factor=rep(env_names,length(vec_env)))

    envcorplot <- ggplot2::ggplot(data=envcorplot_df,ggplot2::aes(x=lag,y=AIC,group=factor,color=factor)) +
      ggplot2::theme_bw() + ggplot2::geom_point() + ggplot2::geom_line(linetype = "dashed") + ggplot2::ylim(minimo,maximo) +
      ggplot2::labs(title="AIC comparison", subtitle=knobi_results$input$Stock,
                    y="AIC",x="") +
      ggplot2::geom_hline(yintercept = fit_base$aic) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

    print(envcorplot)

  }


  if (plot_out==TRUE){
    p  <-  grDevices::recordPlot()
    grDevices::jpeg("corplot.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  multicovar <- ifelse( is.null(control$multicovar), FALSE, control$multicovar)

  # Fit  ---------------

  if(multicovar == FALSE){

    if( control$ar_cor == FALSE){ selected_var <- env_names[which.max(abs(res_env$selected_lag[,2]))]
      } else { selected_var <- env_names[which.min(abs(res_env$selected_lag[,2]))]}

  } else {

    selected_var <- env_names

    if(!is.null(control$plot3d)) control$plot3d <- NULL

  }


  envsel <- df_env[,selected_var]
  res_env$selected_var <- colnames( envsel) <- selected_var

  if( multicovar){

    cor <- round(cor(as.matrix(envsel),use="na.or.complete"),4)
    for(i in 1:ncol(cor)) cor[i,i] <- 0
    cor <- max(abs(cor))
    if(cor >= 0.35){
      p.mat  <-  corrplot::cor.mtest(envsel)$p
      corrplot::corrplot(round(cor(as.matrix(envsel),use="na.or.complete"),2),method='circle',
                         title="Covariables correlation",mar=c(0,0,2,0),addCoef.col = "black",
                         type="lower",diag=T,p.mat = p.mat, sig.level = 0.05)
      warning("The covariables are highly correlated (see corrplot). To avoid estimation issues, consider removing variables or setting multicovar=FALSE.")
    }

  }

  Data <- list( data=df, start_r=as.numeric(knobi_results$params[['r']]),
    start_K=as.numeric(knobi_results$params[['K']]), start_c=start_c[1], start_p=1)

  kbpm_mult  <- kbpm_fit( Data, envsel, pella, 'Mult')
  res_env$mult$params  <-  kbpm_mult
  res_env$mult$BRPs  <-  BRP( res_env$mult, pella)

  kbpm_add  <- kbpm_fit( Data, envsel, pella, 'Add')
  res_env$add$params  <-  kbpm_add
  res_env$add$BRPs  <-  BRP( res_env$add, pella)

  envdf <- df
  envdf <- cbind( envdf, envsel)

  bv  <- envdf$base
  bv1  <-  envdf$mult <- predict_model( res_env$mult, df$x, pella, 'Mult', envsel)
  bv2  <-  envdf$add <- predict_model( res_env$add, df$x, pella, 'Add', envsel)

  res_env$scaled_var <- envsel

  if( multicovar == F) rownames(res_env$scaled_var) <- df_env$Year-as.numeric(res_env$selected_lag[selected_var,1])

  if(is.null(control$plot3d)){plots3d <- FALSE} else {plots3d <- control$plot3d}


  # 3D plots ---------------------------

  if(plots3d==TRUE){

    x  <-  envdf$x
    y  <-  envdf[[selected_var]]
    ysc <- envsel
    z  <-  envdf$y

    r_a <- kbpm_add['r']; K_a <- kbpm_add['K']; c_a <- kbpm_add['c']
    p_a <- ifelse( pella, kbpm_add['p'], 1)

    r_m <- kbpm_mult['r']; K_m <- kbpm_mult['K']; c_m <- kbpm_mult['c']
    p_m <- ifelse( pella, kbpm_mult['p'], 1)


    grid.lines  <-  400
    cut_a <- max( K_a+c_a*K_a*max(y), K_a+c_a*K_a*min(y))
    x.pred_a  <-  seq(0, cut_a, length.out = grid.lines)
    x.pred_m  <-  seq(0, K_m, length.out = grid.lines)

    y.pred  <-  c( seq(min(y), 0, length.out = grid.lines/2), seq( 0, max(y), length.out = grid.lines/2))
    xy_a  <-  expand.grid(x = x.pred_a, y = y.pred)
    xy_m  <-  expand.grid(x = x.pred_m, y = y.pred)

    z.pred_a  <-  (r_a/p_a) * xy_a$x * (1-(xy_a$x/K_a)^(p_a)) + c_a * xy_a$y * xy_a$x
    z.pred_m  <-  exp(c_m*xy_m$y) * ((r_m/p_m) * xy_m$x * (1-(xy_m$x/K_m)^p_m))
    z.pred_a  <-  matrix(z.pred_a,nrow=grid.lines,ncol=grid.lines)
    z.pred_m  <-  matrix(z.pred_m,nrow=grid.lines,ncol=grid.lines)

    y2  <-  y * attr(ysc, 'scaled:scale') + attr(ysc, 'scaled:center')
    y.pred2  <-  y.pred * attr(ysc, 'scaled:scale') + attr(ysc, 'scaled:center')

    plot3D::scatter3D(x, y2, z, bty="b2", pch = 19, cex = 1, cex.axis = 0.6,
                      colkey = list(length = 0.5, width = 0.5, cex.clab = 0.8, cex.axis = 0.8),
                      xlim = c(0,max(x.pred_a)), zlim = c(0,max(z.pred_a)*1.5),
                      theta = 28, phi = 20, ticktype = "detailed", clim = c(0,max(z.pred_a,z)),
                      xlab = "SSB", ylab = selected_var, zlab = "SP", clab = "Surplus production",
                      surf = list(x = x.pred_a, y = y.pred2, z = z.pred_a, facets = NA, fit = z),
                      main = "Additive model: Production curve", sub = knobi_results$input$Stock)

    plot3d_add <- grDevices::recordPlot()

    if (plot_out == TRUE){
      grDevices::jpeg("additive_model.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(plot3d_add)
      grDevices::dev.off()
    }

    plot3D::scatter3D(x, y2, z, bty="b2", pch = 19, cex = 1, cex.axis = 0.6,
                      colkey = list(length = 0.5, width = 0.5, cex.clab = 0.8, cex.axis = 0.8),
                      xlim = c(0, max(x.pred_m)), zlim = c(0,max(z.pred_m)*1.5),
                      theta = 28, phi = 18, ticktype = "detailed", clim=c(0,max(z.pred_m,z)),
                      xlab = "SSB", ylab = selected_var, zlab = "SP", clab = "Surplus production",
                      surf = list(x = x.pred_m, y = y.pred2, z = z.pred_m, facets = NA, fit = z),
                      main = "Multiplicative model: Production curve", sub = knobi_results$input$Stock)
    plot3d_mult <- grDevices::recordPlot()

    if (plot_out==TRUE){
      grDevices::jpeg("multiplicative_model.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(plot3d_mult)
      grDevices::dev.off()
    }

    res_env$plots3D <- list(additive_plot=plot3d_add,multiplicative_plot=plot3d_mult)

  }

  envplot_df <- data.frame(SP=c(df$y,bv,bv2,bv1),Year=rep(df$Year,4),
                         factor=c(rep("Observed",length(bv)),rep("Base KBPM",length(bv)),
                                  rep("Environmental Additive",length(bv)),
                                  rep("Environmental Multiplicative",length(bv))))

  env_plot <- ggplot2::ggplot(data=envplot_df,ggplot2::aes(x=Year,y=SP,color=factor)) + ggplot2::theme_bw() +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min(envplot_df$SP),max(envplot_df$SP)) +
    ggplot2::labs(title="Environmental fits", subtitle=knobi_results$input$Stock,
                  y="Surplus Production") +
    ggplot2::theme(legend.position = c(0.15,0.85), plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"))

  print(env_plot)

  if (plot_out==TRUE){
    p  <-  grDevices::recordPlot()
    grDevices::jpeg("fits_env.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  errors  <-  kbpm_error( knobi_results, res_env, plot_out)

  res_env$performance_metrics <- errors$error_table
  res_env$residuals <- errors$residuals

  if(multicovar){
    mscale <- mcenter <- 1:ncol(res_env$scaled_var)
    for (i in 1:ncol(res_env$scaled_var)){
      mscale[i] <- attr( df_env[,selected_var[i]], 'scaled:scale')
      mcenter[i] <- attr( df_env[,selected_var[i]], 'scaled:center')}
    attr( res_env$scaled_var, 'scaled:scale') <- mscale
    attr( res_env$scaled_var, 'scaled:center') <- mcenter}

  res_env$input <- list( control = control, data = data, basecontrol = knobi_results$control)
  res_env$base <- list( params = knobi_results$params, BRPs = knobi_results$BRPs)

  res_env$BRPs <- rbind( res_env$base$BRPs, res_env$add$BRPs, res_env$mult$BRPs)
  rownames( res_env$BRPs) <- c( 'Base', 'Add', 'Mult')

  envdf <- envdf

  res_env$df <- envdf

  class( res_env) <- 'knobi'

  return( res_env)

}
