utils::globalVariables(c("Year", "Residuals", "Model"))

kbpm_error <- function( knobi_results, env_results=NULL, plot_out) {

  X <- knobi_results$df$x
  Y <- knobi_results$df$y
  total_n <- length(X)
  pella <- knobi_results$control$pella

  r_b <- knobi_results$params['r']
  K_b <- knobi_results$params['K']
  p_b <- ifelse( pella, knobi_results$params['p'], 1)

  n_params <- length( knobi_results$params)
  df_b <- total_n-n_params

  rss_b <- sum((Y-((r_b/p_b)*X*(1-(X/K_b)^(p_b))))^2)
  tss_b <- sum((Y-mean(Y))^2)

  ser_b <- sqrt(rss_b/(df_b))
  r2_b <- 1-rss_b/tss_b
  r2_adj_b <- 1-(rss_b/(df_b))/(tss_b/(total_n-1))
  AIC_b <- total_n*log(rss_b/total_n)+2*n_params
  RMSE_b <- sqrt(rss_b/total_n)
  MAPE_b <- 1/total_n*(sum(abs(((Y-((r_b/p_b)*X*(1-(X/K_b)^(p_b)))))/Y)))

  residuals_b <- (Y-((r_b/p_b)*X*(1-(X/K_b)^(p_b))))/stats::sd(Y)

  if(is.null(env_results)){

    error_table <- data.frame(ser_b,r2_b,r2_adj_b,AIC_b, RMSE_b, MAPE_b)
    names <- names(error_table) <- c("SER","R-squared", "adj-R-squared", "AIC", "RMSE", "MAPE")
    rownames(error_table) <- "KBPM_error"

    residuals <- residuals_b

    res_df <- data.frame(Residuals=residuals_b,Year=knobi_results$input$years)

    res_plot <- ggplot2::ggplot(data=res_df,ggplot2::aes(x=Year,y=Residuals)) + ggplot2::theme_bw() +
      ggplot2::geom_hline(yintercept=0) + ggplot2::geom_point() +
      ggplot2::ylim(min(res_df$Residuals),max(res_df$Residuals)) +
      ggplot2::labs(title="Pearson Residuals", subtitle=knobi_results$input$Stock) +
      ggplot2::theme(legend.position = c(0.15,0.85), plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

  } else {

    base_model<-c(ser_b,r2_b,r2_adj_b,AIC_b, RMSE_b, MAPE_b, NA, NA)

    kbpm_add<-env_results$add$params
    kbpm_mult<-env_results$mult$params

    env<-as.matrix(env_results$scaled_var)

    r_a<-kbpm_add['r']; K_a<-kbpm_add['K']
    p_a<-ifelse( pella, kbpm_add['p'], 1)
    c_a<-kbpm_add[(length(kbpm_add)-ncol(env)+1):(length(kbpm_add))]

    n_params_env<-length(kbpm_add)
    df_env<-total_n-n_params_env

    residuals_a <- (Y-as.numeric((r_a/p_a)*X*(1-(X/K_a)^(p_a))+as.numeric(env %*% c_a)*X))/stats::sd(Y)

    rss_a <- sum((Y-as.numeric((r_a/p_a)*X*(1-(X/K_a)^(p_a))+as.numeric(env %*% c_a)*X))^2)
    tss_a <- sum((Y-mean(Y))^2)

    ser_a<-sqrt(rss_a/(df_env))
    r2_a <- 1-rss_a/tss_a
    r2_adj_a <- 1-(rss_a/(df_env))/(tss_a/(total_n-1))
    AIC_a <- total_n*log(rss_a/total_n)+2*n_params_env
    RMSE_a<-sqrt(rss_a/total_n)
    MAPE_a<-1/total_n*(sum(abs((Y-as.numeric((r_a/p_a)*X*(1-(X/K_a)^(p_a))+as.numeric(env %*% c_a)+X))/Y)))
    F_test_a<-((rss_b-rss_a)/(df_b-df_env))/(rss_b/(df_b))
    pF_a<-stats::pf(F_test_a,df_b-df_env,df_env,lower.tail=F)

    env_a_model<-c(ser_a,r2_a,r2_adj_a,AIC_a, RMSE_a, MAPE_a,F_test_a,pF_a)


    r_m<-kbpm_mult['r']; K_m<-kbpm_mult['K']
    p_m<-ifelse( pella, kbpm_mult['p'], 1)
    c_m<-kbpm_mult[(length(kbpm_add)-ncol(env)+1):(length(kbpm_add))]

    residuals_m <- (Y-as.numeric(exp(1)^{as.numeric(env %*% c_m)}*((r_m/p_m)*X*(1-(X/K_m)^(p_m)))))/stats::sd(Y)

    rss_m <- sum((Y-as.numeric(exp(1)^{as.numeric(env %*% c_m)}*((r_m/p_m)*X*(1-(X/K_m)^(p_m)))))^2)
    tss_m <- sum((Y-mean(Y))^2)

    ser_m<-sqrt(rss_m/(df_env))
    r2_m <- 1-rss_m/tss_m
    r2_adj_m <- 1-(rss_m/(df_env))/(tss_m/(total_n-1))
    AIC_m <- total_n*log(rss_m/total_n)+2*n_params_env
    RMSE_m<-sqrt(rss_m/total_n)
    MAPE_m<-1/total_n*(sum(abs((Y-as.numeric(exp(1)^{as.numeric(env %*% c_m)}*((r_m/p_m)*X*(1-(X/K_m)^(p_m)))))/Y)))
    F_test_m<-((rss_b-rss_m)/(df_b-df_env))/(rss_b/(df_b))
    pF_m<-stats::pf(F_test_m,df_b-df_env,df_env,lower.tail=F)

    env_m_model<-c(ser_m,r2_m,r2_adj_m,AIC_m, RMSE_m,
                  MAPE_m,F_test_m,pF_m)

    error_table<-rbind(base_model, env_a_model, env_m_model)
    rownames(error_table)<-c("base model", "additive model", "multiplicative model")
    colnames(error_table)<-c("SER","R-squared", "adj-R-squared", "AIC", "RMSE",
                            "MAPE"," F-value", "Pr(>F)")

    residuals<-list(base=residuals_b,add=residuals_a,mult=residuals_m)

    res_df<-data.frame(Residuals=c(residuals_b,residuals_a,residuals_m),
                      Model=c(rep("Base KBPM",length(residuals_b)),
                              rep("Environmental Additive",length(residuals_a)),
                              rep("Environmental Multiplicative",length(residuals_m))),
                      Year=rep(knobi_results$df$Year,3))

    res_plot<-ggplot2::ggplot(data=res_df,ggplot2::aes(x=Year,y=Residuals,color=Model)) + ggplot2::theme_bw() +
      ggplot2::geom_hline(yintercept=0) + ggplot2::geom_point() +
      ggplot2::ylim(min(res_df$Residuals),max(res_df$Residuals)) +
      ggplot2::labs(title="Pearson Residuals", subtitle=knobi_results$input$Stock) +
      ggplot2::theme(legend.position = c(0.15,0.85), plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

  }
  
  error<-list(error_table=error_table,residuals=residuals)
  
  print(res_plot)
  
  if (plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg(ifelse(is.null(env_results),"residuals.jpeg","residuals_env.jpeg"),width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  return(error)

}



