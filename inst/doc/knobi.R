## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(knobi)

## -----------------------------------------------------------------------------
data( knobi_dataset)
hake_n <- knobi_dataset$hake_n

## -----------------------------------------------------------------------------
data <- list(
  SSB = hake_n$SSB, 
  Catch = hake_n$catches,
  F_input = hake_n$F, 
  Recruitment = hake_n$recruitment,
  RP = list( F_MSY = 0.26),           # Provided by ICES        
  years = hake_n$Year )

## -----------------------------------------------------------------------------
control <- list( pella = "TRUE")

## ----eval=FALSE---------------------------------------------------------------
# hake_n_results <- knobi_fit( data = data,
#                              control = control,
#                              plot_out = FALSE)

## ----echo=FALSE, fig.show='hide'----------------------------------------------
hake_n_results <- knobi_fit( data, control, plot_out=FALSE)

## ----echo=FALSE, warning=FALSE,  fig.width=6, fig.height=4, fig.align = 'center'----
hake_n_results <- knobi_fit( data, control, plot_out=FALSE)

## -----------------------------------------------------------------------------
hake_n_results

## -----------------------------------------------------------------------------
hake_n_results$BRPs

## -----------------------------------------------------------------------------
Env <- knobi_dataset$Env
nlag <- 5
years <- hake_n_results$df$Year

ind <- which(Env[,1]==years[1])
ind1 <- which(Env[,1]==years[length(years)])

Env <- Env[(ind-nlag):ind1,]

## -----------------------------------------------------------------------------
data <- list(
  env = data.frame( AMO=Env$AMO, NAO=Env$NAO),
  years = Env$years)

## -----------------------------------------------------------------------------
control <- list( nlag = nlag)

## ----eval=FALSE---------------------------------------------------------------
# hake_n_environmental <- knobi_env(knobi_results = hake_n_results,
#                                   data = data,
#                                   control = control,
#                                   plot_out = FALSE)

## ----echo=FALSE,  fig.width=6, fig.height=4, fig.align = 'center'-------------
hake_n_environmental <- knobi_env(hake_n_results,data,control)

## -----------------------------------------------------------------------------
hake_n_environmental

## -----------------------------------------------------------------------------
hake_n_environmental$BRPs

## ----eval=FALSE---------------------------------------------------------------
# control$plot3d = TRUE
# knobi_env( hake_n_results, data, control)

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
control <- list( lag=c(2,3), multicovar=TRUE)
hake_n_multi <- knobi_env( hake_n_results, data, control)

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
control_ar <- list( nlag=3, ar_cor=TRUE)
hake_env_ar <- knobi_env( hake_n_results, data = data, control = control_ar)

## -----------------------------------------------------------------------------
hake_env_ar$env_aic

## -----------------------------------------------------------------------------
hake_env_ar$selected_lag

## ----eval=FALSE---------------------------------------------------------------
# hake_n_retros <- knobi_retro( knobi_results = hake_n_results,
#                               nR = 5,
#                               plot_out = FALSE)

## ----echo=FALSE,  fig.width=6, fig.height=4, fig.align = 'center'-------------
hake_n_retros <- knobi_retro( hake_n_results, nR=5, plot_out=FALSE)

## -----------------------------------------------------------------------------
hake_n_retros

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
knobi_retro( hake_n_results, 
             yR = c(2005,2010,2015))

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
knobi_retro( hake_n_results,
             yR = c(2005,2010,2015),
             yR0 = c(1990,1995,1995))

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
knobi_retro(hake_n_results, hake_n_environmental, nR = 3); hake_n_retros

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
knobi_retro( hake_n_results, hake_n_multi,
             yR = c(2005,2010,2015),
             yR0 = c(1990,1995,1995))

## -----------------------------------------------------------------------------
catch <- rep(hake_n_results$input$Catch[length(hake_n_results$input$Catch)],8)
C <- data.frame(catch=catch, catch08=0.8*catch, catch12=1.2*catch)

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
projections <- knobi_proj( knobi_results=hake_n_results, c=C)

## -----------------------------------------------------------------------------
projections

## -----------------------------------------------------------------------------
last_AMO <- Env$AMO[length(Env$AMO)]
env <- data.frame( AMOi=rep(last_AMO,5),
                   AMOii=rep(last_AMO*1.5,5),
                   AMOiii=rep(last_AMO*0.5,5))

C <- C[(1:5),]

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
env_projections <- knobi_proj(hake_n_results, hake_n_environmental, c=C, env=env)
env_projections

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
fmsy <- hake_n_results$BRPs['F_MSY']
ff <- rep(fmsy,5)
f <- data.frame( f=ff, f12=ff*1.2, f08=ff*0.8)

f_projections <- knobi_proj( hake_n_results, f=f, env_results=hake_n_environmental, env=env)
f_projections

## ----fig.width=6, fig.height=4, fig.align = 'center'--------------------------
env <- list( climate_1 = data.frame( AMO=c(0.2,0.2,0.3,0.3,0.4),
                                     NAO=c(0.2,0.2,0.3,0.3,0.4)),
             climate_2 = data.frame( AMO=c(0.2,0.3,0.4,0.5,0.6),
                                     NAO=c(0.2,0.3,0.4,0.5,0.6)))

multiproj <- knobi_proj( hake_n_results, hake_n_multi, c=C, env=env)
multiproj

