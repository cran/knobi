kbpm_fit <- function(Data, env = NULL, pella = TRUE, model = "Base") {

  data <- Data$data
  start_r <- Data$start_r
  start_K <- Data$start_K
  start_p <- if (pella) Data$start_p else NULL

  n_env <- ifelse(model != 'Base', ncol(env), 0)

  params <- c(r = start_r, K = start_K)
  uppers <- c(2, Inf)
  lowers <- c(0.05, 0)

  if (pella) params <- c(params, p = start_p); uppers <- c(uppers, 3.5); lowers <- c(lowers, 0.25)

  my_model <- switch( paste0(model, "_", ifelse(pella, "pella", "no_pella")),
    "Base_pella" = function(x, r, K, p) (r / p) * x * (1 - (x / K)^p),
    "Base_no_pella" = function(x, r, K) r * x * (1 - (x / K)),
    "Mult_pella" = function(x, env, r, K, p, c) (r / p) * x * (1 - (x / K)^p) * exp(env),
    "Mult_no_pella" = function(x, env, r, K, c) r * x * (1 - (x / K)) * exp(env),
    "Add_pella" = function(x, env, r, K, p, c) (r / p) * x * (1 - (x / K)^p) + env * x,
    "Add_no_pella" = function(x, env, r, K, c) r * x * (1 - (x / K)) + env * x)

  if (n_env > 0) {
    start_c <- Data$start_c
    params <- c(params, rep(start_c, n_env))
    uppers <- c(uppers, rep(20, n_env))
    lowers <- c(lowers, rep(-20, n_env))
    cname <- if (n_env > 1) paste0("c", 1:n_env) else "c"
    names(params)[(length(params) - n_env + 1):length(params)] <- cname
  }

  error_fun <- if (model == 'Base'){
    if (pella) { function(p, x, y) sum((y - my_model(x, p["r"], p["K"], p["p"]))^2)
      } else { function(p, x, y) sum((y - my_model(x, p["r"], p["K"]))^2)}} else {
    if (pella) { function(p, x, y, env){
        env <- if (n_env > 1) as.matrix(env) %*% p[4:length(p)] else p["c"] * env
        sum((y - my_model(x, env, p["r"], p["K"], p["p"], p[-(1:4)]))^2)}
    } else { function(p, x, y, env){
        env <- if (n_env > 1) as.matrix(env) %*% p[3:length(p)] else p["c"] * env
        sum((y - my_model(x, env, p["r"], p["K"], p[-(1:3)]))^2)}}
  }

  if (model == "Base") { out <- optimx::optimr(params, error_fun, x = data$x, y = data$y, method = "L-BFGS-B", upper = uppers, lower = lowers)
  } else { out <- optimx::optimr(params, error_fun, x = data$x, y = data$y, method = "L-BFGS-B", upper = uppers, lower = lowers, env = env)}

  if (model != "Base") {
    out <- out$par
    attr(out, "status") <- NULL
  }

  return(out)
}
