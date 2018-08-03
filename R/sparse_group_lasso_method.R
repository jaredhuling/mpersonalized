sparse_group_lasso_method = function(modelYlist, modelXlist, Ybarlist,
                                     Xbarlist, Xsdlist, lambda, alpha)
  {

  q = length(modelXlist)
  p = dim(modelXlist[[1]])[2]

  x = as.matrix(bdiag(modelXlist))
  y = unlist(modelYlist)
  total_n = length(y)

  x = sqrt(total_n) * x; y = sqrt(total_n) * y

  data = list(x = x, y = y)
  SGL_model = SGL(data = data, index = rep(1 : p, q),
                  lambdas = lambda, alpha = alpha, standardize = FALSE)

  nlambda = length(lambda)
  betalist = vector("list", nlambda)
  interceptlist = vector("list", nlambda)

    for (ind in 1:nlambda){

    beta = matrix(SGL_model$beta[, ind], nrow = q, ncol = p, byrow = TRUE)
    intercept = unlist(Ybarlist) - apply(beta * matrix(unlist(Xbarlist), nrow = q, ncol = p, byrow = TRUE), 1, sum)
    beta = beta / matrix(unlist(Xsdlist), nrow = q, ncol = p, byrow = TRUE)

    betalist[[ind]] = beta
    interceptlist[[ind]] = intercept
  }

  return(list(interceptlist = interceptlist, betalist = betalist,
              lambda = lambda))
}


make_design_matrix_sgl_fused <- function(modelXlist, tau0)
{
  ss.vec     <- unlist(lapply(modelXlist, nrow))
  ss.tot     <- sum(ss.vec)
  tau_s      <- sqrt(ss.tot / ss.vec) * tau0
  studies    <- length(modelXlist)
  weight.vec <- 1 / unlist(lapply( 1:studies, function(s) rep(tau_s[s], ss.vec[s]) ))

  list(design = cbind(do.call(rbind, modelXlist), weight.vec * as.matrix(bdiag(modelXlist)) ),
       tau_s = tau_s)
}

recover_beta_sgl_fused <- function(coefs, p, studies, tau_s)
{
  beta.global <- coefs[1:p]
  beta.locals <- coefs[-(1:p)]

  beta <- numeric(p * studies)

  for (s in 1:studies)
  {
    idx.cur <- ((s - 1) * p + 1):(s * p)
    beta[idx.cur] <- beta.locals[idx.cur] / tau_s[s] + beta.global
  }
  beta
}

gen_tau0 <- function(ntau, min.tau = 1e-2)
{
  # generate values for tau0 on a log-linear scale
  # when 0 < tau < 1
  # and on an exp-linear scale
  # when tau > 1
  if (ntau %% 2 == 0)
  {
    tau1 <- exp(seq(log(min.tau), log(1), length.out = ntau / 2 ))
    tau2 <- exp(seq(log(1), log(1 / min.tau), length.out = ntau / 2 + 1 ))

    tau  <- c(tau1, tau2[-1])
  } else
  {
    tau1 <- exp(seq(log(min.tau), log(1), length.out = (ntau + 1) / 2 ))
    tau2 <- exp(seq(log(1), log(1 / min.tau), length.out = (ntau + 1) / 2 ))

    tau  <- c(tau1, tau2[-1])
  }
  tau
}

sparse_group_fused_lasso_method = function(modelYlist, modelXlist, Ybarlist,
                                           Xbarlist, Xsdlist, lambda = NULL, alpha,
                                           fused = FALSE, tau0 = 1, nlambda = 50)
{

  q = length(modelXlist)
  p = dim(modelXlist[[1]])[2]

  y = unlist(modelYlist)
  total_n = length(y)
  y       = sqrt(total_n) * y

  if (!is.null(lambda))
  {
    nlambda = length(lambda)
  }

  ntau    = length(tau0)

  betalist      = vector("list", nlambda * ntau)
  interceptlist = vector("list", nlambda * ntau)

  penalty_parameter_sequence = matrix(0, ncol = 2, nrow = nlambda * ntau)
  colnames(penalty_parameter_sequence) = c("lambda", "tau")

  for (t in 1:ntau)
  {
    dm    <- make_design_matrix_sgl_fused(modelXlist, tau0[t])
    tau_s <- dm$tau_s
    x <- sqrt(total_n) * dm$design

    data = list(x = x, y = y)

    if (!is.null(lambda))
    {
      SGL_model = SGL(data = data, index = rep(1 : p, q + 1), min.frac = 0.05,
                      lambdas = lambda, alpha = alpha, standardize = FALSE)
    } else
    {
      SGL_model = SGL(data = data, index = rep(1 : p, q + 1), min.frac = 0.05,
                      nlam = nlambda, alpha = alpha, standardize = FALSE)
      lambda <- SGL_model$lambda
    }

    beta_all <- matrix(0, nrow = q * p, ncol = ncol(SGL_model$beta))

    for (l in 1:ncol(SGL_model$beta))
    {
      beta_all[,l] <- recover_beta_sgl_fused(SGL_model$beta[,l], p, q, tau_s)
      penalty_parameter_sequence[(l - 1) * ntau + t,] = c(lambda[l], tau0[t])
    }

    for (ind in 1:nlambda)
    {
      beta = matrix(beta_all[, ind], nrow = q, ncol = p, byrow = TRUE)
      intercept = unlist(Ybarlist) - apply(beta * matrix(unlist(Xbarlist), nrow = q, ncol = p, byrow = TRUE), 1, sum)
      beta = beta / matrix(unlist(Xsdlist), nrow = q, ncol = p, byrow = TRUE)

      betalist[[(ind - 1) * ntau + t]]      <- beta
      interceptlist[[(ind - 1) * ntau + t]] <- intercept
      #betalist[[(ind1 - 1) * nlambda2 + ind2]]      = beta
      #interceptlist[[(ind1 - 1) * nlambda2 + ind2]] = intercept

    }

  }


  return(list(interceptlist              = interceptlist,
              betalist                   = betalist,
              lambda                     = lambda,
              penalty_parameter_sequence = penalty_parameter_sequence))
}
