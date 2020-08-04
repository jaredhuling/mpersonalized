sparse_group_lasso_method = function(modelYlist, modelXlist, Wlist,
                                     Ybarlist,
                                     Xbarlist, Xsdlist, lambda, alpha,
                                     surrogate = c("squared_error", "logistic"), standardize = TRUE)
{

  surrogate <- match.arg(surrogate)

  q = length(modelXlist)
  p = dim(modelXlist[[1]])[2]

  x = as.matrix(bdiag(modelXlist))
  y = unlist(modelYlist)
  W = unname(unlist(Wlist))
  total_n = length(y)

  if (surrogate == "squared_error")
  {
    if (standardize) x = sqrt(total_n) * x; y = sqrt(total_n) * y
    type <- "linear"
  } else
  {
    type <- "logit"
  }

  data = list(x = x, y = y)

  SGL_model = SGL2(data = data, index = rep(1 : p, q), type = type,
                   weights = W,
                   lambdas = lambda, alpha = alpha, standardize = FALSE)

  nlambda = length(lambda)
  betalist = vector("list", nlambda)
  interceptlist = vector("list", nlambda)

  for (ind in 1:nlambda)
  {

      beta = matrix(SGL_model$beta[, ind], nrow = q, ncol = p, byrow = TRUE)

      lin_pred <- apply(beta * matrix(unlist(Xbarlist), nrow = q, ncol = p, byrow = TRUE), 1, sum)
      mu_hat <- unlist(Ybarlist)
      if (surrogate == "logistic")
      {
        intercept = rep(0, length(mu_hat))
      } else
      {
        intercept = mu_hat - lin_pred
        beta = beta / matrix(unlist(Xsdlist), nrow = q, ncol = p, byrow = TRUE)
      }

      betalist[[ind]]      = beta
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


make_penalty_factor <- function(modelXlist, tau0)
{
  nv.vec      <- unlist(lapply(modelXlist, ncol))
  nv.tot      <- sum(nv.vec)
  ss.vec      <- unlist(lapply(modelXlist, nrow))
  ss.tot      <- sum(ss.vec)
  tau_s       <- sqrt(ss.tot / ss.vec) * tau0
  studies     <- length(modelXlist)
  penfact.vec <- c(rep(1, ncol(modelXlist[[1]])),
                   unlist(lapply( 1:studies, function(s) rep(tau_s[s], nv.vec[s]) ))  )

  list(penfact = penfact.vec, tau_s = tau_s)
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

sparse_group_fused_lasso_method = function(modelYlist, modelXlist, Wlist,
                                           Ybarlist,
                                           Xbarlist, Xsdlist, lambda = NULL, alpha,
                                           fused = FALSE, tau0 = 1, nlambda = 50,
                                           surrogate = c("squared_error", "logistic"),
                                           standardize = TRUE)
{

  surrogate <- match.arg(surrogate)
  q = length(modelXlist)

  standardize <- standardize & !(surrogate == "logistic")

  if (!standardize & surrogate == "squared_error")
  {
    modelXlist <- lapply(modelXlist, function(x) cbind(1, x))
    nc.vec <- sapply(modelXlist, ncol)

    pf.mult <- c(0, rep(1, nc.vec[1]-1), unlist(lapply(nc.vec, function(nv) c(0, rep(1, nv-1)) )) )
  } else
  {

    nc.vec <- sapply(modelXlist, ncol)

    pf.mult <- c(rep(1, nc.vec[1]), unlist(lapply(nc.vec, function(nv) rep(1, nv))) )
  }

  #nc.vec <- sapply(modelXlist, ncol)
  #pf.mult <- c(rep(1, nc.vec[1]), unlist(lapply(nc.vec, function(nv) rep(1, nv))) )

  p = dim(modelXlist[[1]])[2]

  y = unlist(modelYlist)
  W = unname(unlist(Wlist))
  total_n = length(y)

  if (!is.null(lambda))
  {
    nlambda = length(lambda)
  }

  ntau    = length(tau0)

  betalist      = vector("list", nlambda * ntau)
  interceptlist = vector("list", nlambda * ntau)

  penalty_parameter_sequence = matrix(0, ncol = 2, nrow = nlambda * ntau)
  colnames(penalty_parameter_sequence) = c("lambda", "tau")


  if (surrogate == "squared_error")
  {
    if (standardize) y = sqrt(total_n) * y
    type <- "linear"
    lambda_list <- lambda <- NULL
  } else
  {
    type <- "logit"
    lambda_list <- vector(mode = "list", length = ntau)
  }


  for (t in 1:ntau)
  {
    #dm    <- make_design_matrix_sgl_fused(modelXlist, tau0[t])
    dm    <- make_design_matrix_sgl_fused(modelXlist, 1)
    pflist <- make_penalty_factor(modelXlist, tau0[t])



    tau_s   <- dm$tau_s
    penfact <- drop(pflist$penfact) * pf.mult
    if (surrogate == "squared_error")
    {
      if (standardize)
      {
        x <- sqrt(total_n) * dm$design
      } else
      {
        x <- dm$design
      }
    } else
    {
      x <- dm$design
    }

    data = list(x = x, y = y)

    #penfact <- rep(1, length(penfact))

    if (!is.null(lambda))
    {
      SGL_model = SGL2(data = data, index = rep(1 : p, q + 1), min.frac = 0.05,
                       type = type, weights = W, penaltyweights = penfact,
                       lambdas = lambda, alpha = alpha, standardize = FALSE)
    } else
    {
      SGL_model = SGL2(data = data, index = rep(1 : p, q + 1), min.frac = 0.05,
                       type = type, weights = W, penaltyweights = penfact,
                       nlam = nlambda, alpha = alpha, standardize = FALSE)
      lambda <- lambda_list[[t]] <- SGL_model$lambda
    }

    beta_all <- matrix(0, nrow = q * p, ncol = ncol(SGL_model$beta))

    for (l in 1:ncol(SGL_model$beta))
    {
      #beta_all[,l] <- recover_beta_sgl_fused(SGL_model$beta[,l], p, q, tau_s)
      beta_all[,l] <- recover_beta_sgl_fused(SGL_model$beta[,l], p, q, rep(1, length(tau_s)))
      penalty_parameter_sequence[(l - 1) * ntau + t,] = c(lambda[l], tau0[t])
    }

    for (ind in 1:nlambda)
    {
      beta = matrix(beta_all[, ind], nrow = q, ncol = p, byrow = TRUE)



      mu_hat <- unlist(Ybarlist)
      if (surrogate == "logistic")
      {
        intercept = rep(0, length(mu_hat))
        if (!is.null(SGL_model$intercept))
        {
          intercept <- rep(SGL_model$intercept[ind], length(mu_hat))
        }
      } else
      {
        lin_pred <- apply(beta * matrix(unlist(Xbarlist), nrow = q, ncol = p, byrow = TRUE), 1, sum)
        intercept = mu_hat - lin_pred
        #intercept = rep(0, length(mu_hat))
        beta = beta / matrix(unlist(Xsdlist), nrow = q, ncol = p, byrow = TRUE)
      }

      if (!standardize & surrogate == "squared_error")
      {
        intercept <- beta[,1]
        beta      <- beta[,-1]
      }


      betalist[[(ind - 1) * ntau + t]]      <- beta
      interceptlist[[(ind - 1) * ntau + t]] <- intercept
      #betalist[[(ind1 - 1) * nlambda2 + ind2]]      = beta
      #interceptlist[[(ind1 - 1) * nlambda2 + ind2]] = intercept

    }

  }


  return(list(interceptlist              = interceptlist,
              betalist                   = betalist,
              lambda                     = lambda,
              lambda_list                = lambda_list,
              penalty_parameter_sequence = penalty_parameter_sequence))
}




