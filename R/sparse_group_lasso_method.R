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








fit_sparse_group_fused_lasso = function(modelYlist, modelXlist,
                                        lambda = NULL,
                                        lambda_list = NULL,
                                        alpha = 0.85,
                                        tau0 = 1, nlambda = 50,
                                        type = c("linear", "logit"),
                                        standardize = TRUE)
{

  type <- match.arg(type)

  q = length(modelXlist)

  if (!is.null(lambda_list))
  {
    lambda_list_given <- TRUE
  } else
  {
    lambda_list_given <- FALSE
  }


  standardize <- standardize & !(type == "logit")

  if (!standardize & type == "linear")
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
  total_n = length(y)

  if (!is.null(lambda))
  {
    nlambda = length(lambda)
  }

  if (!is.null(lambda_list))
  {
    nlambda = length(lambda_list[[1]])
  }

  ntau    = length(tau0)

  betalist      = vector("list", nlambda * ntau)
  interceptlist = vector("list", nlambda * ntau)

  penalty_parameter_sequence = matrix(0, ncol = 2, nrow = nlambda * ntau)
  colnames(penalty_parameter_sequence) = c("lambda", "tau")


  if (type == "linear")
  {
    type <- "linear"
    lambda <- NULL
  } else
  {
    type <- "logit"
    if (is.null(lambda_list))
    {
      lambda_list <- vector(mode = "list", length = ntau)
    }
  }



  for (t in 1:ntau)
  {
    #dm    <- make_design_matrix_sgl_fused(modelXlist, tau0[t])
    dm    <- make_design_matrix_sgl_fused(modelXlist, 1)
    pflist <- make_penalty_factor(modelXlist, tau0[t])



    tau_s   <- dm$tau_s
    penfact <- drop(pflist$penfact) * pf.mult
    if (type == "linear")
    {
      x <- dm$design
    } else
    {
      x <- dm$design
    }

    data = list(x = x, y = y)

    #penfact <- rep(1, length(penfact))

    if (!lambda_list_given)
    {
      if (!is.null(lambda))
      {
        SGL_model = SGL2(data = data, index = rep(1 : p, q + 1), min.frac = 0.05,
                         type = type, penaltyweights = penfact,
                         lambdas = lambda, alpha = alpha, standardize = standardize)

        lambda_list <- rep(list(lambda), ntau)
      } else
      {
        SGL_model = SGL2(data = data, index = rep(1 : p, q + 1), min.frac = 0.05,
                         type = type, penaltyweights = penfact,
                         nlam = nlambda, alpha = alpha, standardize = standardize)
        lambda <- lambda_list[[t]] <- SGL_model$lambda
      }
    } else
    {
      SGL_model = SGL2(data = data, index = rep(1 : p, q + 1), min.frac = 0.05,
                      type = type, penaltyweights = penfact,
                      lambdas = lambda_list[[t]], alpha = alpha, standardize = standardize)
      lambda <- lambda_list[[t]]
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



      if (type == "logit")
      {
        intercept = rep(0, length(modelXlist))
        if (!is.null(SGL_model$intercept))
        {
          intercept <- rep(SGL_model$intercept[ind], length(modelXlist))
        }
      } else if (type == "linear")
      {
        intercept = rep(0, length(modelXlist))
        if (!is.null(SGL_model$intercept))
        {
          intercept <- rep(SGL_model$intercept[1], length(modelXlist))
        }
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


fit_sparse_group_fused_lasso_cv = function(modelYlist, modelXlist,
                                           lambda = NULL,
                                           nlambda = 50,
                                           type = c("linear", "logit"),
                                           standardize = TRUE,
                                           tau0 = NULL,
                                           num_tau0    = ifelse(!is.null(tau0), length(tau0), 11),
                                           min_tau     = 1e-2,
                                           alpha = 0.85, cv_folds = 5)
{

  type <- match.arg(type)

  q <- length(modelXlist)
  p <- dim(modelXlist[[1]])[2]

  if (!is.null(alpha))
  {
    if (alpha < 0 | alpha > 1)
    {
      warning("salpha must be between 0 and 1. The default is 0.95!")
      alpha = 0.95
    }
  } else alpha = 0.95

  if (is.null(tau0))
  {
    tau0 <- gen_tau0(num_tau0, min_tau)
  }

  fit_all_data <- fit_sparse_group_fused_lasso(modelYlist = modelYlist,
                                               modelXlist = modelXlist,
                                               lambda = lambda, alpha = alpha,
                                               tau0 = tau0, nlambda = nlambda,
                                               type = type,
                                               standardize = standardize)

  lambda_list <- fit_all_data$lambda_list


  tune_cost = matrix(0, nrow = if(is.null(lambda_list)){nlambda}else{length(lambda_list[[1]])}, ncol = length(tau0))
  rownames(tune_cost) <- paste0("lam=", round(lambda_list[[1]], 2))
  colnames(tune_cost) <- paste0("tau0=", round(tau0, 2))

  #build cross-validation folds
  folds_index <- lapply(modelXlist, function(x) createFolds(1:dim(x)[1], k = cv_folds))


  #carry out method for each cross validation fold
  for (k in 1:cv_folds)
  {

    cv_Xlist = vector("list", q); left_Xlist = vector("list", q)
    cv_Ylist = vector("list", q); left_Ylist = vector("list", q)

    for (j in 1:q) {
      cv_Xlist[[j]] = modelXlist[[j]][-folds_index[[j]][[k]],,drop = FALSE]
      cv_Ylist[[j]] = modelYlist[[j]][-folds_index[[j]][[k]]]

      left_Xlist[[j]] = modelXlist[[j]][folds_index[[j]][[k]],,drop = FALSE]
      left_Ylist[[j]] = modelYlist[[j]][folds_index[[j]][[k]]]
    }


    cv_model <- fit_sparse_group_fused_lasso(modelYlist = cv_Ylist,
                                             modelXlist = cv_Xlist,
                                             lambda_list = lambda_list,
                                             alpha = alpha,
                                             tau0 = tau0,
                                             nlambda = nlambda,
                                             type = type,
                                             standardize = standardize)

    cv_interceptlist = cv_model$interceptlist
    cv_betalist = cv_model$betalist

    nlam <- length(lambda_list[[1]])


    for (ind1 in 1:nlam)
    {
      for (ind2 in 1:length(tau0))
      {
        if (type == "linear")
        {
          tune_cost[ind1, ind2] = tune_cost[ind1, ind2] + sum(unlist(mapply(function(y, x, intercept, beta) {sum((y - drop(intercept) - drop(x %*% beta) ) ^ 2)},
                                                                            y = left_Ylist, x = left_Xlist,
                                                                            intercept = as.list(cv_interceptlist[[(ind1 - 1) * length(tau0) + ind2]]),
                                                                            beta = split(cv_betalist[[(ind1 - 1) * length(tau0) + ind2]], row(cv_betalist[[(ind1 - 1) * length(tau0) + ind2]])), SIMPLIFY = FALSE)))

        } else if (type == "logit")
        {
          tune_cost[ind1, ind2] = tune_cost[ind1, ind2] + mean(unlist(mapply(function(y, x, intercept, beta) glmnet:::auc(y, 1 / (1 + exp(- drop(intercept) - drop(x %*% beta)))),
                                                                            y = left_Ylist, x = left_Xlist,
                                                                            intercept = as.list(cv_interceptlist[[(ind1 - 1) * length(tau0) + ind2]]),
                                                                            beta = split(cv_betalist[[(ind1 - 1) * length(tau0) + ind2]], row(cv_betalist[[(ind1 - 1) * length(tau0) + ind2]])), SIMPLIFY = FALSE)),
                                                               na.rm = TRUE)

        }
      }
    }


  } ## end CV loop

  tune_cost <- tune_cost / cv_folds

  if (type == "linear")
  {
    opt_ind = which(tune_cost == min(tune_cost), arr.ind = TRUE)
    opt_ind1 = opt_ind[1]; opt_ind2 = opt_ind[2]
  } else if (type == "logit")
  {
    opt_ind = which(tune_cost == max(tune_cost), arr.ind = TRUE)
    if (!is.null(dim(opt_ind)))
    {
      opt_ind1 = opt_ind[1,1]; opt_ind2 = opt_ind[1,2]
    } else
    {
      opt_ind1 = opt_ind[1]; opt_ind2 = opt_ind[2]
    }
  }


  penalty_parameter_sequence = fit_all_data$penalty_parameter_sequence


  intercept <- fit_all_data$interceptlist[[(opt_ind1 - 1) * length(tau0) + opt_ind2]]
  beta      <- fit_all_data$betalist[[(opt_ind1 - 1) * length(tau0) + opt_ind2]]

  opt_penalty_parameter = penalty_parameter_sequence[opt_ind1,]

  opt_lambda <- fit_all_data$lambda_list[[opt_ind2]][opt_ind1]
  opt_tau0   <- tau0[opt_ind2]

  list(intercept = intercept,
       beta = beta,
       cvm = tune_cost,
       opt_lambda = opt_lambda,
       opt_tau0 = opt_tau0)
}
