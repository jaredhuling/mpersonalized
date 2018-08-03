#' @title Cross Validation for \code{mpersonalized}
#'
#' @description  This function implments \code{mpersonalized} and use cross validatation to tune penalty parameter.
#'  The optimal penalty parameter is selected by minimizing \deqn{\sum_{i=1}^{n_k}\frac{|\hat{C}_k(X_{i})|}{\sum_{i=1}^{n_k}|\hat{C}_k(X_{i})|}\bigl [1\{\hat{C}_k(X_{i})>0\}-g_k(X_{i})\bigr]^2}
#'  in the leave-out fold, where \eqn{\hat{C}_k(X_{i})} in the leave-out fold is independently estimated from the training set.
#'
#' @param problem A character string specifiy whether user want to solve "meta-analysis" or
#' "multiple outcomes" problem. For \code{problem = "meta-analysis"}, user should also supply
#' \code{Xlist}, \code{Ylist}, \code{Trtlist}. For \code{problem = "multiple outcomes"},
#' user should supply \code{X}, \code{Ylist}, \code{Trt}.
#' @param X Covariate matrix that should be supplied when \code{problem = "multiple outcomes"}
#' with rows indicating subjects and columns indicating covariates.
#' @param Trt Treatment vector that should be supplied when \code{problem = "multiple outcomes"},
#' which should be coded as 0 or 1.
#' @param P Propensity score vector when \code{problem = "multiple outcomes"}. If not supplied,
#' then study is treated as randomzied trial and the propensity score is estimated as the proportion
#' of 1's in \code{Trt} for every subject.
#' @param Xlist A list object that should be supplied when \code{problem = "meta-analysis"},
#' with \eqn{k}th element denoting the covariate matrix of study \eqn{k}.
#' @param Ylist When \code{problem = "meta-analysis"}, \code{Ylist} should be a list object with \eqn{k}th element
#' denoting the response vector of study \eqn{k}. When \code{problem = "multiple outcomes"}, \code{Ylist} should
#' be a list object with \eqn{k}th element denoting the \eqn{k}th outcome.
#' @param Trtlist A list object that should be supplied when \code{problem = "meta-analysis"},
#' with \eqn{k}th element denoting the treatment vector of study \eqn{k} (coded as 0 or 1).
#' @param Plist A list object that should be supplied when \code{problem = "meta-analysis"},
#' with \eqn{k}the element denoting the propensity score vector of study \eqn{k}.
#' If not supplied, then each study is treated as randomized trial and the corresponding propensity score
#' is estimated as the proportion of 1's in the \eqn{k}th element of \code{Trtlist} for all subjects.
#' @param typelist A list object with \eqn{k}th element denoting the type of outcome corresponding
#' to the \eqn{k}th element in \code{Ylist}. Each element could be "continuous" or "binary".
#' @param penalty For different rules, the penalty could be "lasso", "GL", "SGL", "fused",
#' "lasso+fused", "GL+fused", "SGL+fused", or "SGL+SL". For single rule, the penalty could only be "lasso".
#' For \code{penalty = "none"}, use function \code{mpersonalized} instead.
#' User should always input \code{penalty} and then supply correponding penalty parameters sequence
#' if needed.
#' @param lambda1 \eqn{\lambda_1} in the framework of different rules. If not supplied, a default
#' sequence will be computed.
#' @param lambda2 \eqn{\lambda_2} in the framework of different rules. If not supplied, a default
#' sequence will be computed.
#' @param tau0 Parameter \eqn{\tau_0} for the \code{"SGL+SL"} penalty in the framework of different rules.
#' If not supplied, a default sequence will be computed.
#' @param alpha \eqn{\alpha} in the framework of different rules. If not supplied, a default value
#' will be used depending on \code{penalty}.
#' @param single_rule_lambda \eqn{\lambda_{single}} in the framework of single rule.
#' @param single_rule A logical value, whether the single treatment framework is used. Deafult is \code{FALSE}.
#' @param cv_folds Number of folds needed for cross-validation. Default is 5
#' @param admm_control A list of parameters which user can specify to control the admm algorithm.
#' In \code{admm_control}, the following parameters can be supplied:
#' \code{abs.tol}, absolute tolerance; \code{rel.tol}, relative tolerance; \code{maxit}, maximum number of iterations;
#' \code{rho}, Lagrangian parameter.
#' @param contrast_builder_control A list of parameters which user can specify to control estimation of
#' contrast function. In \code{contrast_builder_control},
#' the following parameters could be supplied: \code{eff_aug}, a logical value whether efficiency augmentation
#' should be implemented; \code{response_model}, a character string specify what outcome model to use
#' if \code{eff_aug = TRUE}, \code{response_model} could be "lasso" or "linear";
#' \code{contrast_builder_folds}, the number of folds used in cross validation when \code{response_model = "lasso"}.
#' @param num_lambda1 If \code{lambda1} is not specified by user, user could still specify the length of the
#' \code{lambda1} sequence. The default length is 10.
#' @param num_lambda2 If \code{lambda2} is not specified by user, user could still specify the length of the
#' \code{lambda2} sequence. The default length is 10.
#' @param num_tau0 If \code{tau0} is not specified by user, the user can still specify the length of the
#' \code{tau0} sequence. The default length is 11.
#' @param min_tau If \code{tau0} is not specified by user, \code{min_tau} specifies the minimum value
#' for \eqn{\tau_0}. The largest value for \eqn{\tau_0} will be \code{1 / min_tau}.
#' @param num_single_rule_lambda If \code{single_rule_lambda} is not specified, user could still specify the length
#' of the \code{single_rule_lambda} sequence. The default length is 50.
#'
#' @import glmnet SGL caret Matrix
#'
#' @return An S3 object of class "mp_cv", which contains the information of the model with the optimal lambda. It can be supplied
#' to some other functions in mperosnalized package for further analysis or prediction.
#'
#' \item{penalty_parameter_sequence}{A matrix object with each row denoting a configuration of the penalty parameters.}
#' \item{opt_penalty_parameter}{Optimal penalty parameter chosen by minimizing the cross validation error.}
#' \item{intercept}{The vector of intercepts corresponding to the optimal penalty parameter.}
#' \item{beta}{The coefficient matrix corresponding to the optimal penalty parameter.}
#' \item{number_covariates}{Number of candidate covariates considered.}
#' \item{number_studies_or_outcomes}{Number of studies if \code{problem = "meta-analysis"} or number of outcomes if \code{problem = "multiple outcomes"}.}
#'
#' @examples
#' set.seed(123)
#' sim_dat = simulated_dataset(n = 200, problem = "meta-analysis")
#' Xlist = sim_dat$Xlist; Ylist = sim_dat$Ylist; Trtlist = sim_dat$Trtlist
#'
#' # fit different rules with group lasso penalty
#' mp_cvmod_diff = mpersonalized_cv(problem = "meta-analysis",
#'                                  Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist,
#'                                  penalty = "GL", single_rule = FALSE)
#'
#' mp_cvmod_diff$intercept
#' mp_cvmod_diff$beta
#'
#' # fit a single rule with lasso penalty
#' mp_cvmod_single = mpersonalized_cv(problem = "meta-analysis",
#'                                    Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist,
#'                                    penalty = "lasso", single_rule = TRUE)
#'
#' mp_cvmod_single$intercept
#' mp_cvmod_single$beta
#' set.seed(NULL)
#' @export

mpersonalized_cv = function(problem = c("meta-analysis", "multiple outcomes"),
                            X, Trt, P = NULL,
                            Xlist, Ylist, Trtlist, Plist = replicate(length(Xlist), NULL, simplify = FALSE),
                            typelist = replicate(length(Xlist), "continuous", simplify = FALSE),
                            penalty = c("lasso", "GL", "SGL", "fused",
                                      "lasso+fused", "GL+fused", "SGL+fused",
                                      "SGL+SL"),
                            lambda1 = NULL, lambda2 = NULL, tau0 = NULL,
                            single_rule_lambda = NULL,
                            num_lambda1 = ifelse(!is.null(lambda1), length(lambda1),10),
                            num_lambda2 = ifelse(!is.null(lambda2), length(lambda2),10),
                            num_tau0    = ifelse(!is.null(tau0), length(tau0), 11),
                            min_tau     = 1e-2,
                            num_single_rule_lambda = ifelse(!is.null(single_rule_lambda), length(single_rule_lambda), 50),
                            alpha = NULL, single_rule = FALSE, cv_folds = 5,
                            admm_control = NULL,
                            contrast_builder_control = NULL){

  penalty = match.arg(penalty)
  problem = match.arg(problem)

  if (problem == "multiple outcomes"){

    if (missing(X) | missing(Ylist) | missing(Trt))
      stop("For multiple outcomes, X, Ylist, Trt need to be supplied!")

    q = length(Ylist)
    Xlist = replicate(q, X, simplify = FALSE)
    Trtlist = replicate(q, Trt, simplify = FALSE)
    #a default estimate for the propensity score
    if (is.null(P)){
      P = rep(sum(Trt) / length(Trt), length(Trt))
    } else {
      if (length(P) == 1)
        P = rep(P, length(Trt))
    }

    Plist = replicate(q, P, simplify = FALSE)

  } else if (problem == "meta-analysis"){

    if (missing(Xlist) | missing(Ylist) | missing(Trtlist))
      stop("For meta-analysis, Xlist, Ylist, Trtlist need to be supplied!")
    #a default estimate for the propensity score
    Plist = mapply(
      function(P, Trt){
        if (is.null(P)){
          P = rep(sum(Trt) / length(Trt), length(Trt))
        } else {
          if (length(P) == 1)
            P = rep(P, length(Trt))
        }
        return(P)
      },
      P = Plist,
      Trt = Trtlist,
      SIMPLIFY = FALSE)
  }

  q = length(Xlist)
  p = dim(Xlist[[1]])[2]


  #construct contrast for the data
  Conlist = vector("list", q)
  for (j in 1:q){
    Conlist[[j]] = do.call(contrast_builder, c(list(X = Xlist[[j]],
                                                    Y = Ylist[[j]],
                                                    ori_Trt = Trtlist[[j]],
                                                    P = Plist[[j]],
                                                    type = typelist[[j]]),
                                               contrast_builder_control))
  }


  standardized_data = contrast_standardize(Conlist = Conlist, Xlist = Xlist,
                                           single_rule = single_rule)
  modelYlist = standardized_data$modelYlist
  modelXlist = standardized_data$modelXlist

  #check whether the information provided is correct and set up the value for lambdas if not provided
  if (single_rule == TRUE){

    Ybar = standardized_data$Ybar
    Xbar = standardized_data$Xbar
    Xsd = standardized_data$Xsd

    if (penalty != "lasso")
      stop("When single rule is required, the penalty must be lasso!(for penalty = none, use function 'mpersonalized' instead.")

    if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
      warning("When single rule = TRUE, the value for lambda1, lambda2, alpha are ignored!")

    if (is.null(single_rule_lambda)){
      lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                       penalty = penalty, single_rule = single_rule,
                                       num_single_rule_lambda = num_single_rule_lambda)

      single_rule_lambda = lambda_default$single_rule_lambda
    }

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (!is.null(single_rule_lambda))
      warning("When single rule = FALSE, the value for single_rule_lambda is ignored!")

    if (penalty == "none")
      stop("For penalty = none, use function 'mpersonalized' instead.")

    if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")) {

      if (penalty == "fused"){
        if (!is.null(alpha)){
          warning("When penalty = fused, values of alpha is ignored!")
          alpha = NULL
        }

        if (!is.null(lambda1))
          if (sum(lambda1 != 0) > 0)
            warning("When penalty = fused, value of lambda1 is automatically set to be 0!")

        lambda1 = 0
      }

      if (penalty == "lasso+fused"){

        if (!is.null(alpha))
          if(alpha != 1)
            warning("When penalty = lasso+fused, alpha is automatically set to be 1!")

        alpha = 1

      } else if (penalty == "GL+fused"){

        if (!is.null(alpha))
          if(alpha != 0)
            warning("When penalty = GL+fused, alpha is automatically set to be 0!")

        alpha = 0

      } else if (penalty == "SGL+fused"){

        if (!is.null(alpha)){
          if (alpha == 0 | alpha == 1){
            warning("When penalty = SGL+fused, alpha cannot be set as 0 or 1, and default is 0.95!")
            alpha = 0.95
          }
        } else alpha = 0.95
      }

      if (is.null(lambda1) | is.null(lambda2)){
        lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                         penalty = penalty, single_rule = single_rule, alpha = alpha,
                                         num_lambda1 = num_lambda1, num_lambda2 = num_lambda2,
                                         lambda1 = lambda1, lambda2 = lambda2)

        if (is.null(lambda1))
          lambda1 = lambda_default$lambda1

        if (is.null(lambda2))
          lambda2 = lambda_default$lambda2
      }

    } else if (penalty %in% c("lasso", "GL", "SGL", "SGL+SL")){

      if (!is.null(lambda2)){
        if (sum(lambda2 != 0) > 0){
          warning("When penalty = lasso/GL/SGL/SGL+SL, the value for lambda2 is ignored and automatically set to be 0!")
        }
      }

      lambda2 = 0

      if (penalty == "lasso"){

        if (!is.null(alpha))
          if (alpha != 1)
            warning("When penalty = lasso, alpha is automatically set to be 1!")

        alpha = 1

      } else if (penalty == "GL"){

        if (!is.null(alpha))
          if (alpha != 0)
            warning("When penalty = GL, alpha is automatically set to be 0!")

        alpha = 0

      } else if (penalty == "SGL"){

        if (!is.null(alpha)){
          if (alpha == 0 | alpha == 1){
            warning("When penalty = SGL, alpha cannot be set as 0 or 1, and default is 0.95!")
            alpha = 0.95
          }
        } else alpha = 0.95
      } else if (penalty == "SGL+SL"){

        if (!is.null(alpha)){
          if (alpha < 0 | alpha > 1){
            warning("When penalty = SGL+SL, alpha must be between 0 and 1. The default is 0.95!")
            alpha = 0.95
          }
        } else alpha = 0.95
      }

      if (is.null(lambda1)){ # & penalty != "SGL+SL"){
        lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                         penalty = penalty, single_rule = single_rule, alpha = alpha,
                                         num_lambda1 = num_lambda1, lambda1 = lambda1)

        lambda1 = lambda_default$lambda1
      }

      if (penalty == "SGL+SL")
      {
        if (is.null(tau0))
        {
          tau0 <- gen_tau0(num_tau0, min_tau)
        }
      }
    }
  }

  #build cross-validation folds
  folds_index = lapply(Xlist, function(x) createFolds(1:dim(x)[1], k = cv_folds))

  #determine the dimension for the cost for tuning penalty parameters
  if (single_rule == TRUE){
    tune_cost = numeric(length(single_rule_lambda))
  } else {
    if (penalty %in% c("lasso", "GL", "SGL")){
      tune_cost = numeric(length(lambda1))
      names(tune_cost) <- paste0("lam=", round(lambda1, 2))
    } else if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){
      tune_cost = matrix(0, nrow = length(lambda1), ncol = length(lambda2))
      rownames(tune_cost) <- paste0("lam1=", round(lambda1, 2))
      colnames(tune_cost) <- paste0("lam2=", round(lambda2, 2))
    } else if (penalty %in% c("SGL+SL")){
      tune_cost = matrix(0, nrow = if(is.null(lambda1)){num_lambda1}else{length(lambda1)}, ncol = length(tau0))
      rownames(tune_cost) <- paste0("lam=", round(lambda1, 2))
      colnames(tune_cost) <- paste0("tau0=", round(tau0, 2))
    }
  }

  #carry out method for each cross validation fold
  for (k in 1:cv_folds){

    cv_Xlist = vector("list", q); left_Xlist = vector("list", q)
    cv_Ylist = vector("list", q); left_Ylist = vector("list", q)
    cv_Trtlist = vector("list", q); left_Trtlist = vector("list", q)
    cv_Plist = vector("list", q); left_Plist = vector("list", q)

    for (j in 1:q){
      cv_Xlist[[j]] = Xlist[[j]][-folds_index[[j]][[k]],]
      cv_Ylist[[j]] = Ylist[[j]][-folds_index[[j]][[k]]]
      cv_Trtlist[[j]] = Trtlist[[j]][-folds_index[[j]][[k]]]
      cv_Plist[[j]] = Plist[[j]][-folds_index[[j]][[k]]]

      left_Xlist[[j]] = Xlist[[j]][folds_index[[j]][[k]],]
      left_Ylist[[j]] = Ylist[[j]][folds_index[[j]][[k]]]
      left_Trtlist[[j]] = Trtlist[[j]][folds_index[[j]][[k]]]
      left_Plist[[j]] = Plist[[j]][folds_index[[j]][[k]]]
    }

    cv_Conlist = vector("list", q)
    left_Conlist = vector("list", q)

    for (j in 1:q){

      cv_Conlist[[j]] = do.call(contrast_builder, c(list(X = cv_Xlist[[j]],
                                                         Y = cv_Ylist[[j]],
                                                         ori_Trt = cv_Trtlist[[j]],
                                                         P = cv_Plist[[j]],
                                                         type = typelist[[j]]),
                                                    contrast_builder_control))

      left_Conlist[[j]] = do.call(contrast_builder, c(list(X = left_Xlist[[j]],
                                                           Y = left_Ylist[[j]],
                                                           ori_Trt = left_Trtlist[[j]],
                                                           P = left_Plist[[j]],
                                                           type = typelist[[j]]),
                                                      contrast_builder_control))
    }

    cv_standardized_data = contrast_standardize(Conlist = cv_Conlist, Xlist = cv_Xlist, single_rule = single_rule)
    cv_modelYlist = cv_standardized_data$modelYlist
    cv_modelXlist = cv_standardized_data$modelXlist

    #transform contrast into binary data for the left out fold
    left_sConlist = lapply(left_Conlist, function(y) as.numeric(y > 0))
    left_Wlist = lapply(left_Conlist, abs)

    left_dataWlist = lapply(left_Wlist, sum)
    left_adj_Wlist = mapply(function(w, dataw) w / dataw, w = left_Wlist,
                       dataw = left_dataWlist, SIMPLIFY = FALSE)

    if (single_rule == TRUE){
      cv_Ybar = cv_standardized_data$Ybar
      cv_Xbar = cv_standardized_data$Xbar
      cv_Xsd = cv_standardized_data$Xsd

      cv_model = single_rule_lasso_method(modelYlist = cv_modelYlist, modelXlist = cv_modelXlist,
                                          Ybar = cv_Ybar, Xbar = cv_Xbar, Xsd = cv_Xsd, lambda = single_rule_lambda)
      cv_interceptlist = cv_model$interceptlist
      cv_betalist = cv_model$betalist

      for (ind in 1:length(single_rule_lambda))
        tune_cost[ind] = tune_cost[ind] + sum(unlist(mapply(function(w, y, x) sum(w * (y - cv_interceptlist[[ind]] - x %*% cv_betalist[[ind]]) ^ 2),
                                                            w = left_adj_Wlist, y = left_sConlist, x = left_Xlist, SIMPLIFY = FALSE)))

    } else {
      cv_Ybarlist = cv_standardized_data$Ybarlist
      cv_Xbarlist = cv_standardized_data$Xbarlist
      cv_Xsdlist = cv_standardized_data$Xsdlist

      if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){

        cv_model = meta_method(modelYlist = cv_modelYlist, modelXlist = cv_modelXlist,
                               Ybarlist = cv_Ybarlist, Xbarlist = cv_Xbarlist, Xsdlist = cv_Xsdlist,
                               lambda1 = lambda1, lambda2 = lambda2, alpha = alpha, admm_control = admm_control)
        cv_interceptlist = cv_model$interceptlist
        cv_betalist = cv_model$betalist

        for (ind1 in 1:length(lambda1))
          for (ind2 in 1:length(lambda2))
            tune_cost[ind1, ind2] = tune_cost[ind1, ind2] + sum(unlist(mapply(function(w, y, x, intercept, beta) sum(w * (y - intercept - x %*% beta) ^ 2),
                                                                              w = left_adj_Wlist, y = left_sConlist, x = left_Xlist,
                                                                              intercept = as.list(cv_interceptlist[[(ind1 - 1) * length(lambda2) + ind2]]),
                                                                              beta = split(cv_betalist[[(ind1 - 1) * length(lambda2) + ind2]], row(cv_betalist[[(ind1 - 1) * length(lambda2) + ind2]])), SIMPLIFY = FALSE)))

      } else if (penalty %in% c("lasso", "SGL", "GL")) {

        cv_model = sparse_group_lasso_method(modelYlist = cv_modelYlist, modelXlist = cv_modelXlist,
                                             Ybarlist = cv_Ybarlist, Xbarlist = cv_Xbarlist, Xsdlist = cv_Xsdlist,
                                             lambda = lambda1, alpha = alpha)
        cv_interceptlist = cv_model$interceptlist
        cv_betalist = cv_model$betalist

        for (ind in 1:length(lambda1))
          tune_cost[ind] = tune_cost[ind] + sum(unlist(mapply(function(w, y, x, intercept, beta) sum(w * (y - intercept - x %*% beta) ^ 2),
                                                              w = left_adj_Wlist, y = left_sConlist, x = left_Xlist,
                                                              intercept = as.list(cv_interceptlist[[ind]]),
                                                              beta = split(cv_betalist[[ind]], row(cv_betalist[[ind]])), SIMPLIFY = FALSE)))

      }  else if (penalty %in% c("SGL+SL")) {
        cv_model = sparse_group_fused_lasso_method(modelYlist = cv_modelYlist, modelXlist = cv_modelXlist,
                                                   Ybarlist = cv_Ybarlist, Xbarlist = cv_Xbarlist, Xsdlist = cv_Xsdlist,
                                                   lambda = lambda1, alpha = alpha, tau0 = tau0, nlambda = num_lambda1)
        cv_interceptlist = cv_model$interceptlist
        cv_betalist = cv_model$betalist

        penalty_parameter_sequence_all <- cv_model$penalty_parameter_sequence

        if (is.null(lambda1))
        {
          nlam1 <- num_lambda1
        } else
        {
          nlam1 <- length(lambda1)
        }

        for (ind1 in 1:nlam1)
          for (ind2 in 1:length(tau0))
            tune_cost[ind1, ind2] = tune_cost[ind1, ind2] + sum(unlist(mapply(function(w, y, x, intercept, beta) sum(w * (y - intercept - x %*% beta) ^ 2),
                                                                              w = left_adj_Wlist, y = left_sConlist, x = left_Xlist,
                                                                              intercept = as.list(cv_interceptlist[[(ind1 - 1) * length(tau0) + ind2]]),
                                                                              beta = split(cv_betalist[[(ind1 - 1) * length(tau0) + ind2]], row(cv_betalist[[(ind1 - 1) * length(tau0) + ind2]])), SIMPLIFY = FALSE)))

      }
    }
  }

  #for the complete data set

  if (single_rule == TRUE){

    Ybar = standardized_data$Ybar
    Xbar = standardized_data$Xbar
    Xsd = standardized_data$Xsd

    full_model = single_rule_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                          Ybar = Ybar, Xbar = Xbar, Xsd = Xsd, lambda = single_rule_lambda)

    opt_ind = which.min(tune_cost)

    penalty_parameter_sequence = as.matrix(single_rule_lambda)
    colnames(penalty_parameter_sequence) = "single_rule_lambda"

    penalty_parameter_sequence_all = penalty_parameter_sequence

    model_info = list(intercept = full_model$interceptlist[[opt_ind]], beta = full_model$betalist[[opt_ind]],
                      penalty_parameter_sequence = penalty_parameter_sequence,
                      opt_penalty_parameter= penalty_parameter_sequence[opt_ind,],
                      cv_error = tune_cost,
                      penalty = "lasso", single_rule = TRUE,
                      number_covariates = p, number_studies_or_outcomes = q,
                      Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                      Ybar = Ybar, Xbar = Xbar, Xsd = Xsd,
                      problem = problem)

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){

      opt_ind = which(tune_cost == min(tune_cost), arr.ind = TRUE)
      opt_ind1 = opt_ind[1]; opt_ind2 = opt_ind[2]
      full_model = meta_method(modelYlist = modelYlist, modelXlist = modelXlist,
                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                               lambda1 = lambda1[opt_ind1], lambda2 = lambda2[opt_ind2], alpha = alpha, admm_control = admm_control)

      penalty_parameter_sequence = matrix(0, ncol = 2, nrow = length(lambda1) * length(lambda2))
      colnames(penalty_parameter_sequence) = c("lambda1", "lambda2")

      for (ind1 in 1:length(lambda1)){
        for (ind2 in 1:length(lambda2)){
          penalty_parameter_sequence[(ind1 - 1) * length(lambda2) + ind2,] = c(lambda1[ind1], lambda2[ind2])
        }
      }

      penalty_parameter_sequence_all = penalty_parameter_sequence

      opt_penalty_parameter = penalty_parameter_sequence[(opt_ind1 - 1) * length(lambda2) + opt_ind2,]

      model_info = list(intercept = full_model$interceptlist[[1]], beta = full_model$betalist[[1]],
                        penalty_parameter_sequence = penalty_parameter_sequence,
                        opt_penalty_parameter = opt_penalty_parameter,
                        cv_error = tune_cost,
                        alpha = alpha, penalty = penalty, single_rule = FALSE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        Ybar = Ybarlist, Xbar = Xbarlist, Xsd = Xsdlist,
                        problem = problem)

    } else if (penalty %in% c("lasso", "SGL", "GL")){

      full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                             Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                             lambda = lambda1, alpha = alpha)

      opt_ind = which.min(tune_cost)

      penalty_parameter_sequence = as.matrix(lambda1)
      colnames(penalty_parameter_sequence) = "lambda1"

      penalty_parameter_sequence_all = penalty_parameter_sequence

      model_info = list(intercept = full_model$interceptlist[[opt_ind]], beta = full_model$betalist[[opt_ind]],
                        penalty_parameter_sequence = penalty_parameter_sequence,
                        opt_penalty_parameter = penalty_parameter_sequence[opt_ind,],
                        cv_error = tune_cost,
                        alpha = alpha, penalty = penalty, single_rule = FALSE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        Ybar = Ybarlist, Xbar = Xbarlist, Xsd = Xsdlist,
                        problem = problem)

    } else if (penalty %in% c("SGL+SL")){

      opt_ind = which(tune_cost == min(tune_cost), arr.ind = TRUE)
      opt_ind1 = opt_ind[1]; opt_ind2 = opt_ind[2]

      full_model = sparse_group_fused_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                                   Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                                   lambda = lambda1, alpha = alpha, tau0 = tau0[opt_ind2],
                                                   nlambda = num_lambda1)

      penalty_parameter_sequence = full_model$penalty_parameter_sequence

      #opt_penalty_parameter = penalty_parameter_sequence[(opt_ind1 - 1) * length(lambda2) + opt_ind2,]
      opt_penalty_parameter = penalty_parameter_sequence[opt_ind1, ]


      model_info = list(intercept = full_model$interceptlist[[opt_ind1]],
                        beta = full_model$betalist[[opt_ind1]],
                        penalty_parameter_sequence = penalty_parameter_sequence_all,
                        opt_penalty_parameter = opt_penalty_parameter,
                        cv_error = tune_cost,
                        alpha = alpha, penalty = penalty, single_rule = FALSE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        Ybar = Ybarlist, Xbar = Xbarlist, Xsd = Xsdlist,
                        problem = problem)

    }
  }

  class(model_info) = "mp_cv"
  return(model_info)
}
