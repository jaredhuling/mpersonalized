#' @title Meta-analysis/Multiple Outcomes for Personalized Medicine with Cross Validation
#'
#' @details This function implments basically implements \code{mpersonalized} but use cross validatation for the tuning of penalty parameter.
#'  The optimal penalty parameter is selected by minimizing \deqn{\sum_{i=1}^{n_k}\frac{|\hat{C}_k(X_{i})|}{\sum_{i=1}^{n_k}|\hat{C}_k(X_{i})|}\bigl [1\{\hat{C}_k(X_{i})>0\}-g_k(X_{i})\bigr]^2}
#'  where \eqn{\hat{C}_k(X_{i})} in the leave-out fold is separately estimated from the training set.
#'
#' @param problem a character specifiy whether you want to solve "meta-analysis" or "multiple outcomes" problem. For "meta-analysis" problem,
#'  the user should supply \code{Xlist}, \code{Ylist}, \code{Trtlist} and \code{Plist}. For "multiple outcomes" problem,
#'  the user should supply \code{X}, \code{Ylist}, \code{Trt} and \code{P}.
#' @param X the covariate matrix that should be supplied when the problem is "multiple outcomes" with rows indicating subjects and columns indicating covariates.
#' @param Trt the treatment vector that should be supplied when the problem is "multiple outcomes". It should be coded as 0 or 1.
#' @param P the propensity score vector when the problem is "multiple outcomes".
#' @param Xlist a list object with \eqn{k}th element denoting the covariate matrix of study \eqn{k}. This should be supplied when the problem is
#'  "meta-analysis".
#' @param Ylist When the problem is "meta-analysis", \code{Ylist} should be a list object with \eqn{k}th element denoting the response vector of study \eqn{k}. When the
#'  problem is "multiple outcomes", \code{Ylist} should be a list object with \eqn{k}th element denoting the \eqn{k}th outcome.
#' @param Trtlist  a list object with \eqn{k}th element denoting the treatment vector of study \eqn{k} (coded as 0 or 1). This should be supplied when the problem is
#'  "meta-analysis".
#' @param Plist a list object with \eqn{k}the element denoting the propensity score vector of study \eqn{k}.
#' @param typlelist a list object with \eqn{k}th element denoting the type of response corresponding to the \eqn{k}th element in the list \code{Ylist}.
#'  Each element should be "continuous" or "binary".
#' @param penalty For different rules, the penalty could be "none", "lasso", "GL", "SGL", "fused",
#'  "lasso+fused", "GL+fused", "SGL+fused". For unique rule, the penalty could be "none" or "lasso".
#' @param lambda1 lambda1 supplied in the framework when different rules are used.
#' @param lambda2 lambda2 supplied in the framework when different rules are used.
#' @param alpha alpha in the framework when different rules are used.
#' @param unique_rule_lambda \eqn{\lambda} when unique rule is used.
#' @param unique_rule a logical value, whether a unique treatment rule is required
#' @param cv_folds number of folds needed for cross-validation, default is 5
#' @import glmnet SGL caret Matrix
#' @return an S3 object of class "mp_cv", which contains the information of the model with the best fitted lambda. It can be supplied to the predict function.
#' @export

mpersonalized_cv = function(problem = c("meta-analysis", "multiple outcomes"),
                            X, Trt, P = NULL,
                            Xlist, Ylist, Trtlist, Plist = replicate(length(Xlist), NULL, simplify = FALSE),
                            typelist = replicate(length(Xlist), "continuous", simplify = FALSE),
                            penalty = c("lasso", "GL", "SGL", "fused",
                                      "lasso+fused", "GL+fused", "SGL+fused"),
                            lambda1 = NULL, lambda2 = NULL, unique_rule_lambda = NULL,
                            alpha = NULL, unique_rule = FALSE, cv_folds = 5){

  penalty = match.arg(penalty)
  problem = match.arg(problem)

  if (problem == "multiple outcomes"){

    if (missing(X) | missing(Ylist) | missing(Trt))
      stop("For multiple outcomes, X, Ylist, Trt need to be supplied!")

    q = length(Ylist)
    Xlist = replicate(q, X, simplify = FALSE)
    Trtlist = replicate(q, Trt, simplify = FALSE)
    Plist = replicate(q, P, simplify = FALSE)

  } else if (problem == "meta-analysis"){

    if (missing(Xlist) | missing(Ylist) | missing(Trtlist))
      stop("For meta-analysis, Xlist, Ylist, Trtlist need to be supplied!")
  }

  q = length(Xlist)
  p = dim(Xlist[[1]])[2]


  #construct contrast for the data
  Conlist = mapply(contrast_builder, X = Xlist, Y = Ylist,
                   ori_Trt = Trtlist, P = Plist, type = typelist,
                   MoreArgs = list(response_model = "lasso"), SIMPLIFY = FALSE)

  standardized_data = contrast_standardize(Conlist = Conlist, Xlist = Xlist,
                                           unique_rule = unique_rule)
  modelYlist = standardized_data$modelYlist
  modelXlist = standardized_data$modelXlist

  #check model is correct and set up the value for lambdas if not provided
  if (unique_rule == TRUE){

    Ybar = standardized_data$Ybar
    Xbar = standardized_data$Xbar
    Xsd = standardized_data$Xsd

    if (penalty != "lasso")
      stop("When unique rule is required, the penalty must be lasso!(for no penalty, use 'mersonalized' instead.")

    if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
      warning("When unique rule = TRUE, the value for lambda1, lambda2, alpha are ignored!")

    if (is.null(unique_rule_lambda)){
      full_model = unique_rule_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                            Ybar = Ybar, Xbar = Xbar, Xsd = Xsd, lambda = unique_rule_lambda)
      unique_rule_lambda = full_model$lambda
    }

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (!is.null(unique_rule_lambda))
      warning("When unique rule = FALSE, the value for unique_rule_lambda is ignored!")

    if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")) {

      if (penalty != "fused"){
        if (is.null(lambda1) | is.null(lambda2))
          stop("When penalty = lasso+fused/GL+fused/SGL+fused, both values of lambda1 and lambda2 must be supplied!")
        if (sum(lambda1 == 0) > 0 | sum(lambda2 == 0) > 0)
          stop("When penalty = lasso+fused/GL+fused/SGL+fused, do not supply 0 in lambda1 and lambda2!")
      }

      if (penalty == "fused"){
        if (is.null(lambda2))
          stop("When penalty = fused, lambda2 must be supplied!")
        if (sum(lambda2 == 0) > 0)
          stop("When penalty = fused, do not supply 0 in lambda2!")
      }

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

        if (!is.null(alpha))
          if (alpha == 0 | alpha == 1){
            warning("When penalty = SGL+fused, alpha cannot be set as 0 or 1, and default is 0.95!")
            alpha = 0.95
          } else if (is.null(alpha)){
            alpha = 0.95
          }
      }
    } else if (penalty %in% c("lasso", "GL", "SGL")){


      if (!is.null(lambda2)){
        if (sum(lambda2 != 0) > 0){
          warning("When penalty = lasso/GL/SGL, the value for lambda2 is ignored and automatically set to be 0!")
          lambda2 = 0
        }
      } else lambda2 = 0

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
        } else {
          alpha = 0.95
        }
      }

      if (is.null(lambda1)){
        full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                               lambda = lambda1, alpha = alpha)
        lambda1 = full_model$lambda
      }
    }
  }

  #build cross-validation folds
  folds_index = lapply(Xlist, function(x) createFolds(1:dim(x)[1], k = cv_folds))

  #determine the dimension for the cost for tuning penalty parameters
  if (unique_rule == TRUE){
    tune_cost = numeric(length(unique_rule_lambda))
  } else {
    if (penalty %in% c("lasso", "GL", "SGL")){
      tune_cost = numeric(length(lambda1))
    } else if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){
      tune_cost = matrix(0, nrow = length(lambda1), ncol = length(lambda2))
    }
  }

  #carry out method for each cross validation fold
  for (k in 1:cv_folds){

    cv_Xlist = vector("list", q); left_Xlist = vector("list", q)
    cv_Ylist = vector("list", q); left_Ylist = vector("list", q)
    cv_Trtlist = vector("list", q); left_Trtlist = vector("list", q)

    for (j in 1:q){
      cv_Xlist[[j]] = Xlist[[j]][-folds_index[[j]][[k]],]
      cv_Ylist[[j]] = Ylist[[j]][-folds_index[[j]][[k]]]
      cv_Trtlist[[j]] = Trtlist[[j]][-folds_index[[j]][[k]]]

      left_Xlist[[j]] = Xlist[[j]][folds_index[[j]][[k]],]
      left_Ylist[[j]] = Ylist[[j]][folds_index[[j]][[k]]]
      left_Trtlist[[j]] = Trtlist[[j]][folds_index[[j]][[k]]]
    }

    cv_Conlist = mapply(contrast_builder, X = cv_Xlist, Y = cv_Ylist,
                     ori_Trt = cv_Trtlist, P = Plist, type = typelist,
                     MoreArgs = list(response_model = "lasso"), SIMPLIFY = FALSE)
    left_Conlist = mapply(contrast_builder, X = left_Xlist, Y = left_Ylist,
                      ori_Trt = left_Trtlist, P = Plist, type = typelist,
                      MoreArgs = list(response_model = "lasso"), SIMPLIFY = FALSE)

    cv_standardized_data = contrast_standardize(Conlist = cv_Conlist, Xlist = cv_Xlist, unique_rule = unique_rule)
    cv_modelYlist = cv_standardized_data$modelYlist
    cv_modelXlist = cv_standardized_data$modelXlist

    #transform contrast into binary data for the left out fold
    left_sConlist = lapply(left_Conlist, function(y) as.numeric(y > 0))
    left_Wlist = lapply(left_Conlist, abs)

    left_dataWlist = lapply(left_Wlist, sum)
    left_adj_Wlist = mapply(function(w, dataw) w / dataw, w = left_Wlist,
                       dataw = left_dataWlist, SIMPLIFY = FALSE)

    if (unique_rule == TRUE){
      cv_Ybar = cv_standardized_data$Ybar
      cv_Xbar = cv_standardized_data$Xbar
      cv_Xsd = cv_standardized_data$Xsd

      cv_model = unique_rule_lasso_method(modelYlist = cv_modelYlist, modelXlist = cv_modelXlist,
                                          Ybar = cv_Ybar, Xbar = cv_Xbar, Xsd = cv_Xsd, lambda = unique_rule_lambda)
      cv_interceptlist = cv_model$interceptlist
      cv_betalist = cv_model$betalist

      for (ind in 1:length(unique_rule_lambda))
        tune_cost[ind] = tune_cost[ind] + sum(unlist(mapply(function(w, y, x) sum(w * (y - cv_interceptlist[[ind]] - x %*% cv_betalist[[ind]]) ^ 2),
                                                            w = left_adj_Wlist, y = left_sConlist, x = left_Xlist, SIMPLIFY = FALSE)))

    } else {
      cv_Ybarlist = cv_standardized_data$Ybarlist
      cv_Xbarlist = cv_standardized_data$Xbarlist
      cv_Xsdlist = cv_standardized_data$Xsdlist

      if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){

        cv_model = meta_method(modelYlist = cv_modelYlist, modelXlist = cv_modelXlist,
                               Ybarlist = cv_Ybarlist, Xbarlist = cv_Xbarlist, Xsdlist = cv_Xsdlist,
                               lambda1 = lambda1, lambda2 = lambda2, alpha = alpha)
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

      }
    }
  }

  #for the whole data set

  if (unique_rule == TRUE){

    Ybar = standardized_data$Ybar
    Xbar = standardized_data$Xbar
    Xsd = standardized_data$Xsd

    full_model = unique_rule_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                          Ybar = Ybar, Xbar = Xbar, Xsd = Xsd, lambda = unique_rule_lambda)

    opt_ind = which.min(tune_cost)
    model_info = list(intercept = full_model$interceptlist[[opt_ind]], beta = full_model$betalist[[opt_ind]],
                      unique_lambda = full_model$lambda,
                      opt_unique_lambda = unique_rule_lambda[opt_ind], penalty = "lasso", unique_rule = TRUE,
                      number_covariates = p, number_studies = q,
                      Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){

      opt_ind = which(tune_cost == min(tune_cost), arr.ind = TRUE)
      opt_ind1 = opt_ind[1]; opt_ind2 = opt_ind[2]
      full_model = meta_method(modelYlist = modelYlist, modelXlist = modelXlist,
                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                               lambda1 = lambda1[opt_ind1], lambda2 = lambda2[opt_ind2], alpha = alpha)

      model_info = list(intercept = full_model$interceptlist[[1]], beta = full_model$betalist[[1]],
                        iters = full_model$iterslist[[1]],
                        lambda1 = lambda1, lambda2 = lambda2,
                        opt_lambda = list(opt_lambda1 = lambda1[opt_ind1], opt_lambda2 = lambda2[opt_ind2]),
                        alpha = alpha, penalty = penalty, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    } else if (penalty %in% c("lasso", "SGL", "GL")){

      full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                             Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                             lambda = lambda1, alpha = alpha)

      opt_ind = which.min(tune_cost)
      model_info = list(intercept = full_model$interceptlist[[opt_ind]], beta = full_model$betalist[[opt_ind]],
                        lambda1 = full_model$lambda, lambda2 = 0,
                        opt_lambda1 = lambda1[opt_ind], opt_lambda2 = 0,
                        alpha = alpha, penalty = penalty, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    }
  }

  class(model_info) = "mp_cv"
  return(model_info)
}
