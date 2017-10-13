#' Meta-analysis/Multiple Outcomes for Personalized Medicine with Cross Validation
#' @param Xlist a list object with \eqn{k}th element denoting the covariate matrix of study k
#' @param Ylist a list object with \eqn{k}th element denoting the response vector of study k
#' @param Trtlist  a list object with \eqn{k}th element denoting the treatment vector of study k (coded as 0 or 1)
#' @param typlelist a list object with \eqn{k}th element denoting the type of response in study k, can be continuous or binary, default is continuous
#' @param model the model to be used for the above framework, can be meta-analysis,
#' sparse group lasso, group lasso or lasso(linear does not need tuning)
#' @param lambda1 lambda1 in the framework above
#' @param lambda2 lambda2 in the framework above
#' @param alpha alpha in the framework above
#' @param unique_rule_lambda \eqn{\lambda_{uni}} when unique treatment rule is required
#' @param unique_rule a logical value, whether a unique treatment rule is required
#' @param cv_folds number of folds needed for cross-validation, default is 5
#' @import glmnet SGL caret Matrix
#' @return an S3 object of class "mp_cv", which contains the information of the model with the best fitted lambda. It can be supplied to the predict function.
#' @export

MetaPersonalized_cv = function(Xlist, Ylist, Trtlist, Plist, typelist = NULL,
                               model = c("lasso", "GL", "SGL", "fused",
                                         "lasso+fused", "GL+fused", "SGL+fused"),
                              lambda1 = NULL, lambda2 = NULL, unique_rule_lambda = NULL,
                              alpha = NULL, unique_rule = FALSE, cv_folds = 5){

  model = match.arg(model)

  q = length(Xlist)
  p = dim(Xlist[[1]])[2]

  if(is.null(typelist))
    typelist = replicate(q, list("continuous"))

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

    if (model != "lasso")
      stop("When unique rule is required, the model must be lasso!(for linear model, use 'MetaPersonalized' instead.")

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

    if (model %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")) {

      if (model != "fused"){
        if (is.null(lambda1) | is.null(lambda2))
          stop("When model = lasso+fused/GL+fused/SGL+fused, both values of lambda1 and lambda2 must be supplied!")
        if (sum(lambda1 == 0) > 0 | sum(lambda2 == 0) > 0)
          stop("When model = lasso+fused/GL+fused/SGL+fused, do not supply 0 in lambda1 and lambda2!")
      }

      if (model == "fused"){
        if (is.null(lambda2))
          stop("When model = fused, lambda2 must be supplied!")
        if (sum(lambda2 == 0) > 0)
          stop("When model = fused, do not supply 0 in lambda2!")
      }

      if (model == "fused"){
        if (!is.null(alpha)){
          warning("When model = fused, values of alpha is ignored!")
          alpha = NULL
        }

        if (!is.null(lambda1))
          if (sum(lambda1 != 0) > 0)
            warning("When model = fused, value of lambda1 is automatically set to be 0!")

        lambda1 = 0
      }

      if (model == "lasso+fused"){

        if (!is.null(alpha))
          if(alpha != 1)
            warning("When model = lasso+fused, alpha is automatically set to be 1!")

        alpha = 1

      } else if (model == "GL+fused"){

        if (!is.null(alpha))
          if(alpha != 0)
            warning("When model = GL+fused, alpha is automatically set to be 0!")

        alpha = 0

      } else if (model == "SGL+fused"){

        if (!is.null(alpha))
          if (alpha == 0 | alpha == 1){
            warning("When model = SGL+fused, alpha cannot be set as 0 or 1, and default is 0.95!")
            alpha = 0.95
          } else if (is.null(alpha)){
            alpha = 0.95
          }
      }
    } else if (model %in% c("lasso", "GL", "SGL")){


      if (!is.null(lambda2)){
        if (sum(lambda2 != 0) > 0){
          warning("When model = lasso/GL/SGL, the value for lambda2 is ignored and automatically set to be 0!")
          lambda2 = 0
        }
      } else lambda2 = 0

      if (model == "lasso"){

        if (!is.null(alpha))
          if (alpha != 1)
            warning("When model = lasso, alpha is automatically set to be 1!")

        alpha = 1

      } else if (model == "GL"){

        if (!is.null(alpha))
          if (alpha != 0)
            warning("When model = GL, alpha is automatically set to be 0!")

        alpha = 0

      } else if (model == "SGL"){

        if (!is.null(alpha)){
          if (alpha == 0 | alpha == 1){
            warning("When model = SGL, alpha cannot be set as 0 or 1, and default is 0.95!")
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
    if (model %in% c("lasso", "GL", "SGL")){
      tune_cost = numeric(length(lambda1))
    } else if (model %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){
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

      if (model %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){

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

      } else if (model %in% c("lasso", "SGL", "GL")) {

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
                      opt_unique_lambda = unique_rule_lambda[opt_ind], model = "lasso", unique_rule = TRUE,
                      number_covariates = p, number_studies = q,
                      Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (model %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){

      opt_ind = which(tune_cost == min(tune_cost), arr.ind = TRUE)
      opt_ind1 = opt_ind[1]; opt_ind2 = opt_ind[2]
      full_model = meta_method(modelYlist = modelYlist, modelXlist = modelXlist,
                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                               lambda1 = lambda1[opt_ind1], lambda2 = lambda2[opt_ind2], alpha = alpha)

      model_info = list(intercept = full_model$interceptlist[[1]], beta = full_model$betalist[[1]],
                        iters = full_model$iterslist[[1]],
                        lambda1 = lambda1, lambda2 = lambda2,
                        opt_lambda = list(opt_lambda1 = lambda1[opt_ind1], opt_lambda2 = lambda2[opt_ind2]),
                        alpha = alpha, model = model, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    } else if (model %in% c("lasso", "SGL", "GL")){

      full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                             Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                             lambda = lambda1, alpha = alpha)

      opt_ind = which.min(tune_cost)
      model_info = list(intercept = full_model$interceptlist[[opt_ind]], beta = full_model$betalist[[opt_ind]],
                        lambda1 = full_model$lambda, lambda2 = 0,
                        opt_lambda1 = lambda1[opt_ind], opt_lambda2 = 0,
                        alpha = alpha, model = model, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    }
  }

  class(model_info) = "mp_cv"
  return(model_info)
}
