#' @title Meta-analysis/Multiple Outcomes for Personalized Medicine
#'
#' @details Assume the total number of studies is K. This function implements meta-analysis for personalized medicine based on the following framework:
#' \deqn{ \min_{g_1,\dots,g_K} \frac{1}{2}\sum_{k=1}^K \sum_{i=1}^{n_k}\frac{|\hat{C}_k(X_{i})|}{\sum_{i=1}^{n_k}|\hat{C}_k(X_{i})|}\bigl [1\{\hat{C}_k(X_{i})>0\}-g_k(X_{i})\bigr]^2 + h(g_1,\dots,g_K)}
#' Here the regularization function \eqn{h} is of the form of a sum of sparse group lasso and fused lasso penalty
#' \deqn{h = (1-\alpha)\lambda_1\sqrt{q} \sum_{j=1}^p \|\boldsymbol{\beta_j}\|_2+\alpha \lambda_1  \sum_{j=1}^p \|\boldsymbol{\beta_j}\|_1+ \lambda_2 \sum_{j=1}^p \sum_{1\le a < b \le K}|\beta_{ja}-\beta_{jb}|}
#' where \eqn{\boldsymbol{\beta_j}=(\beta_{j1},\dots,\beta_{jK})}
#'
#' If we would like a unique rule to be obtained, we should let \eqn{g_1 = \dots= g_K} and
#' \deqn{h = \lambda_{uni} \|\beta\|_1}
#'
#' By setting \eqn{\lambda_1, \lambda_2, \alpha} differently, different models could be obtained.
#' \itemize{
#' \item If \eqn{\lambda_1,\lambda2 \ne 0} and \eqn{alpha \ne 0} or \eqn{1}, it is sparse group lasso + fused model.
#' \item If \eqn{\lambda_2 = 0} and \eqn{\alpha \ne 0} or \eqn{1}, a sparse group lasso model is fitted.
#' \item If \eqn{\lambda_2 = 0} and \eqn{\alpha = 0}, a group lasso model is fitted.
#' \item If \eqn{\lambda_2 = 0} and \eqn{\alpha = 1}, a lasso model is fitted.
#' \item If \eqn{\lambda_1, \lambda_2 = 0}, a linear model is fitted.
#' }
#'
#' @param problem a character specifiy whether you want to solve "meta-analysis" or "multiple outcomes" problem
#' @param Xlist a list object with \eqn{k}th element denoting the covariate matrix of study k
#' @param Ylist a list object with \eqn{k}th element denoting the response vector of study k
#' @param Trtlist  a list object with \eqn{k}th element denoting the treatment vector of study k (coded as 0 or 1)
#' @param typlelist a list object with \eqn{k}th element denoting the type of response in study k, can be continuous or binary, default is continuous
#' @param model the model to be used for the above framework, can be linear, meta-analysis,
#' sparse group lasso, group lasso or lasso
#' @param lambda1 lambda1 in the framework above
#' @param lambda2 lambda2 in the framework above
#' @param alpha alpha in the framework above
#' @param unique_rule_lambda \eqn{\lambda_{uni}} when unique treatment rule is required
#' @param unique_rule a logical value, whether a unique treatment rule is required
#' @import glmnet SGL Matrix
#'
#' @return an S3 object of class "mp", which contains the information of the fitted model. It could be supplied
#' to the predict function
#' @export

MetaPersonalized = function(problem = c("meta-analysis", "multiple outcomes"),
                            X, Trt, P,
                            Xlist, Ylist, Trtlist, Plist, typelist = NULL,
                            penalty = c("linear", "lasso", "GL", "SGL", "fused",
                                      "lasso+fused", "GL+fused", "SGL+fused"),
                            lambda1 = NULL, lambda2 = NULL, unique_rule_lambda = NULL,
                            alpha = NULL, unique_rule = FALSE){

  model = match.arg(model)
  problem = match.arg(problem)

  if (problem == "multiple outcomes"){

    if (is.null(X) | is.null(Ylist) | is.null(Trt) | is.null(P))
      stop("For multiple outcomes, X, Ylist, Trt, P need to be supplied!")

    q = length(Ylist)
    Xlist = replicate(q, X, simplify = FALSE)
    Trtlist = replicate(q, Trt, simplify = FALSE)
    Plist = replicate(q, P, simplify = FALSE)

  } else if (problem == "meta-analysis"){

    if (is.null(Xlist) | is.null(Ylist) | is.null(Trtlist) | is.null(Plist))
      stop("For meta-analysis, Xlist, Ylist, Trtlist, Plist need to be supplied!")
  }

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

  if (unique_rule == TRUE){

    Ybar = standardized_data$Ybar
    Xbar = standardized_data$Xbar
    Xsd = standardized_data$Xsd

    if (model != "lasso" & model != "linear")
      stop("When unique rule = TRUE, the model must be lasso or linear with default as linear!")

    if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
      warning("When unique rule = TRUE, the value for lambda1, lambda2, alpha are ignored!")

    if (model == "lasso"){

      full_model = unique_rule_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                            Ybar = Ybar, Xbar = Xbar, Xsd = Xsd, lambda = unique_rule_lambda)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        lambda = full_model$lambda,
                        model = model, unique_rule = TRUE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    } else if (model == "linear"){

      if (!is.null(unique_rule_lambda))
        warning("When unique rule = TRUE and model = linear, the value for unique_rule_lambda are ignored!")

      full_model = unique_rule_linear_method(Conlist = Conlist, Xlist = Xlist)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        model = model, unique_rule = TRUE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)
    }

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (!is.null(unique_rule_lambda))
        warning("When unique rule = FALSE, the value for unique_rule_lambda is ignored!")

    if (model == "linear"){

      if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
        warning("When model = linear, the value for lambda1, lambda2, alpha are ignored!")

      full_model = linear_method(Conlist = Conlist, Xlist = Xlist)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        model = model, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    } else if (model %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")) {

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

      full_model = meta_method(modelYlist = modelYlist, modelXlist = modelXlist,
                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                               lambda1 = lambda1, lambda2 = lambda2, alpha = alpha)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        iterslist = full_model$iterslist,
                        lambda1 = lambda1, lambda2 = lambda2,
                        alpha = alpha, model= model, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

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


      full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                             Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                             lambda = lambda1, alpha = alpha)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        lambda1 = full_model$lambda, lambda2 = lambda2,
                        alpha = alpha, model = model, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)
    }
  }

  class(model_info) = "mp"
  return(model_info)

}
