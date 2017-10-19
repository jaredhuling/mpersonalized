#' @title Meta-analysis/Multiple Outcomes for Personalized Medicine
#'
#' @details Assume the total number of studies is \eqn{K}. This function is aimed to solve meta-analysis/multiple outcomes problems for personalized medicine based on the following framework:
#' \deqn{ \min_{g_1,\dots,g_K} \frac{1}{2}\sum_{k=1}^K \sum_{i=1}^{n_k}\frac{|\hat{C}_k(X_{i})|}{\sum_{i=1}^{n_k}|\hat{C}_k(X_{i})|}\bigl [1\{\hat{C}_k(X_{i})>0\}-g_k(X_{i})\bigr]^2 + h(g_1,\dots,g_K)}
#' Here the regularization function \eqn{h} is of the form of a sum of sparse group lasso and fused lasso penalty
#' \deqn{h = (1-\alpha)\lambda_1\sqrt{q} \sum_{j=1}^p \|\boldsymbol{\beta_j}\|_2+\alpha \lambda_1  \sum_{j=1}^p \|\boldsymbol{\beta_j}\|_1+ \lambda_2 \sum_{j=1}^p \sum_{1\le a < b \le K}|\beta_{ja}-\beta_{jb}|}
#' where \eqn{\boldsymbol{\beta_j}=(\beta_{j1},\dots,\beta_{jK})}
#'
#' If we would like a unique rule to be obtained, we let \eqn{g_1 = \dots= g_K} and solve the following question instead
#' \deqn{\min_{g} \frac{1}{2}\sum_{k=1}^K \sum_{i=1}^{n_k}\frac{|\hat{C}_k(X_{i})|}{\sum_{i=1}^{n_k}|\hat{C}_k(X_{i})|}\bigl [1\{\hat{C}_k(X_{i})>0\}-g(X_{i})\bigr]^2 + h(g_1,\dots,g_K) + \lambda_{uni} \|\beta\|_1}
#'
#' If we want different rules, by setting \eqn{\lambda_1, \lambda_2, \alpha} differently, different penalties can be obtained.
#' \itemize{
#' \item If \eqn{\lambda_1, \lambda_2 \ne 0} and \eqn{\alpha \ne 0} or \eqn{1}, the penalty is "SGL+fused".
#' \item If \eqn{\lambda_1, \lambda_2 \ne 0} and \eqn{\alpha = 0}, the penalty is "GL+fused".
#' \item If \eqn{\lambda_1, \lambda_2 \ne 0} and \eqn{\alpha = 1}, the penalty is "lasso+fused".
#' \item If \eqn{\lambda_1 = 0, \lambda_2 \ne 0}, the penalty is "fused".
#' \item If \eqn{\lambda_1 \ne0, \lambda_2 = 0} and \eqn{\alpha \ne 0} or \eqn{1}, the penalty is "SGL".
#' \item If \eqn{\lambda_1 \ne0, \lambda_2 = 0} and \eqn{\alpha = 0}, the penalty is "GL".
#' \item If \eqn{\lambda_1 \ne0, \lambda_2 = 0} and \eqn{\alpha = 1}, the penalty is "lasso".
#' \item If \eqn{\lambda_1, \lambda_2 = 0}, there is no penalty.
#' }
#'
#' If we want unique rule,
#' \itemize{
#' \item If \eqn{\lambda \ne 0}, the penalty is "lasso".
#' \item If \eqn{\lambda = 0}, there is no penalty.
#' }
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
#' @import glmnet SGL Matrix
#'
#' @return an S3 object of class "mp", which contains the information of the fitted model. It could be supplied
#' to the predict function
#' @export

mpersonalized = function(problem = c("meta-analysis", "multiple outcomes"),
                         X, Trt, P,
                         Xlist, Ylist, Trtlist, Plist, typelist = NULL,
                         penalty = c("none", "lasso", "GL", "SGL", "fused",
                                     "lasso+fused", "GL+fused", "SGL+fused"),
                         lambda1 = NULL, lambda2 = NULL, unique_rule_lambda = NULL,
                         alpha = NULL, unique_rule = FALSE){

  penalty = match.arg(penalty)
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

    if (penalty != "lasso" & penalty != "none")
      stop("When unique rule = TRUE, the penalty must be lasso or none with default as none!")

    if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
      warning("When unique rule = TRUE, the value for lambda1, lambda2, alpha are ignored!")

    if (penalty == "lasso"){

      full_model = unique_rule_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                            Ybar = Ybar, Xbar = Xbar, Xsd = Xsd, lambda = unique_rule_lambda)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        lambda = full_model$lambda,
                        penalty = penalty, unique_rule = TRUE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    } else if (penalty == "none"){

      if (!is.null(unique_rule_lambda))
        warning("When unique rule = TRUE and penalty = none, the value for unique_rule_lambda are ignored!")

      full_model = unique_rule_linear_method(Conlist = Conlist, Xlist = Xlist)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty = penalty, unique_rule = TRUE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)
    }

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (!is.null(unique_rule_lambda))
        warning("When unique rule = FALSE, the value for unique_rule_lambda is ignored!")

    if (penalty == "none"){

      if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
        warning("When penalty = none, the value for lambda1, lambda2, alpha are ignored!")

      full_model = linear_method(Conlist = Conlist, Xlist = Xlist)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty = penalty, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

    } else if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")) {

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

      full_model = meta_method(modelYlist = modelYlist, modelXlist = modelXlist,
                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                               lambda1 = lambda1, lambda2 = lambda2, alpha = alpha)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        iterslist = full_model$iterslist,
                        lambda1 = lambda1, lambda2 = lambda2,
                        alpha = alpha, penalty = penalty, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)

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


      full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                             Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                             lambda = lambda1, alpha = alpha)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        lambda1 = full_model$lambda, lambda2 = lambda2,
                        alpha = alpha, penalty = penalty, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist)
    }
  }

  class(model_info) = "mp"
  return(model_info)

}
