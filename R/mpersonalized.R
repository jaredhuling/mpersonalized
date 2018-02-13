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
#' @param problem a character specifiy whether the user want to solve "meta-analysis" or "multiple outcomes" problem. For "meta-analysis" problem,
#'  the user should supply \code{Xlist}, \code{Ylist}, \code{Trtlist} and \code{Plist}. For "multiple outcomes" problem,
#'  the user should supply \code{X}, \code{Ylist}, \code{Trt} and \code{P}.
#' @param X the covariate matrix that should be supplied when the problem is "multiple outcomes" with rows indicating subjects and columns indicating covariates.
#' @param Trt the treatment vector that should be supplied when the problem is "multiple outcomes". It should be coded as 0 or 1.
#' @param P the propensity score vector when the problem is "multiple outcomes". If not supplied, then study is treated as randomzied trial and the propensity
#'  score is estimated as the proportion of 1's in Trt.
#' @param Xlist a list object with \eqn{k}th element denoting the covariate matrix of study \eqn{k}. This should be supplied when the problem is
#'  "meta-analysis".
#' @param Ylist When the problem is "meta-analysis", \code{Ylist} should be a list object with \eqn{k}th element denoting the response vector of study \eqn{k}. When the
#'  problem is "multiple outcomes", \code{Ylist} should be a list object with \eqn{k}th element denoting the \eqn{k}th outcome.
#' @param Trtlist  a list object with \eqn{k}th element denoting the treatment vector of study \eqn{k} (coded as 0 or 1). This should be supplied when the problem is
#'  "meta-analysis".
#' @param Plist a list object with \eqn{k}the element denoting the propensity score vector of study \eqn{k}. If not supplied, then
#'  each study is treated as randomized trial and the corresponding propensity score is estimated as the proportion of 1's in Trt.
#' @param typlelist a list object with \eqn{k}th element denoting the type of response corresponding to the \eqn{k}th element in the list \code{Ylist}.
#'  Each element should be "continuous" or "binary".
#' @param penalty For different rules, the penalty could be "none", "lasso", "GL", "SGL", "fused",
#'  "lasso+fused", "GL+fused", "SGL+fused". For unique rule, the penalty could be "none" or "lasso".
#' @param lambda1 lambda1 supplied in the framework when different rules are used.
#' @param lambda2 lambda2 supplied in the framework when different rules are used.
#' @param alpha alpha in the framework when different rules are used.
#' @param unique_rule_lambda \eqn{\lambda} when unique rule is used.
#' @param unique_rule a logical value, whether a unique treatment rule is required
#' @param admm_control a list of parameters which control the admm algorithm. In admm_control, the following parameters can be supplied:
#' abs.tol, absolute tolerance; rel.tol, relative tolerance; maxit, maximum number of iterations; rho, Lagrangian parameter.
#' @param contrast_builder_control a list of parameters which control the contrast building process. In contrast_builder_control,
#' the following parameters could be supplied: response_model, this could be "lasso" or "linear"; contrast_builder_folds,
#' the number of folds used in cross validation when response_model = "lasso".
#' @param num_lambda1 length of the lambda1 sequence and default to be 10 if lambda1 is not provided
#' @param num_lambda2 length of the lambda2 sequence and default to be 10 if lambda2 is not provided
#' @param num_unique_rule_lambda length of the unique_rule_lambda sequence and default to be 50 if unique_rule_lambda is not provided
#' @import glmnet SGL Matrix
#'
#' @return an S3 object of class "mp", which contains the information of the fitted model. It could be supplied
#' to the predict function
#' @export

mpersonalized = function(problem = c("meta-analysis", "multiple outcomes"),
                         X, Trt, P = NULL,
                         Xlist, Ylist, Trtlist, Plist = replicate(length(Xlist), NULL, simplify = FALSE),
                         typelist = replicate(length(Xlist), "continuous", simplify = FALSE),
                         penalty = c("none", "lasso", "GL", "SGL", "fused",
                                     "lasso+fused", "GL+fused", "SGL+fused"),
                         lambda1 = NULL, lambda2 = NULL, unique_rule_lambda = NULL,
                         num_lambda1 = ifelse(!is.null(lambda1), length(lambda1),10),
                         num_lambda2 = ifelse(!is.null(lambda2), length(lambda2),10),
                         num_unique_rule_lambda = ifelse(!is.null(unique_rule_lambda), length(unique_rule_lambda), 50),
                         alpha = NULL, unique_rule = FALSE,
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
    #create default value for the propensity score
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
    #create defauly value for the propensity score
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
  #have to use for loop because we need to use do.call to input the
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

      if (is.null(unique_rule_lambda)){
        lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                         penalty = penalty, unique_rule = unique_rule,
                                         num_unique_rule_lambda = num_unique_rule_lambda)

        unique_rule_lambda = lambda_default$unique_rule_lambda
      }

      full_model = unique_rule_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                            Ybar = Ybar, Xbar = Xbar, Xsd = Xsd, lambda = unique_rule_lambda)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        unique_rule_lambda = unique_rule_lambda,
                        penalty = penalty, unique_rule = TRUE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist)

    } else if (penalty == "none"){

      if (!is.null(unique_rule_lambda))
        warning("When unique rule = TRUE and penalty = none, the value for unique_rule_lambda are ignored!")

      full_model = unique_rule_linear_method(Conlist = Conlist, Xlist = Xlist)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty = penalty, unique_rule = TRUE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist)
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
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist)

    } else if (penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")) {

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
                                         penalty = penalty, unique_rule = unique_rule, alpha = alpha,
                                         num_lambda1 = num_lambda1, num_lambda2 = num_lambda2)

        if (is.null(lambda1))
          lambda1 = lambda_default$lambda1

        if (is.null(lambda2))
          lambda2 = lambda_default$lambda2
      }

      full_model = meta_method(modelYlist = modelYlist, modelXlist = modelXlist,
                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                               lambda1 = lambda1, lambda2 = lambda2, alpha = alpha, admm_control = admm_control)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        iterslist = full_model$iterslist,
                        lambda1 = lambda1, lambda2 = lambda2,
                        alpha = alpha, penalty = penalty, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist)

    } else if (penalty %in% c("lasso", "GL", "SGL")){

      if (!is.null(lambda2)){
        if (sum(lambda2 != 0) > 0){
          warning("When penalty = lasso/GL/SGL, the value for lambda2 is ignored and automatically set to be 0!")
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
      }

      if (is.null(lambda1)){
        lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                         penalty = penalty, unique_rule = unique_rule, alpha = alpha,
                                         num_lambda1 = num_lambda1)

        lambda1 = lambda_default$lambda1
      }


      full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                             Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                             lambda = lambda1, alpha = alpha)

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        lambda1 = lambda1, lambda2 = lambda2,
                        alpha = alpha, penalty = penalty, unique_rule = FALSE,
                        number_covariates = p, number_studies = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist)
    }
  }

  class(model_info) = "mp"
  return(model_info)

}
