#' @title A General Framework to Solve Personalized Medicine in the Settings of Meta-analysis/Multiple Outcomes
#'
#' @description This function solves the personalized medicine problem by extending the contrast classificaiton
#' framework (Zhang, 2012). By adding proper penalty to the original classification loss, variable selection could be
#' implemented and estimation efficiency could be improved. Computation algorithm differs based on the penalty function employed,
#' but mainly through ADMM algorithm, glmnet package and SGL package. This function is also flexible enough to let
#' user choose whether different classification rules or a single rule should be estimated for multiple studies/outcomes.
#'
#' @details Assume the total number of studies/outcomes is \eqn{K} and we denote the contrast
#' estimator for the \eqn{k}th study/outcome as \eqn{\hat{C}_k} and the corresponding recommendation
#' rule as \eqn{g_k}.
#'
#' If we want different rules for each study/outcome, this function solves meta-analysis/multiple outcomes
#' problems for personalized medicine based on the framework
#' \deqn{ \min_{g_1,\dots,g_K} \frac{1}{2}\sum_{k=1}^K \sum_{i=1}^{n_k}\frac{|\hat{C}_k(X_{i})|}{\sum_{i=1}^{n_k}|\hat{C}_k(X_{i})|}\bigl [1\{\hat{C}_k(X_{i})>0\}-g_k(X_{i})\bigr]^2 + h(g_1,\dots,g_K)}
#' Here the regularization function \eqn{h} is of the form of a sum of sparse group lasso and fused lasso penalty
#' \deqn{h = (1-\alpha)\lambda_1\sqrt{q} \sum_{j=1}^p \|\boldsymbol{\beta_j}\|_2+\alpha \lambda_1  \sum_{j=1}^p \|\boldsymbol{\beta_j}\|_1+ \lambda_2 \sum_{j=1}^p \sum_{1\le a < b \le K}|\beta_{ja}-\beta_{jb}|}
#' where \eqn{\boldsymbol{\beta_j}=(\beta_{j1},\dots,\beta_{jK})}
#'
#' By setting \eqn{\lambda_1, \lambda_2, \alpha} differently, different penalties can be obtained.
#' \itemize{
#' \item If \eqn{\lambda_1, \lambda_2 \ne 0} and \eqn{\alpha \ne 0} or \eqn{1}, the penalty is "SGL+fused".
#' \item If \eqn{\lambda_1, \lambda_2 \ne 0} and \eqn{\alpha = 0}, the penalty is "GL+fused".
#' \item If \eqn{\lambda_1, \lambda_2 \ne 0} and \eqn{\alpha = 1}, the penalty is "lasso+fused".
#' \item If \eqn{\lambda_1 = 0, \lambda_2 \ne 0}, the penalty is "fused".
#' \item If \eqn{\lambda_1 \ne0, \lambda_2 = 0} and \eqn{\alpha \ne 0} or \eqn{1}, the penalty is "SGL".
#' \item If \eqn{\lambda_1 \ne0, \lambda_2 = 0} and \eqn{\alpha = 0}, the penalty is "GL".
#' \item If \eqn{\lambda_1 \ne0, \lambda_2 = 0} and \eqn{\alpha = 1}, the penalty is "lasso".
#' \item If \eqn{\lambda_1, \lambda_2 = 0}, there is no penalty.
#' \item If \eqn{\tau_0} is very large for \code{penalty = "SGL+SL"}, all coefficients across all studies will be equal.
#' \item If \eqn{\tau_0 = 0} for \code{penalty = "SGL+SL"}, an error will be returned. However, if \eqn{\tau_0} is close to 0,
#' more heterogeneity across studies will be preferred
#' }
#'
#' On the other hand, if we would like to fit a single rule for all studies/outcomes, we let \eqn{g_1 = \dots= g_K} and
#' solve the following problem instead
#' \deqn{\min_{g} \frac{1}{2}\sum_{k=1}^K \sum_{i=1}^{n_k}\frac{|\hat{C}_k(X_{i})|}{\sum_{i=1}^{n_k}|\hat{C}_k(X_{i})|}\bigl [1\{\hat{C}_k(X_{i})>0\}-g(X_{i})\bigr]^2 + h(g_1,\dots,g_K) + \lambda_{single} \|\beta\|_1}
#'
#' Depending on the value of \eqn{\lambda_{single}}
#' \itemize{
#' \item If \eqn{\lambda_{single} \ne 0}, the penalty is "lasso".
#' \item If \eqn{\lambda_{single} = 0}, there is no penalty.
#' }
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
#' @param penalty For different rules, the penalty could be "none", "lasso", "GL", "SGL", "fused",
#' "lasso+fused", "GL+fused", "SGL+fused", or "SGL+SL". For single rule, the penalty could be "none" or "lasso".
#' User should always input \code{penalty} and then supply correponding penalty parameters sequence
#' if needed. Default option is "none".
#' @param lambda1 \eqn{\lambda_1} in the framework of different rules. If not supplied, a default
#' sequence will be computed.
#' @param lambda2 \eqn{\lambda_2} in the framework of different rules. If not supplied, a default
#' sequence will be computed.
#' @param tau0 Parameter \eqn{\tau_0} for the \code{"SGL+SL"} penalty in the framework of different rules.
#' If not supplied, a default sequence will be computed.
#' @param alpha \eqn{\alpha} in the framework of different rules. If not supplied, a default value
#' will be used depending on \code{penalty}.
#' @param single_rule_lambda \eqn{\lambda_{single}} in the framework of single rule.
#' @param single_rule A logical value, whether the single rule framework is used. Deafult is \code{FALSE}.
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
#' @param num_lambda1 If \code{lambda1} is not specified by user, the user can still specify the length of the
#' \code{lambda1} sequence. The default length is 10.
#' @param num_lambda2 If \code{lambda2} is not specified by user, the user can still specify the length of the
#' \code{lambda2} sequence. The default length is 10.
#' @param num_tau0 If \code{tau0} is not specified by user, the user can still specify the length of the
#' \code{tau0} sequence. The default length is 11.
#' @param min_tau If \code{tau0} is not specified by user, \code{min_tau} specifies the minimum value
#' for \eqn{\tau_0}. The largest value for \eqn{\tau_0} will be \code{1 / min_tau}.
#' @param num_single_rule_lambda If \code{single_rule_lambda} is not specified, user could still specify the length
#' of the \code{single_rule_lambda} sequence. The default length is 50.
#' @import glmnet SGL Matrix
#'
#' @return An S3 object of class "mp", which contains the information of the fitted model. It could be supplied
#' to some other functions in mperosnalized package for further analysis or prediction.
#'
#' \item{penalty_parameter_sequence}{A matrix object with each row denoting a
#' configuration of the penalty parameters.}
#' \item{interceptlist}{A list object with each element denoting a vector of intercepts. The \eqn{k}th element
#' corresponds to the \eqn{k}th row in \code{penalty_parameter_sequence}.}
#' \item{betalist}{A list object with each element denoting a coefficient matrix. The \eqn{k}th element
#' corresponds to the \eqn{k}th row in \code{penalty_parameter_sequence}.}
#' \item{number_covariates}{Number of candidate covariates considered.}
#' \item{number_studies_or_outcomes}{Number of studies if \code{problem = "meta-analysis"} or number of outcomes
#' if \code{problem = "multiple outcomes"}.}
#'
#' @references Zhang, B. and Tsiatis, A. A. and Davidian, M. and Zhang, M. and Laber, E.(2012) \emph{
#' Estimating optimal treatment regimes from a classification perspective, Stat, 1(1):103-114.}
#'
#' @examples
#' set.seed(123)
#' sim_dat  = simulated_dataset(n = 200, problem = "meta-analysis")
#' Xlist = sim_dat$Xlist; Ylist = sim_dat$Ylist; Trtlist = sim_dat$Trtlist
#'
#' # fit different rules with SGL penalty for this meta-analysis problem
#' mp_mod_diff = mpersonalized(problem = "meta-analysis",
#'                             Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist,
#'                             penalty = "SGL", single_rule = FALSE)
#'
#' # fir a single rule with lasso penalty
#' mp_mod_single = mpersonalized(problem = "meta-analysis",
#'                               Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist,
#'                               penalty = "lasso", single_rule = TRUE)
#' set.seed(NULL)
#'
#' @export

mpersonalized = function(problem = c("meta-analysis", "multiple outcomes"),
                         X, Trt, P = NULL,
                         Xlist, Ylist, Trtlist, Plist = replicate(length(Xlist), NULL, simplify = FALSE),
                         typelist = replicate(length(Xlist), "continuous", simplify = FALSE),
                         penalty = c("none", "lasso", "GL", "SGL", "fused",
                                     "SGL+SL",
                                     "lasso+fused", "GL+fused", "SGL+fused"),
                         lambda1 = NULL, lambda2 = NULL, tau0 = NULL,
                         single_rule_lambda = NULL,
                         num_lambda1 = ifelse(!is.null(lambda1), length(lambda1), 10),
                         num_lambda2 = ifelse(!is.null(lambda2), length(lambda2), 10),
                         num_tau0    = ifelse(!is.null(tau0), length(tau0), 11),
                         min_tau     = 1e-2,
                         num_single_rule_lambda = ifelse(!is.null(single_rule_lambda), length(single_rule_lambda), 50),
                         alpha = NULL, single_rule = FALSE,
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
                                           single_rule = single_rule)
  modelYlist = standardized_data$modelYlist
  modelXlist = standardized_data$modelXlist

  if (single_rule == TRUE)
  {

    Ybar = standardized_data$Ybar
    Xbar = standardized_data$Xbar
    Xsd = standardized_data$Xsd

    if (penalty != "lasso" & penalty != "none")
      stop("When single rule = TRUE, the penalty must be lasso or none with default as none!")

    if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
      warning("When single rule = TRUE, the value for lambda1, lambda2, alpha are ignored!")

    if (penalty == "lasso"){

      if (is.null(single_rule_lambda)){
        lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                         penalty = penalty, single_rule = single_rule,
                                         num_single_rule_lambda = num_single_rule_lambda)

        single_rule_lambda = lambda_default$single_rule_lambda
      }

      full_model = single_rule_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                            Ybar = Ybar, Xbar = Xbar, Xsd = Xsd, lambda = single_rule_lambda)

      penalty_parameter_sequence = as.matrix(single_rule_lambda)
      colnames(penalty_parameter_sequence) = "single_rule_lambda"

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty_parameter_sequence = penalty_parameter_sequence,
                        penalty = penalty, single_rule = TRUE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        problem = problem)

    } else if (penalty == "none"){

      if (!is.null(single_rule_lambda))
        warning("When single rule = TRUE and penalty = none, the value for single_rule_lambda are ignored!")

      full_model = single_rule_linear_method(Conlist = Conlist, Xlist = Xlist)

      penalty_parameter_sequence = NULL

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty_parameter_sequence = penalty_parameter_sequence,
                        penalty = penalty, single_rule = TRUE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        problem = problem)
    }

  } else {

    Ybarlist = standardized_data$Ybarlist
    Xbarlist = standardized_data$Xbarlist
    Xsdlist = standardized_data$Xsdlist

    if (!is.null(single_rule_lambda))
        warning("When single rule = FALSE, the value for single_rule_lambda is ignored!")

    if (penalty == "none"){

      if (!is.null(lambda1) | !is.null(lambda2) | !is.null(alpha))
        warning("When penalty = none, the value for lambda1, lambda2, alpha are ignored!")

      full_model = linear_method(Conlist = Conlist, Xlist = Xlist)

      penalty_parameter_sequence = NULL

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty_parameter_sequence = penalty_parameter_sequence,
                        penalty = penalty, single_rule = FALSE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        problem = problem)

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

      if (is.null(lambda1) | is.null(lambda2))
      {
        lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                         penalty = penalty, single_rule = single_rule, alpha = alpha,
                                         num_lambda1 = num_lambda1, num_lambda2 = num_lambda2,
                                         lambda1 = lambda1, lambda2 = lambda2)

        if (is.null(lambda1))
          lambda1 = lambda_default$lambda1

        if (is.null(lambda2))
          lambda2 = lambda_default$lambda2
      }

      full_model = meta_method(modelYlist = modelYlist, modelXlist = modelXlist,
                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                               lambda1 = lambda1, lambda2 = lambda2, alpha = alpha, admm_control = admm_control)

      penalty_parameter_sequence = full_model$penalty_parameter_sequence

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty_parameter_sequence = penalty_parameter_sequence,
                        alpha = alpha, penalty = penalty, single_rule = FALSE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        problem = problem)

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

      if (is.null(lambda1)){  #  & penalty != "SGL+SL"
        lambda_default = lambda_estimate(modelXlist = modelXlist, modelYlist = modelYlist,
                                         penalty = penalty, single_rule = single_rule, alpha = alpha,
                                         num_lambda1 = num_lambda1, lambda1 = lambda1)

        lambda1 = lambda_default$lambda1
      }

      if (penalty != "SGL+SL")
      {
        full_model = sparse_group_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                               Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                               lambda = lambda1, alpha = alpha)

        penalty_parameter_sequence = as.matrix(lambda1)
        colnames(penalty_parameter_sequence) = "lambda1"
      } else
      {

        if (is.null(tau0))
        {
          tau0 <- gen_tau0(num_tau0, min_tau)
        }

        full_model = sparse_group_fused_lasso_method(modelYlist = modelYlist, modelXlist = modelXlist,
                                                     Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist,
                                                     lambda = lambda1, alpha = alpha, tau0 = tau0,
                                                     nlambda = num_lambda1)

        penalty_parameter_sequence = full_model$penalty_parameter_sequence
      }

      model_info = list(interceptlist = full_model$interceptlist, betalist = full_model$betalist,
                        penalty_parameter_sequence = penalty_parameter_sequence,
                        alpha = alpha, penalty = penalty, single_rule = FALSE,
                        number_covariates = p, number_studies_or_outcomes = q,
                        Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, Plist = Plist,
                        problem = problem)
    }
  }

  class(model_info) = "mp"
  return(model_info)

}
