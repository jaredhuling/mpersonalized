#' @title Coefficients of a Fitted "mp" Object
#'
#' @description This function provides \code{coef} method for "mp" class objects.
#'
#' @param object A fitted "mp" object returned by "mpersonalized".
#' @param ... not used
#'
#' @return A list object. Each element in the list is the fitted coefficients
#' corresponding to one penalty parameter value in \code{mp$penalty_parameter_sequence}.
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
#' mp_coef = coef(mp_mod_diff)
#' set.seed(NULL)
#' @export
#' @importFrom graphics plot
#' @importFrom methods as
#' @importFrom stats binomial gaussian glm lm predict quantile rbinom rnorm sd
coef.mp = function(object, ...) {

  mp <- object
  single_rule = mp$single_rule

  if (single_rule == TRUE){

    coeflist = mapply(function(intercept, beta) c(intercept, beta),
                      intercept = mp$interceptlist, beta = mp$betalist,
                      SIMPLIFY = FALSE)
  } else {

    coeflist = mapply(function(intercept, beta) cbind(intercept, beta),
                      intercept = mp$interceptlist, beta = mp$betalist,
                      SIMPLIFY = FALSE)
  }

  return(coeflist)
}

#' @title Coefficients of a Fitted "mp_cv" Object
#'
#' @description This function provides \code{coef} method for "mp_cv" class objects.
#'
#' @param object A fitted "mp" object returned by "mpersonalized".
#' @param ... not used
#'
#' @return The fitted coefficients corresponding to the optimal penalty parameter.
#'
#' @examples
#' set.seed(123)
#' sim_dat  = simulated_dataset(n = 200, problem = "meta-analysis")
#' Xlist = sim_dat$Xlist; Ylist = sim_dat$Ylist; Trtlist = sim_dat$Trtlist
#'
#' # fit different rules with lasso penalty for this meta-analysis problem
#' mp_cvmod_diff = mpersonalized_cv(problem = "meta-analysis",
#'                                  Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist,
#'                                  penalty = "lasso", single_rule = FALSE)
#'
#' mp_coef = coef(mp_cvmod_diff)
#' set.seed(NULL)
#' @export
coef.mp_cv = function(object, ...) {

  mp_cv <- object
  single_rule = mp_cv$single_rule

  if (single_rule == TRUE){

    coefs = c(mp_cv$intercept, mp_cv$beta)
  } else {

    coefs = cbind(mp_cv$intercept, mp_cv$beta)
  }

  return(coefs)
}

