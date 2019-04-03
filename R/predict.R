#' @title Prediction for a Fitted "mp" Object
#'
#' @description This function predicts optimal treatment of new subjects for a mpersonalized model.
#' If different rules are used in the fitting procedure, an overall treatment recommendation
#' based on all stuides/outcomes could be provided together with optimal treatments for
#' each study/outcome.
#'
#' @details This function predicts for each penalty parameter in the
#' \code{penalty_parameter_sequence} of the "mp" object. The overall recommended treatment is
#' given as an weighted average of the recommended treatments from each study/outcome, and the weight
#' can be specified by user.
#'
#' @param object A fitted "mp" object returned by "mpersonalized"
#' @param newx  Covariate matrix of new patients. If not supplied, by default the prediction
#' is for the original dataset in the "mp" object. Notice: when \code{problem = "meta-analysis"} and
#' the prediction is for the original dataset, subjects in each study are only predicted
#' using the treatment recommendation rule of the study they belong to.
#' @param overall_rec A logical value. If \code{overall_rec = TRUE}, an overall recommendation will
#' be provided as an weighted average of the optimal treatment from each individual study/outcome. Only useful
#' when \code{newx} is provided.
#' @param weight A weight vector for the overall recommendation, only needed when \code{overall_rec = TRUE}.
#' By default, equal weights are assigned to each study/outcome.
#' @param ... not used
#'
#' @return A list object of two elements.
#' .
#'
#' \item{opt_treatment}{A list object with each element denoting the prediction based on a penalty parameter
#' configuration in \code{mp$penalty_parameter_sequence}. Specifically, if \code{newx} is provided,
#' each element is a recommendation matrix with each row denoting a subject and each column denoting
#' a study/outcome; otherwise, each element is a list of vectors with each vector representing the optimal treatment
#' for each study/outcome. If \code{overall_rec = TRUE},
#' the weighted overall recommended treatment will be further provided as well.
#' If the overall recommened treatment is equal to 0.5, it means the weighted sum is equal for 0 and 1.}
#' \item{benefit_score}{A list object of benefit scores computed from \eqn{g_1, \dots, g_K}. Similar structure as
#' \code{opt_treatment}.}
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
#' newx = matrix(rnorm(100 * mp_mod_diff$number_covariates), nrow = 100)
#'
#' # predict on newx
#' pred_new = predict(object = mp_mod_diff, newx = newx, overall_rec = TRUE)
#'
#' # predict on old dataset
#' pred_old = predict(object = mp_mod_diff)
#' set.seed(NULL)
#' @export
predict.mp = function(object, newx = NULL, weight = NULL, overall_rec = TRUE, ...) {

  mp <- object
  q = mp$number_studies_or_outcomes
  single_rule = mp$single_rule

  if (!is.null(newx)){
    #When the newx is provided
    if (single_rule == TRUE){

      benefit_score = mapply(function(intercept, beta) intercept + newx %*% beta,
                             intercept = mp$interceptlist, beta = lapply(mp$betalist, function(x) replicate(q, x)),
                             SIMPLIFY = FALSE)
      treatment = lapply(benefit_score, function(score) ifelse(score > 0.5, 1, 0))

    } else {

      benefit_score = mapply(function(intercept, beta) sweep(newx %*% t(beta), 2, intercept, '+'),
                             intercept = mp$interceptlist, beta = mp$betalist, SIMPLIFY = FALSE)
      treatment = lapply(benefit_score, function(score) ifelse(score > 0.5, 1, 0))

    }

    #compute the overall weighted recomendation
    if (overall_rec == TRUE){

      if (is.null(weight)){
        weight = rep(1 / q, q)
      } else{
        weight = weight / sum(weight)
      }

      treatment = lapply(treatment, function(trt) {
        wt_trt = apply(sweep(trt, 2, weight, '*'), 1, sum)
        trt = cbind(trt, ifelse(wt_trt < 0.5 - 10 ^ (-6), 0, ifelse(wt_trt > 0.5 + 10 ^ (-6), 1, 0.5)))
        colnames(trt) = c(paste("Study", 1:q), "Overall Rec")
        return(trt)
      })

    } else {

      treatment = lapply(treatment, function(trt) {
        colnames(trt) = paste("Study", 1:q)
        return(trt)
      })

    }

    benefit_score = lapply(benefit_score, function(score) {
      colnames(score) = paste("Study", 1:q)
      return(score)
    })

  } else {
    #When newx is not provided
    Xlist = mp$Xlist
    benefit_score = replicate(length(mp$betalist), list())
    treatment = replicate(length(mp$betalist), list())

    if (single_rule == TRUE){

      for (i in 1:length(mp$betalist)){

        benefit_score[[i]] = mapply(function(study_intercept, study_x, study_beta) study_intercept + study_x %*% study_beta,
                                    study_x = Xlist,
                                    MoreArgs = list(study_intercept = mp$interceptlist[[i]], study_beta = mp$betalist[[i]]),
                                    SIMPLIFY = FALSE)

        treatment[[i]] = lapply(benefit_score[[i]], function(score) ifelse(score > 0.5, 1, 0))

        names(benefit_score[[i]]) = paste("Study", 1:q)
        names(treatment[[i]]) = paste("Study", 1:q)
      }

    } else {

      for (i in 1:length(mp$betalist)){

        benefit_score[[i]] = mapply(function(study_intercept, study_x, study_beta) study_intercept + study_x %*% study_beta,
                                    study_x = Xlist,
                                    study_intercept = as.list(mp$interceptlist[[i]]), study_beta = split(mp$betalist[[i]], row(mp$betalist[[i]])),
                                    SIMPLIFY = FALSE)

        treatment[[i]] = lapply( benefit_score[[i]], function(score) ifelse(score > 0.5, 1, 0))

        names(benefit_score[[i]]) = paste("Study", 1:q)
        names(treatment[[i]]) = paste("Study", 1:q)
      }

    }
  }

  pred_info = list(opt_treatment = treatment, benefit_score = benefit_score)

  return(pred_info)
}


#' @title Prediction for a Fitted "mp_cv" Object
#'
#' @description This function predicts optimal treatment of new subjects for a cross-validated mpersonalized model.
#'
#' @param object A fitted "mp_cv" object returned by "mpersonalized_cv" function
#' @param newx Covariate matrix of new patients. If not supplied, by default the prediction
#' is for the original dataset in the "mp_cv" object. Prediction results will differ
#' based on whether \code{newx} is provided or not. Similar to \code{predict.mp}.
#' @param overall_rec A logical value. If \code{overall_rec = TRUE}, an overall recommendation will
#' be provided as an weighted average of the optimal treatment from each individual study/outcome. Only useful
#' when \code{newx} is provided.
#' @param weight A weight vector for the overall recommendation, only needed when \code{overall_rec = TRUE}.
#' By default, equal weights are assigned to each study/outcome.
#' @param ... not used
#'
#' @return A list object with two elements. Similar to the returned value of \code{predict.mp}, but now it only predicts
#' for the optimal parameter penalty.
#' \item{opt_treatment}{If \code{newx} is provided, a recommendation matrix with each row denoting a subject and
#' each column denoting a study/outcome; otherwise, each element is a list of vectors with each vector representing the optimal treatment
#' for each study/outcome. If \code{overall_rec = TRUE},
#' the weighted overall recommended treatment will be further provided as well.
#' If the overall recommened treatment is equal to 0.5, it means the weighted sum is equal for 0 and 1.}
#' \item{benefit_score}{Benefit scores computed from \eqn{g_1, \dots, g_K}. Similar to structure of \code{opt_treatment}.}
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
#' newx = matrix(rnorm(100 * mp_cvmod_diff$number_covariates), nrow = 100)
#'
#' # predict on newx
#' pred_new = predict(object = mp_cvmod_diff, newx = newx, overall_rec = TRUE)
#'
#' # predict on old dataset
#' pred_old = predict(object = mp_cvmod_diff)
#' set.seed(NULL)
#' @export
predict.mp_cv = function(object, newx = NULL, weight = NULL, overall_rec = TRUE, ...) {

  mp_cv <- object
  q = mp_cv$number_studies_or_outcomes
  single_rule = mp_cv$single_rule

  if (!is.null(newx)){
    #when newx is provided
    if (single_rule == TRUE){

      benefit_score = mp_cv$intercept + newx %*% replicate(q, mp_cv$beta)
      treatment = ifelse(benefit_score > 0.5, 1, 0)

    } else {

      benefit_score = sweep(newx %*% t(mp_cv$beta), 2, mp_cv$intercept, '+')
      treatment = ifelse(benefit_score > 0.5, 1, 0)

    }

    if (overall_rec == TRUE){

      if (is.null(weight)){
        weight = rep(1 / q, q)
      } else{
        weight = weight / sum(weight)
      }

      #compute the overall weighted recomendation
      wt_trt = apply(sweep(treatment, 2, weight, '*'), 1, sum)
      treatment = cbind(treatment, ifelse(wt_trt < 0.5 - 10 ^ (-6), 0, ifelse(wt_trt > 0.5 + 10 ^ (-6), 1, 0.5)))
      colnames(treatment) = c(paste("Study", 1:q), "Overall Rec")

    } else {

      colnames(treatment) = paste("Study", 1:q)

    }

    colnames(benefit_score) = paste("Study", 1:q)

  } else {
    #When newx is not provided
    Xlist = mp_cv$Xlist
    benefit_score = replicate(q, list())
    treatment = replicate(q, list())

    if (single_rule == TRUE){

      benefit_score = mapply(function(study_intercept, study_x, study_beta) study_intercept + study_x %*% study_beta,
                             study_x = Xlist,
                             MoreArgs = list(study_intercept = mp_cv$intercept, study_beta = mp_cv$beta),
                             SIMPLIFY = FALSE)

      treatment = lapply(benefit_score, function(score) ifelse(score > 0.5, 1, 0))

      names(benefit_score) = paste("Study", 1:q)
      names(treatment) = paste("Study", 1:q)
    } else {

      benefit_score = mapply(function(study_intercept, study_x, study_beta) study_intercept + study_x %*% study_beta,
                                  study_x = Xlist,
                                  study_intercept = as.list(mp_cv$intercept), study_beta = split(mp_cv$beta, row(mp_cv$beta)),
                                  SIMPLIFY = FALSE)

      treatment = lapply(benefit_score, function(score) ifelse(score > 0.5, 1, 0))

      names(benefit_score) = paste("Study", 1:q)
      names(treatment) = paste("Study", 1:q)
    }
  }

  pred_info = list(opt_treatment = treatment, benefit_score = benefit_score)

  return(pred_info)
}
