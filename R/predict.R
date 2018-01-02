#' @title Predict for "mp" object
#'
#' @description This function predict the benefit scores and optimal treatment for new patients
#'
#' @param mp the fitted "mp" object returned by "mpersonalized"
#' @param newx the covariate matrix of the new patients. If newx = NULL, then the prediction is made for each data set based on the rule of that study/outcome.
#' @param overall_rec a logical value. If \code{TRUE}, an overall recommendation will be made weighted by the "weight" parameter.
#' @param weight a weight vector for the overall recommendation. If leave as \code{NULL}, a equally weighted recommendation will be made.
#'
#' @return a list of results with each element in the list corresponding to the prediciton by each different penalty parameter. Each element in this list contains:
#' \item{treatment}{recommended treatment for each patient for each study/outcome. If \code{overall_rec = TRUE}, the weighted overall recommended treatment will be computed as well.
#' If the overall recommened treatment is equal to 0.5, it means the sum of weight is equal for 0 and 1.}
#' \item{benefit_score}{the benefit score computed from \eqn{g_1, \dots, g_K}}
#' @export
predict.mp = function(mp, newx = NULL, weight = NULL, overall_rec = TRUE){

  q = mp$number_studies
  unique_rule = mp$unique_rule

  if (!is.null(newx)){
    #When the newx is provided
    if (unique_rule == TRUE){

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

    if (unique_rule == TRUE){

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

  pred_info = list(treatment = treatment, benefit_score = benefit_score)

  return(pred_info)
}


#' @title Predict for "mp_cv" object
#'
#' @description This function predict the benefit scores and optimal treatment for new patients
#'
#' @param mp_cv the fitted "mp_cv" object returned from "mpersonalized_cv" function
#' @param newx the covariate matrix of the new patients. If newx = NULL, then the prediction is made for each data set based on the rule of that study/outcome.
#' @param overall_rec a logical value. If \code{TRUE}, an overall recommendation will be made weighted by the "weight" parameter.
#' @param weight a weight vector for the overall recommendation. If leave as \code{NULL}, a equally weighted recommendation will be made.
#'
#' @return the prediciton by using the optimal penalty parameter selected by cross validation. It contains:
#' \item{treatment}{recommended treatment for each patient for each study/outcome. If \code{overall_rec = TRUE}, the weighted overall recommended treatment will be computed as well.
#' If the overall recommened treatment is equal to 0.5, it means the sum of weight is equal for 0 and 1.}
#' \item{benefit_score}{the benefit score computed from \eqn{g_1, \dots, g_K}}
#' @export
predict.mp_cv = function(mp_cv, newx = NULL, weight = NULL, overall_rec = TRUE){

  q = mp_cv$number_studies
  unique_rule = mp_cv$unique_rule

  if (!is.null(newx)){
    #when newx is provided
    if (unique_rule == TRUE){

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

    if (unique_rule == TRUE){

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

  pred_info = list(treatment = treatment, benefit_score = benefit_score)

  return(pred_info)
}
