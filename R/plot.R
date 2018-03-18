#' @title Plot for a 'mp' Class Object.
#'
#' @details This function plots the results for estimated treatment effects. Depending on the received
#' treatment and recommended treatment, the group means of the outcome are computed and the relations between them are plotted.
#' This plot provides a sanity check of the treatment recommendation rule. By specifiying the index of the penalty parameters, we can
#' obtain the plots of the corresponding treatment recommendation rule.
#'
#' @param mp the 'mp' class object returned by \code{mpersonalzied} function
#' @param ind1 the index of the lambda1 if different rules are used
#' @param ind2 the index of the lambda2 if different rules are used
#' @param single_ind the index of the single_rule_lambda if an single rule is used
#'
#' @import ggplot2 gridExtra
#' @return a list object with \eqn{k}th element denoting the plot of study k
#' @export
plot.mp = function(mp, ind1, ind2, single_ind){

  Ylist = mp$Ylist
  Trtlist = mp$Trtlist
  single_rule = mp$single_rule
  q = mp$number_studies
  penalty = mp$penalty

  if (single_rule == TRUE){

    pred = predict(mp)$treatment[[single_ind]]

  } else {

    if (penalty == "linear"){

      pred = predict(mp)$treatment[[1]]

    } else {

      nlambda1 = length(mp$lambda1)
      nlambda2 = length(mp$lambda2)

      pred = predict(mp)$treatment[[(ind1 - 1) * nlambda2 + ind2]]

    }
  }

  #group 1 defined as receive 1 and recommend 1; group 2 as receive 0 and recommend 1
  #group 3 defined as receive 0 and recommend 0; group 4 as receive 1 and recommend 0
  meanlist = mapply(function(Trt, pred, Y) c(mean(Y[Trt == 1 & pred == 1]), mean(Y[Trt == 0 & pred == 1]),
                                             mean(Y[Trt == 0 & pred == 0]), mean(Y[Trt == 1 & pred == 0])),
                    Trt = Trtlist, pred = pred, Y = Ylist, SIMPLIFY = FALSE)

  plotlist = replicate(q, list())
  for (i in 1:q){
    plot_dat = data.frame(mean = meanlist[[i]], recommend = as.factor(c(1, 1, 0, 0)),
                          received = as.factor(c(1, 0, 0, 1)))
    plotlist[[i]] = ggplot(data = plot_dat, aes(y = mean, x = received, group = recommend)) +
      geom_line(aes(color = recommend)) + ggtitle(label = paste("Study ",i))
  }

  return(plotlist)
}


#' @title Plot for a 'mp' Class Object.
#'
#' @details This function plots the results for estimated treatment effects by using the estimated optimal
#' treatment recommendation rule obtained from corss validation.
#'
#' @param mp_cv the 'mp_cv' class object returned by \code{mpersonalzied_cv} function
#'
#' @import ggplot2 gridExtra
#' @return a list object with \eqn{k}th element denoting the plot of study k
#' @export
plot.mp_cv = function(mp_cv){

  Ylist = mp_cv$Ylist
  Trtlist = mp_cv$Trtlist
  q = mp_cv$number_studies

  pred = predict(mp_cv)$treatment

  #group 1 defined as receive 1 and recommend 1; group 2 as receive 0 and recommend 1
  #group 3 defined as receive 0 and recommend 0; group 4 as receive 1 and recommend 0
  meanlist = mapply(function(Trt, pred, Y) c(mean(Y[Trt == 1 & pred == 1]), mean(Y[Trt == 0 & pred == 1]),
                                             mean(Y[Trt == 0 & pred == 0]), mean(Y[Trt == 1 & pred == 0])),
                    Trt = Trtlist, pred = pred, Y = Ylist, SIMPLIFY = FALSE)

  plotlist = replicate(q, list())
  for (i in 1:q){
    plot_dat = data.frame(mean = meanlist[[i]], recommend = as.factor(c(1, 1, 0, 0)),
                          received = as.factor(c(1, 0, 0, 1)))
    plotlist[[i]] = ggplot(data = plot_dat, aes(y = mean, x = received, group = recommend)) +
      geom_line(aes(color = recommend)) + ggtitle(label = paste("Study ",i))
  }

  return(plotlist)
}
