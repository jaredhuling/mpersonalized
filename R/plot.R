#' @title Interaction Plot for an "mp" Class Object.
#'
#' @description This function plots interaction between received treatment and recommended treatment,
#' which provides an estimate of treatment effect of the identified subgroup.
#'
#' @details In the interaction plot, each point is the group mean given a received treatment
#' and a recommended treatment. Although usually
#' overestimating treatment effect in training set, interaction plots provides a sanity check for treatment
#' recommendation rules. Given a specific index of penalty parameter, the function
#' plots corresponding interaction plots.
#'
#' @param mp A fitted "mp" class object returned by \code{mpersonalzied} function
#' @param penalty_index The index of penalty parameter configuration in \code{mp$penalty_parameter_sequence}.
#' When \code{mp$penalty = "none"}, \code{penalty_index} is automatically set to be 1.
#'
#' @import ggplot2 gridExtra
#' @return A list object with each element as the interaction plots for a penalty parameter configuration.
#'
#' @examples
#' set.seed(123)
#' sim_dat  = simulated_dataset(n = 200, problem = "meta-analysis")
#' Xlist = sim_dat$Xlist; Ylist = sim_dat$Ylist; Trtlist = sim_dat$Trtlist
#'
#' # fit different rules with SGL penalty for this meta-analysis problem
#' mp_mod_diff = mpersonalized(problem = "meta-analysis",
#'                             Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist,
#'                             penalty = "lasso", single_rule = FALSE)
#'
#' # interaction plot of the 5th penalty parameter
#' plots = plot(mp = mp_mod_diff, penalty_index = 5)
#' set.seed(NULL)
#' @export
plot.mp = function(mp, penalty_index){

  Ylist = mp$Ylist
  Trtlist = mp$Trtlist
  Plist = mp$Plist
  single_rule = mp$single_rule
  q = mp$number_studies_or_outcomes
  penalty = mp$penalty
  problem = mp$problem


  if (penalty == "none"){

    pred = predict(mp)$opt_treatment[[1]]

  } else {

    if (missing(penalty_index))
      stop("For penalty not equal to 'none', penalty_index must be inputed!")

    pred = predict(mp)$opt_treatment[[penalty_index]]

  }


  #group 1 defined as receive 1 and recommend 1; group 2 as receive 0 and recommend 1
  #group 3 defined as receive 0 and recommend 0; group 4 as receive 1 and recommend 0
  meanlist = mapply(function(Trt, pred, Y,  P) c(sum(Y * as.numeric(Trt == 1 & pred == 1) / P) / sum(as.numeric(Trt == 1 & pred == 1) / P),
                                                 sum(Y * as.numeric(Trt == 0 & pred == 1) / (1 - P)) / sum(as.numeric(Trt == 0 & pred == 1) / (1 - P)),
                                                 sum(Y * as.numeric(Trt == 0 & pred == 0) / (1 - P)) / sum(as.numeric(Trt == 0 & pred == 0) / (1 - P)),
                                                 sum(Y * as.numeric(Trt == 1 & pred == 0) / P) / sum(as.numeric(Trt == 1 & pred == 1) / P)),
                    Trt = Trtlist, pred = pred, Y = Ylist, P = Plist, SIMPLIFY = FALSE)

  plotlist = replicate(q, list())
  for (i in 1:q){
    plot_dat = data.frame(mean = meanlist[[i]], recommend = as.factor(c(1, 1, 0, 0)),
                          received = as.factor(c(1, 0, 0, 1)))
    plotlist[[i]] = ggplot(data = plot_dat, aes(y = mean, x = recommend, group = received)) +
      geom_line(aes(color = received), size = 2, alpha = 0.4) + geom_point(size = 3, aes(color = received))

    if (problem == "meta-analysis")
      plotlist[[i]] = plotlist[[i]] + ggtitle(label = paste("Study ",i))

    if (problem == "multiple outcomes")
      plotlist[[i]] = plotlist[[i]] + ggtitle(label = paste("Outcome ",i))
  }

  return(plotlist)
}


#' @title Interaction Plot for an "mp_cv" Class Object.
#'
#' @description This function plots interaction between received treatment and recommended treatment,
#' given the optimal penalty parameter.
#'
#' @param mp_cv A fitted 'mp_cv' class object returned by \code{mpersonalzied_cv} function
#'
#' @import ggplot2 gridExtra
#' @return A list object representing the interaction plots for the optimal penalty parameter configuration.
#' Specifically, \eqn{k}th element is the interaction plot for the \eqn{k}th study/outcome.
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
#' plots = plot(mp_cv = mp_cvmod_diff)
#' set.seed(NULL)
#' @export
plot.mp_cv = function(mp_cv){

  Ylist = mp_cv$Ylist
  Trtlist = mp_cv$Trtlist
  Plist = mp_cv$Plist
  q = mp_cv$number_studies_or_outcomes
  problem = mp_cv$problem

  pred = predict(mp_cv)$opt_treatment

  #group 1 defined as receive 1 and recommend 1; group 2 as receive 0 and recommend 1
  #group 3 defined as receive 0 and recommend 0; group 4 as receive 1 and recommend 0
  meanlist = mapply(function(Trt, pred, Y,  P) c(sum(Y * as.numeric(Trt == 1 & pred == 1) / P) / sum(as.numeric(Trt == 1 & pred == 1) / P),
                                                 sum(Y * as.numeric(Trt == 0 & pred == 1) / (1 - P)) / sum(as.numeric(Trt == 0 & pred == 1) / (1 - P)),
                                                 sum(Y * as.numeric(Trt == 0 & pred == 0) / (1 - P)) / sum(as.numeric(Trt == 0 & pred == 0) / (1 - P)),
                                                 sum(Y * as.numeric(Trt == 1 & pred == 0) / P) / sum(as.numeric(Trt == 1 & pred == 1) / P)),
                    Trt = Trtlist, pred = pred, Y = Ylist, P = Plist, SIMPLIFY = FALSE)

  plotlist = replicate(q, list())
  for (i in 1:q){
    plot_dat = data.frame(mean = meanlist[[i]], recommend = as.factor(c(1, 1, 0, 0)),
                          received = as.factor(c(1, 0, 0, 1)))
    plotlist[[i]] = ggplot(data = plot_dat, aes(y = mean, x = recommend, group = received)) +
      geom_line(aes(color = received), size = 2, alpha = 0.4) + geom_point(size = 3, aes(color = received))

    if (problem == "meta-analysis")
      plotlist[[i]] = plotlist[[i]] + ggtitle(label = paste("Study ",i))

    if (problem == "multiple outcomes")
      plotlist[[i]] = plotlist[[i]] + ggtitle(label = paste("Outcome ",i))
  }

  return(plotlist)
}
