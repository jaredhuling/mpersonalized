#' @title Plot for a 'mp' Class Object
#'
#' @param mp the 'mp' class object returned by MetaPersonalzied function
#' @param ind1 the index of the wanted lambda1
#' @param ind2 the index of the wanted lambda2
#' @param unique_ind the index of the wanted unique_rule_lambda
#'
#' @import ggplot2 gridExtra
#' @return a list object with \eqn{k}th element denoting the plot of study k
#' @export
plot.mp = function(mp, ind1, ind2, unique_ind){

  Ylist = mp$Ylist
  Trtlist = mp$Trtlist
  unique_rule = mp$unique_rule
  q = mp$number_studies
  penalty = mp$penalty

  if (unique_rule == TRUE){

    pred = predict(mp)$treatment[[unique_ind]]

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


#' @title Plot for a 'mp_cv' Class Object
#'
#' @param mp_cv the 'mp_cv' class object returned by MetaPersonalzied function
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
