#' @export

single_rule_linear_method = function(Conlist, Xlist){

  sConlist = lapply(Conlist, function(y) as.numeric(y > 0))
  Wlist = lapply(Conlist, function(y) abs(y))

  dataWlist = lapply(Wlist, sum)
  adj_Wlist = mapply(function(w, dataw) w / dataw, w = Wlist,
                          dataw = dataWlist, SIMPLIFY = FALSE)

  x = do.call(rbind, Xlist)
  y = unlist(sConlist)
  w = unlist(adj_Wlist)

  linear_mod = lm(y ~ x, weights = w)

  return(list(interceptlist = list(linear_mod$coef[1]), betalist = list(linear_mod$coef[-1])))
}
