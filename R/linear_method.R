linear_method = function(Conlist, Xlist){

  sConlist = lapply(Conlist, function(y) as.numeric(y > 0))
  Wlist = lapply(Conlist, function(y) abs(y))

  coeflist = mapply(function(x, y, w) lm(y ~ x, weights = w)$coef,
                    y = sConlist, x = Xlist, w = Wlist, SIMPLIFY = FALSE)

  intercept = mapply(function(coef) coef[1], coef = coeflist, SIMPLIFY = TRUE)
  beta = do.call(rbind, mapply(function(coef) coef[-1], coef = coeflist, SIMPLIFY = FALSE))

  return(list(interceptlist = list(intercept), betalist = list(beta)))
}
