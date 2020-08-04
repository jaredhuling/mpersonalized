contrast_standardize = function(Conlist, Xlist, single_rule, standardize = TRUE)
{
  #surrogate <- match.arg(surrogate)
  if (single_rule == FALSE)
  {
    sConlist = lapply(Conlist, function(y) as.numeric(y > 0))
    Wlist = lapply(Conlist, abs)

    dataWlist = lapply(Wlist, sum)
    adj_Wlist = mapply(function(w, dataw) w / dataw, w = Wlist,
                       dataw = dataWlist, SIMPLIFY = FALSE)

    Xsdlist = lapply(Xlist, function(x) apply(x, 2, sd))
    std_Xlist = mapply(function(x, xsd) sweep(x, 2, xsd, '/'),
                       x = Xlist, xsd = Xsdlist, SIMPLIFY = FALSE)
    Ybarlist = mapply(function(y, w) sum(y * w),
                      y = sConlist, w = adj_Wlist, SIMPLIFY = FALSE)
    Xbarlist = mapply(function(x, w) apply(sweep(x, 1, w, '*'), 2, sum),
                      x = std_Xlist, w = adj_Wlist, SIMPLIFY = FALSE)

    std_Ylist = mapply(function(y, mean) y - mean,
                       y = sConlist, mean = Ybarlist, SIMPLIFY = FALSE)


    std_Xlist = mapply(function(x, mean) sweep(x, 2, mean, '-'),
                       x = std_Xlist, mean = Xbarlist, SIMPLIFY = FALSE)

    if (standardize)
    {
      modelYlist = mapply(function(y, w) y * sqrt(w),
                          y = std_Ylist, w = adj_Wlist, SIMPLIFY = FALSE)
      modelXlist = mapply(function(x, w) sweep(x, 1, sqrt(w), '*'),
                          x = std_Xlist, w = adj_Wlist, SIMPLIFY = FALSE)
      Wlist <- lapply(Wlist, function(lst) rep(1, length(lst)))
    } else
    {
      modelYlist <- sConlist
      modelXlist <- Xlist

      Xbarlist <- lapply(Xbarlist, function(lst) rep(0, length(lst)))
      Ybarlist <- lapply(Ybarlist, function(lst) rep(0, length(lst)))
      Xsdlist  <- lapply(Xsdlist, function(lst)  rep(1, length(lst)))
    }



    return(list(modelYlist = modelYlist, modelXlist = modelXlist,
                Wlist = Wlist,
                Ybarlist = Ybarlist, Xbarlist = Xbarlist, Xsdlist = Xsdlist))

  } else
  {
    sConlist = lapply(Conlist, function(y) as.numeric(y > 0))
    Wlist = lapply(Conlist, function(y) abs(y))

    dataWlist = lapply(Wlist, sum)
    adj_Wlist = mapply(function(w, dataw) w / dataw, w = Wlist, dataw = dataWlist,
                       SIMPLIFY = FALSE)

    poolX = do.call(rbind, Xlist)
    poolY = unlist(sConlist)
    adj_W = unlist(adj_Wlist)
    Xsd = apply(poolX, 2, sd)

    std_Xlist = lapply(Xlist, function(x) sweep(x, 2, Xsd, '/'))

    Ybar = sum(adj_W * poolY) / sum(adj_W)
    Xbar = apply(sweep(poolX, 1, adj_W, '*'), 2, sum) / sum(adj_W)

    std_Ylist = lapply(sConlist, function(y) y - Ybar)
    std_Xlist = lapply(std_Xlist, function(x) sweep(x, 2, Xbar, '-'))

    if (standardize)
    {
      modelYlist = mapply(function(y, w) y * sqrt(w),
                          y = std_Ylist, w = adj_Wlist, SIMPLIFY = FALSE)
      modelXlist = mapply(function(x, w) sweep(x, 1, sqrt(w), '*'),
                          x = std_Xlist, w = adj_Wlist, SIMPLIFY = FALSE)

      Wlist <- lapply(Wlist, function(lst) rep(1, length(lst)))
    } else
    {
      modelYlist <- sConlist
      modelXlist <- Xlist

      Xsd  <- rep(1, length(Xsd))
      Xbar <- rep(0, length(Xbar))
      Ybar <- rep(0, length(Ybar))
    }

    return(list(modelYlist = modelYlist, modelXlist = modelXlist,
                Wlist = Wlist,
                Ybar = Ybar, Xbar = Xbar, Xsd = Xsd))
  }
}
