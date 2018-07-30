meta_method = function(modelYlist, modelXlist, Ybarlist, Xbarlist, Xsdlist, lambda1, lambda2, alpha, admm_control){

  q = length(modelXlist)
  p = dim(modelXlist[[1]])[2]

  tempx = bdiag(modelXlist)
  x = matrix(0, nrow = nrow(tempx),ncol = ncol(tempx))
  for (i in 1:p)
    for (j in 1:q)
      x[,((i-1)*q+j)] = tempx[,((j-1)*p+i)]

  y = unlist(modelYlist)

  ngencon = q * (q + 1) / 2 # number of generalized lasso constraint

  nlambda1 = length(lambda1)
  nlambda2 = length(lambda2)

  betalist = vector("list", nlambda1 * nlambda2)
  interceptlist = vector("list", nlambda1 * nlambda2)
  iterslist = vector("list", nlambda1 * nlambda2)


  admm_lambda1 = ifelse(lambda1 != 0, lambda1 * (1 - alpha) * sqrt(q), 0)
  admm_lambda2 = ifelse(lambda1 != 0, lambda1 * alpha, 0)
  admm_lambda3 = lambda2

  penalty_parameter_sequence = matrix(0, ncol = 2, nrow = nlambda1 * nlambda2)
  colnames(penalty_parameter_sequence) = c("lambda1", "lambda2")

  for (ind1 in 1:nlambda1)
    for (ind2 in 1:nlambda2){
      pen1 = admm_lambda1[ind1]
      pen2 = admm_lambda2[ind1]
      pen3 = admm_lambda3[ind2]

      penalty_parameter_sequence[(ind1 - 1) * nlambda2 + ind2,] = c(lambda1[ind1], lambda2[ind2])

      result = do.call(admm_optim, c(admm_control, list(x = x, y = y, p = p, q = q,
                                                        lambda1 = pen1,lambda2 = pen2, lambda3 = pen3)))

      beta = matrix(result$beta, nrow = q, ncol = p, byrow = FALSE)
      intercept = unlist(Ybarlist) - apply(beta * matrix(unlist(Xbarlist), nrow = q, ncol = p, byrow = TRUE), 1, sum)
      beta = beta / matrix(unlist(Xsdlist), nrow = q, ncol = p, byrow = TRUE)

      betalist[[(ind1 - 1) * nlambda2 + ind2]]      = beta
      interceptlist[[(ind1 - 1) * nlambda2 + ind2]] = intercept
      iterslist[[(ind1 - 1) * nlambda2 + ind2]]     = result$iters
    }

  return(list(interceptlist = interceptlist,
              betalist = betalist,
              iterslist = iterslist,
              penalty_parameter_sequence = penalty_parameter_sequence))
}
