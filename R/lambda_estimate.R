#' @import glmnet SGL Matrix genlasso

lambda_estimate =  function(modelXlist, modelYlist, penalty, single_rule, alpha,
                            num_lambda1, num_lambda2, num_single_rule_lambda){
  lambda_estimate = NULL

  #need to estimate lambda1?
  if (single_rule == FALSE &
      penalty != "fused"){

    q = length(modelXlist)
    p = dim(modelXlist[[1]])[2]

    x = as.matrix(bdiag(modelXlist))
    y = unlist(modelYlist)
    total_n = length(y)

    x = sqrt(total_n) * x; y = sqrt(total_n) * y

    data = list(x = x, y = y)
    SGL_lambda = SGL(data = data, index = rep(1 : p, q),
                     alpha = alpha, standardize = FALSE, nlam = num_lambda1)$lambdas
    lambda1 = SGL_lambda

    lambda_estimate$lambda1 = lambda1
  }

  #need to estimate lambda2?
  if (single_rule == FALSE &
      penalty %in% c("fused", "lasso+fused", "GL+fused", "SGL+fused")){

    q = length(modelXlist)
    p = dim(modelXlist[[1]])[2]

    x = as.matrix(bdiag(modelXlist))
    y = unlist(modelYlist)
    D = NULL
    for (i in 2:q)
      for (j in 1:(i-1)){
        newrow = numeric(q)
        newrow[i] = -1
        newrow[j] = 1
        D = rbind(D, newrow)
      }
    D = as.matrix(bdiag(replicate(p, D, simplify = FALSE)))

    genlasso_lambda = genlasso(y = y, X = x, D = D)$lambda
    lambda2 = quantile(genlasso_lambda, probs = seq(1, 0, -1 / (num_lambda2 - 1)), names = FALSE)

    lambda_estimate$lambda2 = lambda2
  }


  #need to estimate single_rule_lambda?
  if (single_rule == TRUE){
    x = do.call(rbind, modelXlist)
    y = unlist(modelYlist)
    total_n = length(y)

    x = sqrt(total_n) * x; y = sqrt(total_n) * y

    lasso_lambda = glmnet(x = x, y = y, family = "gaussian",
                         standardize = FALSE, intercept  = FALSE,
                         nlambda = num_single_rule_lambda)$lambda
    single_rule_lambda = lasso_lambda

    lambda_estimate$single_rule_lambda = single_rule_lambda
  }

  return(lambda_estimate)
}
