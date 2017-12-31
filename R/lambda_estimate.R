#' @title Default Penalty Parameter Sequence if Not Given
#' @param modelXlist the Xlist in contrast framework after standardization
#' @param modelYlist the Ylist in contrast framework after standardization
#' @param penalty penalty type
#' @param unique_rule whether a unique treatment rule is required
#' @import glmnet SGL Matrix genlasso
#'
#' @return estimated lambda for required penalty if not provided
#' @export
lambda_estimate =  function(modelXlist, modelYlist, penalty, unique_rule, alpha){
  lambda_estimate = NULL

  #need to estimate lambda1?
  if (unique_rule == FALSE &
      penalty != "fused"){

    q = length(modelXlist)
    p = dim(modelXlist[[1]])[2]

    x = as.matrix(bdiag(modelXlist))
    y = unlist(modelYlist)
    total_n = length(y)

    x = sqrt(total_n) * x; y = sqrt(total_n) * y

    data = list(x = x, y = y)
    SGL_lambda = SGL(data = data, index = rep(1 : p, q),
                     alpha = alpha, standardize = FALSE)$lambdas
    lambda1 = SGL_lambda

    lambda_estimate$lambda1 = lambda1
  }

  #need to estimate lambda2?
  if (unique_rule == FALSE &
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
    lambda2 = quantile(genlasso_lambda, probs = seq(1, 0.1, -0.1))

    lambda_estimate$lambda2 = lambda2
  }


  #need to estimate unique_rule_lambda?
  if (unique_rule == TRUE){
    x = do.call(rbind, modelXlist)
    y = unlist(modelYlist)
    total_n = length(y)

    x = sqrt(total_n) * x; y = sqrt(total_n) * y

    lasso_lambda = glmnet(x = x, y = y, family = "gaussian",
                         standardize = FALSE, intercept  = FALSE)$lambda
    unique_rule_lambda = lasso_lambda

    lambda_estimate$unique_rule_lambda = unique_rule_lambda
  }

  return(lambda_estimate)
}
