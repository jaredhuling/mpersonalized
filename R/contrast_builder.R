#' @export

contrast_builder = function(X, Y, ori_Trt, P, response_model = c("lasso", "linear"),
                            type = c("continuous", "binary"), contrast_builder_folds = 3){

  type = match.arg(type)
  response_model = match.arg(response_model)

  Trt = 2 * (ori_Trt - 0.5)
  n = dim(X)[1]; p = dim(X)[2]
  IntX = sweep(X, 1, Trt, '*')
  CbX=cbind(X, Trt, IntX)

  CbX0 = cbind(X, rep(-1, n), -X)
  CbX1 = cbind(X, rep(1, n), X)

  if (response_model == "lasso"){
    if (type == "continuous"){
      lasmod = cv.glmnet(y = Y, x = CbX, family = "gaussian",
                         nfolds = contrast_builder_folds, nlambda = 10)
      Trteff0 = predict(lasmod, newx = CbX0, s = "lambda.min")
      Trteff1 = predict(lasmod, newx = CbX1, s = "lambda.min")
    } else if (type == "binary") {
      lasmod = cv.glmnet(y = Y, x = CbX, family = "binomial",
                         nfolds = contrast_builder_folds, nlambda = 10)
      Trteff0 = predict(lasmod, newx = CbX0, s = "lambda.min")
      Trteff0 = exp(Trteff0) / (1 + exp(Trteff0))
      Trteff1 = predict(lasmod, newx = CbX1, s = "lambda.min")
      Trteff1 = exp(Trteff1) / (1 + exp(Trteff1))
    }
  } else if (response_model == "linear"){
    dat = data.frame(y = y, x = CbX)
    if (type == "continuous")
      glmmod = glm(y ~ ., data = dat, family = gaussian())
    if (type == "binary")
      glmmod = glm(y ~ ., data = dat, family = binomial())

    newdat0 = data.frame(x = CbX0)
    Trteff0 = predict(glmmod, newdat0, type = "response")
    newdat1 = data.frame(x = CbX1)
    Trteff1 = predict(glmmod, newdat1, type = "response")
  }

  Y_adj = Y - (1 - P) * Trteff1 - P * Trteff0

  est_cont = ori_Trt * Y_adj / P - (1 - ori_Trt) * Y_adj / (1 - P)

  return(est_cont)

}
