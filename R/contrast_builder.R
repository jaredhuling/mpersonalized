contrast_builder = function(X, Y, ori_Trt, P, eff_aug = TRUE,
                            response_model = c("lasso", "linear", "randomForest"),
                            type = c("continuous", "binary"),
                            contrast_builder_folds = 10)
{


  type = match.arg(type)
  response_model = match.arg(response_model)

  Trt = 2 * (ori_Trt - 0.5)
  n = dim(X)[1]; p = dim(X)[2]

  if (eff_aug == TRUE)
  {
    IntX = sweep(X, 1, Trt, '*')
    CbX=cbind(X, Trt, IntX)

    CbX0 = cbind(X, rep(-1, n), -X)
    CbX1 = cbind(X, rep(1, n), X)

    if (response_model == "lasso")
    {
      if (type == "continuous")
      {
        lasmod = cv.glmnet(y = Y, x = CbX, family = "gaussian",
                           nfolds = contrast_builder_folds, nlambda = 10)
        Trteff0 = predict(lasmod, newx = CbX0, s = "lambda.min")
        Trteff1 = predict(lasmod, newx = CbX1, s = "lambda.min")
      } else if (type == "binary")
      {
        lasmod = cv.glmnet(y = Y, x = CbX, family = "binomial",
                           nfolds = contrast_builder_folds, nlambda = 10)
        Trteff0 = predict(lasmod, newx = CbX0, s = "lambda.min")
        Trteff0 = exp(Trteff0) / (1 + exp(Trteff0))
        Trteff1 = predict(lasmod, newx = CbX1, s = "lambda.min")
        Trteff1 = exp(Trteff1) / (1 + exp(Trteff1))
      }
    } else if (response_model == "linear")
    {
      dat = data.frame(y = Y, x = CbX)
      if (type == "continuous")
        glmmod = glm(y ~ ., data = dat, family = gaussian())
      if (type == "binary")
        glmmod = glm(y ~ ., data = dat, family = binomial())

      newdat0 = data.frame(x = CbX0)
      Trteff0 = predict(glmmod, newdat0, type = "response")
      newdat1 = data.frame(x = CbX1)
      Trteff1 = predict(glmmod, newdat1, type = "response")
    } else if (response_model == "randomForest")
    {
      dat = data.frame(y = Y, x = X, Trt = as.character(Trt))
      if (type == "continuous")
      {
        rfmod = randomForest(y ~ ., data = dat, family = gaussian())
        prd_type <- "response"

        Trt_levels <- sort(unique(Trt))

        CbX0_rf <- data.frame(x = X, Trt = as.character(Trt_levels[1]))
        CbX1_rf <- data.frame(x = X, Trt = as.character(Trt_levels[2]))

        Trteff0 = unname(predict(rfmod, CbX0_rf, type = prd_type))
        Trteff1 = unname(predict(rfmod, CbX1_rf, type = prd_type))
      }

      if (type == "binary")
      {
        prd_type <- "prob"
        rfmod = randomForest(y ~ ., data = dat, family = binomial())

        Trt_levels <- sort(unique(Trt))

        CbX0_rf <- data.frame(x = X, Trt = as.character(Trt_levels[1]))
        CbX1_rf <- data.frame(x = X, Trt = as.character(Trt_levels[2]))

        Trteff0 = unname(predict(rfmod, CbX0_rf, type = prd_type)[,2])
        Trteff1 = unname(predict(rfmod, CbX1_rf, type = prd_type)[,2])
      }



    }

    Y_adj = Y - (1 - P) * Trteff1 - P * Trteff0

    est_cont = ori_Trt * Y_adj / P - (1 - ori_Trt) * Y_adj / (1 - P)
  } else { ## if no augmentation:

    est_cont = ori_Trt * Y / P - (1 - ori_Trt) * Y / (1 - P)
  }

  return(est_cont)
}
