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

      if (type == "continuous")
      {
        dat = data.frame(y = Y, x = X, Trt = as.character(Trt))

        rfmod = randomForest(y ~ ., data = dat)
        prd_type <- "response"

        Trt_levels <- sort(unique(Trt))

        CbX0_rf <- data.frame(x = X, Trt = as.character(Trt_levels[1]))
        CbX1_rf <- data.frame(x = X, Trt = as.character(Trt_levels[2]))

        Trteff0 = unname(predict(rfmod, CbX0_rf, type = prd_type))
        Trteff1 = unname(predict(rfmod, CbX1_rf, type = prd_type))
      }

      if (type == "binary")
      {
        dat = data.frame(y = as.factor(Y), x = X, Trt = as.character(Trt))
        prd_type <- "prob"
        rfmod = randomForest(y ~ ., data = dat)

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





contrast_builder_joint = function(Xlist, Ylist, ori_Trtlist, Plist,
                                  eff_aug = TRUE,
                                  response_model = c("SGL+SL"),
                                  type = c("continuous", "binary"),
                                  contrast_builder_folds = 10)
{


  type = match.arg(type)
  response_model = match.arg(response_model)

  Trtlist <- lapply(ori_Trtlist, function(ori_Trt) 2 * (ori_Trt - 0.5))

  q <- length(Xlist)



  if (eff_aug == TRUE)
  {

    CbXlist <- CbX0list <- CbX1list <- vector(mode = "list", length = q)

    for (s in 1:q)
    {
      IntX <- sweep(Xlist[[s]], 1, Trtlist[[s]], '*')
      CbXlist[[s]]  <- cbind(Xlist[[s]], Trtlist[[s]], IntX)

      n = dim(Xlist[[s]])[1]; p = dim(Xlist[[s]])[2]

      CbX0list[[s]] = cbind(Xlist[[s]], rep(-1, n), -Xlist[[s]])
      CbX1list[[s]] = cbind(Xlist[[s]], rep(1, n),   Xlist[[s]])
    }




    if (response_model == "SGL+SL")
    {
      if (type == "continuous")
      {
        sgmod = fit_sparse_group_fused_lasso_cv(modelYlist = Ylist,
                                                modelXlist = CbXlist, type = "linear",
                                                cv_folds = contrast_builder_folds, nlambda = 25)

        Trteff0list <- Trteff1list <- vector(mode = "list", length = q)

        for (s in 1:q)
        {
          intercept_s <- sgmod$intercept[s]
          beta_s      <- sgmod$beta[s,]

          Trteff0list[[s]] <- intercept_s + drop(CbX0list[[s]] %*% beta_s)
          Trteff1list[[s]] <- intercept_s + drop(CbX1list[[s]] %*% beta_s)
        }

      } else if (type == "binary")
      {
        sgmod = fit_sparse_group_fused_lasso_cv(modelYlist = Ylist,
                                                modelXlist = CbXlist, type = "logit",
                                                cv_folds = contrast_builder_folds, nlambda = 25)

        Trteff0list <- Trteff1list <- vector(mode = "list", length = q)

        for (s in 1:q)
        {
          intercept_s <- sgmod$intercept[s]
          beta_s      <- sgmod$beta[s,]

          xbeta0 <- intercept_s + drop(CbX0list[[s]] %*% beta_s)
          xbeta1 <- intercept_s + drop(CbX1list[[s]] %*% beta_s)

          Trteff0list[[s]] <- 1 / (1 + exp(-xbeta0))
          Trteff1list[[s]] <- 1 / (1 + exp(-xbeta1))
        }

      }
    }

    est_cont_list <- vector(mode= "list", length = q)
    for (s in 1:q)
    {
      Y_adj = Ylist[[s]] - (1 - Plist[[s]]) * Trteff1list[[s]] - Plist[[s]] * Trteff0list[[s]]

      est_cont_list[[s]] = ori_Trtlist[[s]] * Y_adj / Plist[[s]] - (1 - ori_Trtlist[[s]]) * Y_adj / (1 - Plist[[s]])
    }

  } else { ## if no augmentation:

    #est_cont = ori_Trt * Y / P - (1 - ori_Trt) * Y / (1 - P)
    est_cont_list <- vector(mode= "list", length = q)
    for (s in 1:q)
    {
      est_cont_list[[s]] <- ori_Trtlist[[s]] * Ylist[[s]] / Plist[[s]] - (1 - ori_Trtlist[[s]]) * Ylist[[s]] / (1 - Plist[[s]])
    }
  }

  return(est_cont_list)
}
