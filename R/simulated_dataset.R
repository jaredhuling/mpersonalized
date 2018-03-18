#' @title Simulated Dataset Generator
#'
#' @description Generate a simulated dataset, which could be used to demonstrate the
#' features of the mpersonalized package.
#' @param n Sample size of the simulated dataset
#' @param problem A character string specified what problem the simulated dataset is
#' generated for. \code{problem} can be set to "meta-analysis" or "multiple outcomes".
#' @details In the simulated dataset, outcomes are generated from the model
#' \deqn{Y = \delta_0 + \bm{X} \bm{\delta} + A (\theta_0 + \bm{X}\bm{\theta})+\epsilon,}
#' where \eqn{\bm{X}} is the baseline covariates and \eqn{A} is the treatment indicator coded
#' as 0,1. For different outcomes or studies, values of \eqn{\delta_0}, \eqn{\bm{\delta}},
#' \eqn{\theta_0} and \eqn{\bm{\theta}} are also different so as to represent the heterogeneity
#' in real problems.
#'
#' The number of different studies/outcomes is set to be 6 and total number of candidate
#' covariates is 50. Treatment indicator \eqn{A} is generated with equal probability of
#' 0 or 1.
#'
#' This function randomly generates the coefficients for each study/outcome and then
#' generates the baseline covariates and error term for each subject. Depending on the
#' value of \code{problem}, generation of baseline covariates are slightly different.
#' For \code{problem = "meta-analysis"}, baseline covariates are generated independently
#' for each study; for \code{problem = "multiple outcomes"}, baseline covariates are
#' the same across different outcomes.
#'
#' @return A list object of the ingredients from the simulated dataset. The elements of
#' this list depends on value of \code{problem}.
#'
#' For \code{problem = "meta-analysis"},
#' \item{Xlist}{a list object with \eqn{k}th element denoting the baseline covariate matrix of \eqn{k}th study }
#' \item{Ylist}{a list object with \eqn{k}th element denoting the response vector of \eqn{k}th study}
#' \item{Trtlist}{a list object with \eqn{k}th element denoting the treatment vector of \end{k}th study and coded as 0 or 1}
#' \item{B}{the coefficient matrix containing \eqn{\delta_0}, \eqn{\bm{\delta}},
#' \eqn{\theta_0} and \eqn{\bm{\theta}}}
#'
#' For \code{problem = "multiple outcomes"},
#' \item{X}{a matrix object denoting the baseline covariate matrix}
#' \item{Ylist}{a list object with \eqn{k}th element denoting the response vector of \eqn{k}th outcome}
#' \item{Trt}{a vector denoting the treatment and coded as 0 or 1)}
#' \item{B}{the coefficient matrix containing \eqn{\delta_0}, \eqn{\bm{\delta}},
#' \eqn{\theta_0} and \eqn{\bm{\theta}}}
#' @export

simulated_dataset = function(n,
                             problem = c("meta-analysis", "multiple outcomes")){

  q = 6; P = 0.5
  p = 50
  B = matrix(0,nrow = q,ncol = 2*p+2)
  B[,1:12] = matrix(rep(2*sign(rnorm(12)),q),byrow = TRUE,nrow = q);
  B[,1:12] = B[,1:12]+0.2*sign(rnorm(12*q))
  B[,(p+2):(p+7)] = matrix(rep(2*sign(rnorm(6)),q),byrow = TRUE,nrow = q);
  B[,p+2] = B[,p+2]*0.5
  B[,(p+2):(p+7)] = B[,(p+2):(p+7)]+0.5*sign(rnorm(6*q))

  if (problem  == "meta-analysis"){
    Xlist = replicate(q, list())
    Ylist = replicate(q, list())
    Trtlist = replicate(q, list())
    set.seed(sim_seed)
    for (j in 1:q){
      X=matrix(rnorm(n*p,0,1),nrow=n)
      X=cbind(rep(1,n),X) #Add Intercept
      Trt=rbinom(n,1,P) #Assignment
      IntX=matrix(NA,nrow=n,ncol=p+1)
      for (i in 1:n)
        IntX[i,]=X[i,]*Trt[i]
      CbX=cbind(X,IntX) #Combined design matrix including interaction term

      Y=CbX%*%B[j,]+rnorm(n,0,2)

      Xlist[[j]] = X[,2 : (p + 1)]; Ylist[[j]] = Y; Trtlist[[j]] = Trt
    }

    return(list(Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist, B = B))
  }


  if (problem == "multiple outcomes"){
    Ylist = replicate(q, list())
    Trtlist = replicate(q, list())

    X=matrix(rnorm(n*p,0,1),nrow=n)
    X=cbind(rep(1,n),X) #Add Intercept
    Trt=rbinom(n,1,P) #Assignment
    IntX=matrix(NA,nrow=n,ncol=p+1)
    for (i in 1:n)
      IntX[i,]=X[i,]*Trt[i]
    CbX=cbind(X,IntX) #Combined design matrix including interaction term

    for (j in 1:q){
      Y=CbX%*%B[j,]+rnorm(n,0,2)

      Ylist[[j]] = Y
    }

    return(list(X = X[,2 : (p + 1)], Ylist = Ylist, Trt = Trt, B = B))
  }
}
