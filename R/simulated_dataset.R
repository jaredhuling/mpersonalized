#' @title Simulated Dataset Generator
#'
#' @description This function could generate a simulated dataset for the test and usage of this package.
#' @param n sample size
#' @param sim_seed the seed supplied to the generation
#' @param problem a character specifiy whether you want to solve "meta-analysis" or "multiple outcomes" problem.
#' @return A simulated dataset.
#' \item{Xlist}{a list object with \eqn{k}th element denoting the covariate matrix of study k}
#' \item{Ylist}{a list object with \eqn{k}th element denoting the response vector of study k}
#' \item{Trtlist}{a list object with \eqn{k}th element denoting the treatment vector of study k (coded as 0 or 1)}
#' @export

simulated_dataset = function(n, sim_seed,
                             problem = c("meta-analysis", "multiple outcomes")){

  q = 6; P = 0.5
  p = 50
  set.seed(1)
  B=matrix(0,nrow=q,ncol=2*p+2)
  B[,1:12]=matrix(rep(2*sign(rnorm(12)),q),byrow=TRUE,nrow=q);
  B[,1:12]=B[,1:12]+0.2*sign(rnorm(12*q))
  B[,(p+2):(p+7)]=matrix(rep(2*sign(rnorm(6)),q),byrow=TRUE,nrow=q);
  B[,p+2]=B[,p+2]*0.5
  B[,(p+2):(p+7)]=B[,(p+2):(p+7)]+0.5*sign(rnorm(6*q))
  set.seed(NULL)

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
    set.seed(NULL)

    return(list(Xlist = Xlist, Ylist = Ylist, Trtlist = Trtlist))
  }


  if (problem == "multiple outcomes"){
    Ylist = replicate(q, list())
    Trtlist = replicate(q, list())

    set.seed(sim_seed)
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
    set.seed(NULL)

    return(list(X = X[,2 : (p + 1)], Ylist = Ylist, Trt = Trt))
  }
}
