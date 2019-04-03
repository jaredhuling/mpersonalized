genlassoD<-function(lambda2,lambda3,q){
  H = NULL

  for (i in 2:q)
    for (j in 1:(i-1)){
      newrow=numeric(q)
      newrow[i]=-1
      newrow[j]=1
      H=rbind(H,newrow)
    }

  blockD=lambda3*H
  blockD=rbind(diag(rep(lambda2,q)),blockD)

  return(blockD)
}

conblockA<-function(blockD){
  q=ncol(blockD)
  ngencon=nrow(blockD)
  blockA=rbind(blockD,diag(q))

  return(blockA)
}

conblockB<-function(blockD){
  q=ncol(blockD)
  ngencon=nrow(blockD)
  blockB=rbind(-diag(ngencon),matrix(0,nrow=q,ncol=ngencon))

  return(blockB)
}

conblockC<-function(blockD){
  q=ncol(blockD)
  ngencon=nrow(blockD)
  blockC=rbind(matrix(0,nrow=ngencon,ncol=q),-diag(q))

  return(blockC)
}


admm_meta <- function(x, y, p,q, lambda1, lambda2,lambda3, rho, abs.tol = 1e-5, rel.tol = 1e-5, maxit = 500L) {

  blockD=genlassoD(lambda2 = lambda2, lambda3 = lambda3, q=q)

  n <- nrow(x)
  ngencon <- nrow(blockD)
  blockA=conblockA(blockD);blockB=conblockB(blockD);blockC=conblockC(blockD)
  Alist=replicate(p, list());Blist=replicate(p, list());Clist=replicate(p, list());Dlist=replicate(p,list());
  for (i in 1:p){
    Alist[[i]]=blockA;Blist[[i]]=blockB;Clist[[i]]=blockC;Dlist[[i]]=blockD;
  }
  A=bdiag(Alist);B=bdiag(Blist);C=bdiag(Clist);D=bdiag(Dlist)


  xtx <- crossprod(x)
  AtA <- crossprod(A)
  xty <- crossprod(x,y)
  AtB <- crossprod(A,B)
  AtC <- crossprod(A,C)
  DtD <- as(crossprod(D),"Matrix")
  #lambda <- lambda * n
  iters <- maxit

  ## if rho value is not supplied,
  ## compute one that is good
  #   if (is.null(rho)) {
  #     eigs <- eigs_sym(xtx, k = 2,
  #                      which = "BE",
  #                      opts = list(maxitr = 500,
  #                                  tol = 1e-4))$values
  #   eigs=eigen(xtx,symmetric=TRUE,only.values=TRUE)
  #   rho <- eigs$values[1]^(1/ 3) * lambda1^(2/3)


  beta <- numeric(p*q)
  gamma<- numeric(ngencon*p)
  eta<-numeric(p*q)
  npen=ngencon+q
  nu <- numeric(npen*p)
  U <- as(xtx + rho * AtA, "Matrix")

  for (i in 1:maxit) {
    v <- xty - rho * (crossprod(A,nu)+AtB%*%gamma+AtC%*%eta)
    beta <- as.vector(solve(U, v))

    nu.prev=nu; gamma.prev=gamma; eta.prev=eta;

    for (j in 1:p){
      blockbeta=beta[((j-1)*q+1):(j*q)]
      blocknu=nu[((j-1)*npen+1):(j*npen)]

      z1 <- blockD%*%blockbeta-crossprod(blockB,blocknu)
      blockgamma <- soft.thresh.scalar(z1, 1/rho)

      z2 <- blockbeta-crossprod(blockC,blocknu)
      blocketa <- soft.thresh.vector(z2,lambda1/rho)

      blocknu <- blocknu + (blockA %*% blockbeta + blockB %*% blockgamma + blockC %*% blocketa)

      nu[((j-1)*npen+1):(j*npen)]=blocknu
      gamma[((j-1)*ngencon+1):(j*ngencon)]=blockgamma
      eta[((j-1)*q+1):(j*q)]=blocketa
    }

    r_norm = sqrt(sum( (nu-nu.prev)^2 ))
    s_norm = sqrt(sum( (rho * AtB %*% (gamma - gamma.prev)+rho * AtC %*% (eta-eta.prev))^2 ))
    eps_pri = sqrt(n) * abs.tol + rel.tol * max(sqrt(sum((A %*% beta)^2)), sqrt(sum((B %*% gamma)^2)), sqrt(sum((C %*% eta)^2)) )
    eps_dual = sqrt(p*q) * abs.tol + rel.tol * sqrt(sum( (crossprod(A,nu))^2 ))


    if (r_norm < eps_pri & s_norm < eps_dual) {
      iters <- i
      break
    }
  }

  list(beta = beta, gamma = gamma, eta = eta, iters = iters)
}


soft.thresh.scalar <- function(a, kappa) {
  pmax(0, a - kappa) - pmax(0, -a - kappa)
}

soft.thresh.vector <- function(a, kappa) {
  a*max(0,1-kappa/sqrt(sum(a^2)))
}

