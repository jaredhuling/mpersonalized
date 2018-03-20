genlassoD<-function(lambda2, lambda3, q){
  H = NULL

  for (i in 2:q)
    for (j in 1:(i-1)){
      newrow = numeric(q)
      newrow[i] = -1
      newrow[j] = 1
      H = rbind(H, newrow)
    }

  blockD = lambda3 * H

  if (lambda2 != 0){
    blockD = rbind(diag(rep(lambda2, q)), blockD)
  }

  return(blockD)
}


soft.thresh.scalar <- function(a, kappa) {
  pmax(0, a - kappa) - pmax(0, -a - kappa)
}

soft.thresh.vector <- function(a, kappa) {
  a*max(0,1-kappa/sqrt(sum(a^2)))
}


admm_optim = function(x, y, p, q, lambda1, lambda2, lambda3,
                      abs.tol = 1e-5, rel.tol = 1e-5, maxit = 500L, rho = NULL){

  blockD = genlassoD(lambda2 = lambda2, lambda3 = lambda3, q = q)#generalized lasso matrix D
  ngencon = nrow(blockD)#number of rows of generalized lasso constraint

  #whether group lasso is involved or not
  if (lambda1 == 0){
    #algorithm only involving a generalized lasso term, for methods without group lasso term
    #e.g., lasso + fused, fused

    blockA = blockD
    blockB = -diag(ngencon)

    A = bdiag(replicate(p, blockA, simplify = FALSE))
    B = bdiag(replicate(p, blockB, simplify = FALSE))

    xtx <- crossprod(x); xty <- crossprod(x, y)
    AtA <- crossprod(A); AtB <- crossprod(A, B)
    iters <- maxit

    if (is.null(rho)){

      eigs = eigen(xtx, symmetric=TRUE, only.values=TRUE)
      rho <- sqrt(eigs$values[1] * eigs$values[p * q])

    }

    beta <- numeric(p * q)
    gamma<- numeric(ngencon * p)
    npen = dim(blockA)[1]
    nu <- numeric(npen * p)
    U <- as(xtx + rho * AtA, "Matrix")

    for (i in 1:maxit) {
      v <- xty - rho * (crossprod(A, nu) + AtB%*%gamma)
      beta <- as.vector(solve(U, v))

      nu.prev = nu; gamma.prev = gamma

      #update paramter blockwisely
      for (j in 1:p){
        blockbeta = beta[((j - 1) * q + 1):(j * q)]
        blocknu = nu[((j - 1) * npen + 1):(j * npen)]

        z1 <- blockD %*% blockbeta - crossprod(blockB, blocknu)
        blockgamma <- soft.thresh.scalar(z1, 1 / rho)

        blocknu <- blocknu + (blockA %*% blockbeta + blockB %*% blockgamma)

        nu[((j - 1) * npen + 1):(j * npen)] = blocknu
        gamma[((j - 1) * ngencon + 1):(j * ngencon)] = blockgamma
      }

      r_norm = sqrt(sum((nu - nu.prev) ^ 2))
      s_norm = sqrt(sum((rho * AtB %*% (gamma - gamma.prev)) ^ 2))
      eps_pri = sqrt(nrow(A)) * abs.tol + rel.tol * max(sqrt(sum((A %*% beta)^2)), sqrt(sum((B %*% gamma)^2)))
      eps_dual = sqrt(ncol(A)) * abs.tol + rel.tol * sqrt(sum((crossprod(A,nu))^2 ))

      if (r_norm < eps_pri & s_norm < eps_dual) {
        iters <- i
        break
      }
    }

    admm_sol = list(beta = beta, gamma = gamma, iters = iters)

  } else {
    #algorithm involving both a group lasso term and a generalized lasso term, for methods with a group lasso
    #e.g.,sparse group lasso + fused, group lasso + fused

    blockA = rbind(blockD, diag(q))
    blockB = rbind(-diag(ngencon), matrix(0, nrow = q, ncol = ngencon))
    blockC = rbind(matrix(0, nrow = ngencon, ncol = q), -diag(q))

    A = bdiag(replicate(p, blockA, simplify = FALSE))
    B = bdiag(replicate(p, blockB, simplify = FALSE))
    C = bdiag(replicate(p, blockC, simplify = FALSE))

    xtx <- crossprod(x); xty <- crossprod(x, y)
    AtA <- crossprod(A); AtB <- crossprod(A, B); AtC <- crossprod(A, C)
    iters <- maxit

    if (is.null(rho)){

      eigs = eigen(xtx, symmetric=TRUE, only.values=TRUE)
      rho <- sqrt(eigs$values[1] * eigs$values[p * q])

    }

    beta <- numeric(p * q)
    gamma<- numeric(ngencon * p)
    eta <- numeric(p * q)
    npen = dim(blockA)[1]
    nu <- numeric(npen * p)
    U <- as(xtx + rho * AtA, "Matrix")

    for (i in 1:maxit) {
      v <- xty - rho * (crossprod(A, nu) + AtB%*%gamma+AtC%*%eta)
      beta <- as.vector(solve(U, v))

      nu.prev = nu; gamma.prev = gamma; eta.prev = eta;

      #update paramter blockwisely
      for (j in 1:p){
        blockbeta = beta[((j - 1) * q + 1):(j * q)]
        blocknu = nu[((j - 1) * npen + 1):(j * npen)]

        z1 <- blockD %*% blockbeta - crossprod(blockB, blocknu)
        blockgamma <- soft.thresh.scalar(z1, 1 / rho)

        z2 <- blockbeta - crossprod(blockC, blocknu)
        blocketa <- soft.thresh.vector(z2, lambda1 / rho)

        blocknu <- blocknu + (blockA %*% blockbeta + blockB %*% blockgamma + blockC %*% blocketa)

        nu[((j - 1) * npen + 1):(j * npen)] = blocknu
        gamma[((j - 1) * ngencon + 1):(j * ngencon)] = blockgamma
        eta[((j - 1) * q + 1):(j * q)] = blocketa
      }

      r_norm = sqrt(sum((nu - nu.prev) ^ 2))
      s_norm = sqrt(sum((rho * AtB %*% (gamma - gamma.prev) + rho * AtC %*% (eta - eta.prev))^2 ))
      eps_pri = sqrt(nrow(A)) * abs.tol + rel.tol * max(sqrt(sum((A %*% beta)^2)), sqrt(sum((B %*% gamma + C %*% eta)^2)))
      eps_dual = sqrt(ncol(A)) * abs.tol + rel.tol * sqrt(sum( (crossprod(A,nu))^2 ))

      if (r_norm < eps_pri & s_norm < eps_dual) {
        iters <- i
        break
      }
    }

    admm_sol = list(beta = beta, gamma = gamma, eta = eta, iters = iters)
  }

  return(admm_sol)
}
