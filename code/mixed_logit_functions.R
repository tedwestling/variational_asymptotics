
log1pex <- function(x) {
  ret <- rep(NA, length(x))
  large <- (x > 13)
  ret[large] <- x[large] 
  ret[!large] <- log1p(exp(x[!large]))
  ret
}

expit <- function(x, d=0) {
  ex <- 1/(1 + exp(-x))
  if(d==0) return(ex)
  if(d==1) return(ex - ex^2)
  if(d==2) return(ex - 3 * ex^2 + 2 * ex^3 )
  if(d==3) return(ex - 7 * ex^2 + 12 * ex^3 - 6 * ex^4)
  if(d==4) return(ex - 15 * ex^2 + 50 * ex^3 - 60 * ex^4 + 24 * ex^5)
}

logit <- function(x) log(x) - log(1-x)

mean01 <- function(f) {
  g <- function(x) sapply(x, function(xi) f(sqrt(2) * xi))
  (1/sqrt(pi)) * ghQuad(g, rule100) 
}

one_elbo <- function(beta, Sigma, mi, Di, Li, Yi, Zi, Xi) {
  Si <- Li %*% diag(exp(Di), nrow=length(Di)) %*% t(Li)
  delta <- c(Xi %*% beta + Zi %*% mi)
  omega <- sqrt(sapply(1:nrow(Zi), function(j) c(Zi[j, ] %*% Si %*% Zi[j,])))
  sum(Yi * delta) -
    mean01(function(gamma) sum(log1pex(delta + omega * gamma))) - 
    .5 * log(det(Sigma)) + 
    .5 * sum(Di) - 
    .5 * sum(diag(solve(Sigma) %*% (Si + mi %*% t(mi))))
}

elbo <- function(beta, Sigma, m, D, L, data, verbose=FALSE, nodes=1) {
  N <- length(data)
  if(nodes > 1) {
    require(parallel)
    elbos <- mclapply(1:N, function(i) {
      if(verbose & i %% 1000 == 0) cat("Obs", i, "of", length(data), "\n")
      one_elbo(beta, Sigma, mi=m[i,], Di=D[i,], Li=L[[i]], Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
    }, mc.cores=nodes)
    elbos <- unlist(elbos)
  } else {
    elbos <- sapply(1:N, function(i) {
      if(verbose & i %% 1000 == 0) cat("Obs", i, "of", length(data), "\n")
      one_elbo(beta, Sigma, mi=m[i,], Di=D[i,], Li=L[[i]], Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
    })
  }
  sum(elbos)
}

one_elbo_grad_beta <- function(beta, Sigma, mi, Di, Li, Yi, Zi, Xi) {
  Si <- Li %*% diag(exp(Di), nrow=length(Di)) %*% t(Li)
  delta <- c(Xi %*% beta + Zi %*% mi)
  omega <- sqrt(sapply(1:nrow(Zi), function(j) c(Zi[j, ] %*% Si %*% Zi[j,])))
  c(t(Xi) %*% Yi - sapply(1:ncol(Xi), function(k) mean01(function(gamma) sum(Xi[,k] * expit(delta + omega * gamma)))))
}

elbo_grad_beta <- function(beta, Sigma, m, D, L, data, nodes=1) {
  N <- length(data)
  if(nodes > 1) {
    require(parallel)
    grad_list <- mclapply(1:N, function(i) {
      one_elbo_grad_beta(beta, Sigma, mi=m[i,], Di=D[i,], Li=L[[i]], Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
    }, mc.cores=nodes)
    grads <- matrix(unlist(grad_list), nrow=N, byrow=TRUE)
  } else {
    grads <- t(sapply(1:N, function(i) {
      one_elbo_grad_beta(beta, Sigma, mi=m[i,], Di=D[i,], Li=L[[i]], Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
    }))
  }
  colSums(grads)
}

one_elbo_hess_beta <- function(beta, Sigma, mi, Di, Li, Yi, Zi, Xi) {
  Si <- Li %*% diag(exp(Di), nrow=length(Di)) %*% t(Li)
  delta <- c(Xi %*% beta + Zi %*%  mi)
  omega <- sqrt(sapply(1:nrow(Zi), function(j) c(Zi[j, ] %*% Si %*% Zi[j,])))
  -sapply(1:ncol(Xi), function(k) {
    sapply(1:ncol(Xi), function(l) {
      mean01(function(gamma) sum(Xi[,k] * Xi[,l] * expit(delta + omega * gamma, 1)))
    })
  })
}

elbo_hess_beta <- function(beta, Sigma, m, D, L, data) {
  N <- length(data)
  hesses <- t(sapply(1:N, function(i) {
    c(one_elbo_hess_beta(beta, Sigma, mi=m[i,], Di=D[i,], Li=L[[i]], Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x))
  }))
  matrix(colSums(hesses), nrow=length(beta))
}

beta_newton_raphson <- function(beta, Sigma, m, D, L, data) {
  K <- length(beta)
  elbo_new <- elbo(beta, Sigma, m, D, L, data)
  grad <- elbo_grad_beta(beta, Sigma, m, D, L, data)
  while(TRUE) {
    elbo_old <- elbo_new
    hess <- elbo_hess_beta(beta, Sigma, m, D, L, data)
    beta <- beta - solve(hess) %*% grad
    grad <- elbo_grad_beta(beta, Sigma, m, D, L, data)
    elbo_new <- elbo(beta, Sigma, m, D, L, data)
    if(abs(elbo_old - elbo_new) / abs(elbo_old) < .0001 & sum(abs(grad)) < .0001) break
  }
  return(list(beta=beta, elbo=elbo_new))
}

one_elbo_grad_psi <- function(beta, Sigma, mi, Di, Li, Yi, Zi, Xi) {
  Si <- Li %*% diag(exp(Di), nrow=length(Di)) %*% t(Li)
  delta <- c(Xi %*% beta + Zi %*% mi)
  omega <- sqrt(sapply(1:nrow(Zi), function(j) c(Zi[j, ] %*% Si %*% Zi[j,])))
  grad_mi <- c(t(Zi) %*% Yi - sapply(1:ncol(Zi), function(k) mean01(function(gamma) sum(Zi[,k] * expit(delta + omega * gamma))))) - solve(Sigma) %*% mi
  
  core <- -.5 * sapply(1:ncol(Zi), function(k) {
    sapply(1:ncol(Zi), function(l) {
      mean01(function(gamma) sum(Zi[,k] * Zi[,l] * expit(delta + omega * gamma, 1)))
    })
  }) - .5 * solve(Sigma) 
  
  grad_Di <- t(Li) %*% core %*% Li + .5 * diag(exp(-Di), nrow=length(Di)) 
  grad_Li <- 2 * core %*% Li %*% diag(exp(Di), nrow=length(Di))
  c(grad_mi, exp(Di) * diag(grad_Di), grad_Li[lower.tri(grad_Li, diag=FALSE)])
}

ind <- function(i, D) {
  vec <- rep(0, D)
  vec[i] <- 1
  vec
}



psi_opt_fun <- function(psi, beta, Sigma, Yi, Zi, Xi) {
  P <- ncol(Zi)
  Li <- matrix(0, nrow=P, ncol=P)
  diag(Li) <- 1
  Li[lower.tri(Li, diag=FALSE)] <- psi[-(1:(2*P))]
  one_elbo(beta, Sigma, mi=psi[1:P], Di= psi[(P + 1):(2 * P)], Li=Li, Yi, Zi, Xi)
}

psi_grad_fun <- function(psi, beta, Sigma, Yi, Zi, Xi) {
  P <- ncol(Zi)
  Li <- matrix(0, nrow=P, ncol=P)
  diag(Li) <- 1
  Li[lower.tri(Li, diag=FALSE)] <- psi[-(1:(2*P))]
  one_elbo_grad_psi(beta, Sigma, mi=psi[1:P], Di= psi[(P + 1):(2 * P)], Li=Li, Yi, Zi, Xi)
}


optimize_elbo <- function(beta, Sigma, m=NULL, D=NULL, L=NULL, data, verbose=TRUE, plot=FALSE, tol=.0001,nodes=1) {
  if(plot) {
    curve(expit(beta[1] + (beta[2] ) * x + (beta[3]) * x^2), ylim=c(0,1), lwd=3, ylab='prob')
    curve(expit(beta[4] + (beta[5] ) * x + (beta[6]) * x^2), ylim=c(0,1), lwd=3, col='red', add=TRUE)
  }
  n <- length(data)
  P <- ncol(data[[1]]$z)
  K <- length(beta)
  testL <- matrix(0, nrow=P, ncol=P)
  diag(testL) <- 1
  psi_stepsizes <- rep(1, n)
  if(is.null(m) | is.null(D) | is.null(L)) {
    m <- D <- matrix(NA, nrow=n, ncol=P)
    L <- lapply(1:n, function(i) testL)
    if(verbose) print("Initializing variational parameters...")
    
    if(nodes > 1) {
      require(parallel)
    
      updates <- mclapply(1:n, function(i) {
        if(verbose & i %% 1000 == 0) cat(i, " of ", n, " \n")
        psi_i <- rep(0, 2*P + P * (P-1)/2)
        psi_grad_i <- psi_grad_fun(psi_i, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
        fn <- function(psi) psi_opt_fun(psi, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x )
        psi_stepsize <- backtrack(fx = fn(psi_i), gx = psi_grad_i, t = 1.5*psi_stepsizes[i], f=fn, x = psi_i)
        #print(stepsize)
        
        psi_i <- psi_i + psi_stepsize * psi_grad_i
        list(i=i, psi_i=psi_i, stepsize=psi_stepsize)
      }, mc.cores=nodes)
      for(i in 1:n){
        m[i,] <- updates[[i]]$psi_i[1:P]#psi_i$par[1:P]
        D[i,] <- updates[[i]]$psi_i[(P+1):(2*P)] #psi_i$par[(P+1):(2*P)]
        L[[i]][lower.tri(testL, diag=FALSE)] <- updates[[i]]$psi_i[-(1:(2*P))] #psi_i$par[-(1:(2*P))]
        psi_stepsizes[i] <- updates[[i]]$stepsize
      }
    } else{
      for(i in 1:n) {
        if(verbose & i %% 1000 == 0) print(paste0(i, " of ", n))
        psi_i <- rep(0, 2*P + P * (P-1)/2)
        psi_grad_i <- psi_grad_fun(psi_i, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
        fn <- function(psi) psi_opt_fun(psi, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x )
        psi_stepsizes[i] <- backtrack(fx = fn(psi_i), gx = psi_grad_i, t = 1.5*psi_stepsizes[i], f=fn, x = psi_i)
        #print(stepsize)
        
        psi_i <- psi_i + psi_stepsizes[i] * psi_grad_i#optim(rep(0, 2*P + P * (P-1)/2), fn=psi_opt_fun, gr=psi_grad_fun, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x, method='BFGS', control=list(fnscale=-1))
        m[i,] <- psi_i[1:P]#psi_i$par[1:P]
        D[i,] <- psi_i[(P+1):(2*P)] #psi_i$par[(P+1):(2*P)]
        L[[i]][lower.tri(testL, diag=FALSE)] <- psi_i[-(1:(2*P))] #psi_i$par[-(1:(2*P))]
      }
    }
  }
  
  elbo_new <- elbo(beta, Sigma, m, D, L, data=data, nodes=nodes)
  beta_stepsize <- 1
  repeat {
    elbo_old <- elbo_new
    if(verbose) print(paste("elbo =", elbo_old))
    
    if(verbose) print("Updating beta...")
    #opt <- betamu_newton_raphson(rep(0,K), rep(0,P), Sigma, m, D, L, data)
    if(verbose) print("Computing gradient...")
    grad <- elbo_grad_beta(beta, Sigma, m, D, L, data, nodes=nodes)#beta_newton_raphson(beta, Sigma, m, D, L, data)
    fn <- function(beta) elbo(beta, Sigma=Sigma, m=m, D=D, L=L, data=data, nodes=nodes)
    if(verbose) print("Backtracking to find stepsize...")
    beta_stepsize <- backtrack(fx = fn(beta), gx = grad, t = 1.5*beta_stepsize, f=fn, x = beta, verbose=verbose)
    if(verbose) print(paste("beta stepsize = ", beta_stepsize))
    beta <- beta + beta_stepsize * grad #opt$beta
    if(verbose) print(paste("beta =",paste(beta, collapse = ", ")))
    if(plot) {s
      curve(expit(beta[1] + (beta[2] ) * x + (beta[3]) * x^2), ylim=c(0,1), lwd=3, ylab='prob')
      curve(expit(beta[4] + (beta[5] ) * x + (beta[6]) * x^2), ylim=c(0,1), lwd=3, col='red', add=TRUE)
    }
    
    if(verbose) print("Updating Sigma...")
    Sigma <- matrix(rowMeans(sapply(1:n, function(i) c(m[i,] %*% t(m[i,]) + L[[i]] %*% diag(exp(D[i,]), nrow=length(D[i,])) %*% t(L[[i]])))), nrow=P)
    if(verbose) { print("Sigma = "); print(Sigma) }
    
    if(verbose) print("Updating variational parameters...")    
    if(nodes > 1) {
      updates <- mclapply(1:n, function(i) {
        if(verbose & i %% 1000 == 0) cat(i, " of ", n, " \n")
        psi_i <- c(m[i,], D[i,], L[[i]][lower.tri(testL, diag=FALSE)])
        psi_grad_i <- psi_grad_fun(psi_i, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
        fn <- function(psi) psi_opt_fun(psi, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x )
        psi_stepsize <- backtrack(fx = fn(psi_i), gx = psi_grad_i, t = 1.5*psi_stepsizes[i], f=fn, x = psi_i)
        #print(stepsize)
        
        psi_i <- psi_i + psi_stepsizes[i] * psi_grad_i
        list(psi_i=psi_i, stepsize=psi_stepsize)
      }, mc.cores=nodes)
      for(i in 1:n){
        m[i,] <- updates[[i]]$psi_i[1:P]#psi_i$par[1:P]
        D[i,] <- updates[[i]]$psi_i[(P+1):(2*P)] #psi_i$par[(P+1):(2*P)]
        L[[i]][lower.tri(testL, diag=FALSE)] <- updates[[i]]$psi_i[-(1:(2*P))] #psi_i$par[-(1:(2*P))]
        psi_stepsizes[i] <- updates[[i]]$stepsize
      }
    } else{
      for(i in 1:n) {
        #if(verbose) if(i %% 100 == 0) print(paste0(i, " of ", n))
        psi_i <- c(m[i,], D[i,], L[[i]][lower.tri(testL, diag=FALSE)])
        psi_grad_i <- psi_grad_fun(psi_i, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
        fn <- function(psi) psi_opt_fun(psi, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x )
        psi_stepsizes[i] <- backtrack(fx = fn(psi_i), gx = psi_grad_i, t = psi_stepsizes[i], f=fn, x = psi_i)
        #print(stepsize)
        
        psi_i <- psi_i + psi_stepsizes[i] * psi_grad_i#optim(rep(0, 2*P + P * (P-1)/2), fn=psi_opt_fun, gr=psi_grad_fun, beta=beta, Sigma=Sigma, Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x, method='BFGS', control=list(fnscale=-1))
        m[i,] <- psi_i[1:P]#psi_i$par[1:P]
        D[i,] <- psi_i[(P+1):(2*P)] #psi_i$par[(P+1):(2*P)]
        L[[i]][lower.tri(testL, diag=FALSE)] <- psi_i[-(1:(2*P))] #psi_i$par[-(1:(2*P))]
      }
    }
    
    
    elbo_new <- elbo(beta, Sigma, m, D, L, data=data, nodes=nodes)
    if(abs(elbo_old - elbo_new) / abs(elbo_old) < tol) break
  }
  print(paste("final elbo:", elbo_new))
  return(list(beta=beta, Sigma=Sigma, m=m, D=D, L=L, elbo=elbo_new))
}

one_elbo_vec <- function(vec, Yi, Zi, Xi) {
  K <- ncol(Xi)
  P <- ncol(Zi)
  beta <- vec[1:K]
  Sigma_D <- vec[(K+1):(K+P)]
  Sigma_L <- diag(rep(1, P))
  Sigma_L[lower.tri(Sigma_L, diag=FALSE)] <- vec[(K+P+1):(K+P+P*(P-1)/2)]
  Sigma <- Sigma_L %*% diag(exp(Sigma_D)) %*% t(Sigma_L)
  
  mi <- vec[(K+P+P*(P-1)/2+1):(K+P+P*(P-1)/2+P)]
  Di <- vec[(K+P+P*(P-1)/2+P+1):(K+P+P*(P-1)/2+2*P)]
  Li <- diag(rep(1,P))
  Li[lower.tri(Li, diag=FALSE)] <- vec[(K+P+P*(P-1)/2+2*P+1):(K+P+P*(P-1)/2+2*P+P*(P-1)/2)]
  one_elbo(beta, Sigma, mi, Di, Li, Yi, Zi, Xi)
}

one_elbo_thetagrad <- function(beta, Sigma, mi, Di, Li, Yi, Zi, Xi) {
  Si <- Li %*% diag(exp(Di), nrow=length(Di)) %*% t(Li)
  
  grad_beta <- one_elbo_grad_beta(beta, Sigma, mi, Di, Li, Yi=Yi, Zi=Zi, Xi=Xi)
  Sigma_inv <- solve(Sigma)
  
  Sigma_L <- LDL(Sigma)$L
  Sigma_L_inv <- solve(Sigma_L)
  Sigma_D <- log(LDL(Sigma)$D)
  
  grad_Sigma_D <- -.5 + .5 * diag(diag(exp(-Sigma_D), nrow=length(Sigma_D)) %*% Sigma_L_inv %*% (Si + mi %*% t(mi)) %*% t(Sigma_L_inv))
  grad_Sigma_L <- (Sigma_inv %*% (Si + mi %*% t(mi)) %*% t(Sigma_L_inv))[lower.tri(Sigma_L, diag=FALSE)]
  
  # grad_Sigma <- -.5 * Sigma_inv + .5 * Sigma_inv %*% (Si + mi %*% t(mi)) %*% Sigma_inv
  # grad_Sigma <- 2 * grad_Sigma - diag(diag(grad_Sigma))
  # grad_Sigma <- grad_Sigma[lower.tri(grad_Sigma, diag=TRUE)]
  #c(grad_beta, grad_Sigma)
  c(grad_beta, grad_Sigma_D, grad_Sigma_L)
}



backtrack <- function(b=.7, alpha=.25, fx, gx, t, f, x, verbose=FALSE) {
  repeat {
    if(verbose) print(paste("Checking stepsize ", t))
    if(f(x + t * gx) > fx + alpha * t * sum(gx^2)) return(t)
    t <- b * t
  }
} 


par_vec_to_list <- function(vec, n, K, P) {
  beta <- vec[1:K]
  #Sigma_D <- vec[(K+1):(K+P)]
  #Sigma_L <- diag(rep(1, P))
  #Sigma_L[lower.tri(Sigma_L, diag=FALSE)] <- vec[(K+P+1):(K+P+P*(P-1)/2)]
  #Sigma <- Sigma_L %*% diag(exp(Sigma_D)) %*% t(Sigma_L)
  
  m <- D <- matrix(NA, nrow=n, ncol=P)
  L <- vector('list', n)
  Sigma <- matrix(0, nrow=P, ncol=P)
  startind <- K#K+P+P*(P-1)/2
  for(i in 1:n) {
    m[i,] <- vec[(startind+1):(startind+P)]
    D[i,] <- vec[(startind+P+1):(startind+2*P)]
    L[[i]] <- diag(rep(1,P))
    L[[i]][lower.tri(L[[i]], diag=FALSE)] <- vec[(startind+2*P+1):(startind+2*P+P*(P-1)/2)]
    startind <- startind+2*P+P*(P-1)/2
    Sigma <- Sigma + m[i,] %*% t(m[i,]) + L[[i]] %*% diag(exp(D[i,]), nrow=ncol(D)) %*% t(L[[i]])
  }
  Sigma <- Sigma / n
  
  
  list(beta=beta, Sigma=Sigma, m=m, D=D, L=L)
}

par_list_to_vec <- function(beta, m, D, L) {
  # Sigma_D <- log(LDL(Sigma)$D)
  # Sigma_L <- LDL(Sigma)$L
  c(beta, 
    #Sigma_D, Sigma_L[lower.tri(Sigma_L, diag=FALSE)],
    c(sapply(1:nrow(m), function(i) c(m[i,], D[i,], L[[i]][lower.tri(L[[i]], diag=FALSE)]))))
}

vec_elbo <- function(vec, data, nodes=1, verbose=FALSE) {
  n <- length(data)
  K <- ncol(data[[1]]$x)
  P <- ncol(data[[1]]$z)
  l <- par_vec_to_list(vec, n, K, P)
  # print(l$beta)
  # print(l$Sigma)
  elbo(l$beta, l$Sigma, l$m, l$D, l$L, data=data, verbose=verbose, nodes=nodes)
}

vec_elbo_grad <- function(vec, data, nodes=1, verbose=FALSE) {
  n <- length(data)
  K <- ncol(data[[1]]$x)
  P <- ncol(data[[1]]$z)
  l <- par_vec_to_list(vec, n, K, P)
  
  beta_grad <- elbo_grad_beta(l$beta, l$Sigma, l$m, l$D, l$L, data=data, nodes=nodes)
  all_grads <- mclapply(1:n, function(i) {
    if(verbose & i %% 1000 == 0) cat("Obs", i, "of", length(data), "\n")
    #th_grad <- one_elbo_thetagrad(l$beta, l$Sigma, l$m[i,], l$D[i,], l$L[[i]], data[[i]]$y, data[[i]]$z, data[[i]]$x)
    one_elbo_grad_psi(l$beta, l$Sigma, l$m[i,], l$D[i,], l$L[[i]], data[[i]]$y, data[[i]]$z, data[[i]]$x)
    #c(th_grad, psi_grad)
  }, mc.cores=nodes)
  all_grads <- matrix(unlist(all_grads), nrow=length(data), byrow=TRUE)
  #theta_ln <- K + P  + P * (P-1)/2
  #c(colSums(all_grads[,1:theta_ln]),c(t(all_grads[,-(1:theta_ln)])))
  c(beta_grad, c(t(all_grads)))
}

one_elbo_theta0 <- function(psi, beta, Sigma, Yi, Zi, Xi) {
  #print(psi)
  K <- ncol(Xi)
  P <- ncol(Zi)
  mi <- psi[1:P]
  Di <- psi[(P+1):(2*P)]
  Li <- diag(rep(1,P))
  Li[lower.tri(Li, diag=FALSE)] <- psi[-(1:(2*P))]
  one_elbo(beta, Sigma, mi, Di, Li, Yi, Zi, Xi)
}


one_elbo_grad_theta0 <- function(psi, beta, Sigma, Yi, Zi, Xi) {
  K <- ncol(Xi)
  P <- ncol(Zi)
  mi <- psi[1:P]
  Di <- psi[(P+1):(2*P)]
  Li <- diag(rep(1,P))
  Li[lower.tri(Li, diag=FALSE)] <- psi[-(1:(2*P))]
  one_elbo_grad_psi(beta, Sigma, mi, Di, Li, Yi, Zi, Xi)
}

one_elbo_thetapsi <- function(thetapsi, Yi, Zi, Xi) {
  K <- ncol(Xi)
  P <- ncol(Zi)
  theta_ln <- K + P * (P+1) / 2
  psi_ln <- 2*P + P * (P-1)/2
  beta <- thetapsi[1:K]
  Sigma_D <- diag(exp(thetapsi[(K+1):(K+P)]), nrow=P)
  Sigma_L <- diag(rep(1, P))
  Sigma_L[lower.tri(Sigma_L, diag=FALSE)] <- thetapsi[(K+P+1):(theta_ln)]
  Sigma <- Sigma_L %*% Sigma_D %*% t(Sigma_L)
  #Sigma[upper.tri(Sigma, diag=FALSE)] <- Sigma[lower.tri(Sigma, diag=FALSE)]
  mi <- thetapsi[(theta_ln+1):(theta_ln + P)]
  Di <- thetapsi[(theta_ln + P+1):(theta_ln + 2*P)]
  Li <- matrix(0, P, P)
  diag(Li) <- 1
  Li[lower.tri(Li, diag=FALSE)] <- thetapsi[(theta_ln + 2*P + 1):(theta_ln + psi_ln)]
  one_elbo(beta, Sigma, mi, Di, Li, Yi, Zi, Xi)
}

one_elbo_theta <- function(theta, mi, Di, Li, Yi, Zi, Xi) {
  K <- ncol(Xi)
  P <- ncol(Zi)
  theta_ln <- K + P * (P+1) / 2
  psi_ln <- 2*P + P * (P-1)/2
  beta <- theta[1:K]
  Sigma_D <- diag(exp(theta[(K+1):(K+P)]), nrow=P)
  Sigma_L <- diag(rep(1, P))
  Sigma_L[lower.tri(Sigma_L, diag=FALSE)] <- theta[(K+P+1):(theta_ln)]
  Sigma <- Sigma_L %*% Sigma_D %*% t(Sigma_L)
  #Sigma[upper.tri(Sigma, diag=FALSE)] <- Sigma[lower.tri(Sigma, diag=FALSE)]
  one_elbo(beta, Sigma, mi, Di, Li, Yi, Zi, Xi)
}


optimize_psi <- function(beta, Sigma, data, psi_start, nodes=1, verbose=FALSE) {
  P <- nrow(Sigma)
  psi_hats <- mclapply(1:length(data), function(i) {
    if(i %% 500 == 0 & verbose) cat("Obs", i, "of", Nsub, "\n")
    one_data <- data[[i]]
    this_opt <- optim(psi_start[[i]], fn=one_elbo_theta0, gr=one_elbo_grad_theta0, beta=beta, Sigma=Sigma, Yi=one_data$y, Zi=one_data$z, Xi=one_data$x, control=list(fnscale=-1), method='L-BFGS-B')
    return(this_opt$par)
  }, mc.cores=nodes)
}

profile_elbo <- function(theta, data, nodes=1, verbose=FALSE) {
  P <- ncol(data[[1]]$z)
  K <- ncol(data[[1]]$x)
  beta <- theta[1:K]
  Sigma_D <- exp(theta[(K+1):(K+P)])
  Sigma_L <- diag(rep(1,P), ncol=P)
  Nsub <- length(data)
  if(P > 1) Sigma_L[lower.tri(Sigma_L, diag=FALSE)] <- theta[-1:(K+P)]
  Sigma <- Sigma_L %*% Sigma_D %*% t(Sigma_L)
  sum(unlist(mclapply(1:length(data), function(i) {
    if(i %% 500 == 0 & verbose) cat("Obs", i, "of", Nsub, "\n")
    one_data <- data[[i]]
    this_opt <- optim(rep(0,2 * P + P * (P-1)/2), fn=one_elbo_theta0, gr=one_elbo_grad_theta0, beta=beta, Sigma=Sigma, Yi=one_data$y, Zi=one_data$z, Xi=one_data$x, control=list(fnscale=-1), method='L-BFGS-B')
    if(this_opt$convergence != 0) cat("WARNING: PSI OPTIMIZATION EXIT CODE NOT 0 AT OBSERVATION ", i, "\n")
    psihat <- this_opt$par
    Li <- diag(rep(1, P), nrow=P)
    if(P > 1) Li[lower.tri(Li, diag=FALSE)] <- psihats[-1:(2*P)]
    one_elbo(beta=beta, Sigma=Sigma, mi=psihat[1:P], Di=psihat[(P+1):(2*P)], Li=Li, Yi=data[[i]]$y, Xi=data[[i]]$x, Zi=data[[i]]$z)
  }, mc.cores=nodes)))
}

profile_elbo_grad <- function(beta, Sigma, data, nodes=1, verbose=FALSE) {
  P <- ncol(data[[1]]$z)
  K <- ncol(data[[1]]$x)
  beta <- theta[1:K]
  Sigma_D <- exp(theta[(K+1):(K+P)])
  Sigma_L <- diag(rep(1,P), ncol=P)
  Nsub <- length(data)
  if(P > 1) Sigma_L[lower.tri(Sigma_L, diag=FALSE)] <- theta[-1:(K+P)]
  Sigma <- Sigma_L %*% Sigma_D %*% t(Sigma_L)
  grads <- mclapply(1:length(data), function(i) {
    if(i %% 500 == 0 & verbose) cat("Obs", i, "of", Nsub, "\n")
    one_data <- data[[i]]
    this_opt <- optim(rep(0,2 * P + P * (P-1)/2), fn=one_elbo_theta0, gr=one_elbo_grad_theta0, beta=beta, Sigma=Sigma, Yi=one_data$y, Zi=one_data$z, Xi=one_data$x, control=list(fnscale=-1), method='L-BFGS-B')
    if(this_opt$convergence != 0) cat("WARNING: PSI OPTIMIZATION EXIT CODE NOT 0 AT OBSERVATION ", i, "\n")
    psihat <- this_opt$par
    Li <- diag(rep(1, P), nrow=P)
    if(P > 1) Li[lower.tri(Li, diag=FALSE)] <- psihats[-1:(2*P)]
    one_elbo_thetgrad(beta=beta, Sigma=Sigma, mi=psihat[1:P], Di=psihat[(P+1):(2*P)], Li=Li, Yi=data[[i]]$y, Xi=data[[i]]$x, Zi=data[[i]]$z)
  }, mc.cores=nodes)
  gradmat <- matrix(unlist(grads), nrow=Nsub, byrow=TRUE)
  colSums(gradmat)
}

elbo_psi_fixed <- function(beta, Sigma, data, psihats, nodes=1, verbose=FALSE) {
  P <- ncol(data[[1]]$z)
  K <- ncol(data[[1]]$x)
  Nsub <- length(data)
  sum(unlist(mclapply(1:length(data), function(i) {
    if(i %% 500 == 0 & verbose) cat("Obs", i, "of", Nsub, "\n")
    one_data <- data[[i]]
    Li <- diag(rep(1, P), nrow=P)
    if(P > 1) Li[lower.tri(Li, diag=FALSE)] <- psihats[[i]][-(1:(2*P))]
    one_elbo(beta=beta, Sigma=Sigma, mi=psihats[[i]][1:P], Di=psihats[[i]][(P+1):(2*P)], Li=Li, Yi=data[[i]]$y, Xi=data[[i]]$x, Zi=data[[i]]$z)
  }, mc.cores=nodes)))
}

elbo_grad_psi_fixed <- function(beta, Sigma, data, psihats, nodes=1, verbose=FALSE) {
  P <- ncol(data[[1]]$z)
  K <- ncol(data[[1]]$x)
  Nsub <- length(data)
  grads <- mclapply(1:length(data), function(i) {
    if(i %% 500 == 0 & verbose) cat("Obs", i, "of", Nsub, "\n")
    one_data <- data[[i]]
    Li <- diag(rep(1, P), nrow=P)
    if(P > 1) Li[lower.tri(Li, diag=FALSE)] <- psihats[[i]][-(1:(2*P))]
    one_elbo_grad_beta(beta=beta, Sigma=Sigma, mi=psihats[[i]][1:P], Di=psihats[[i]][(P+1):(2*P)], Li=Li, Yi=data[[i]]$y, Xi=data[[i]]$x, Zi=data[[i]]$z)
  }, mc.cores=nodes)
  gradmat <- matrix(unlist(grads), nrow=Nsub, byrow=TRUE)
  colSums(gradmat)
}

mixed_logit_vem <- function(data, beta_start=rep(0, ncol(data[[1]]$x)), Sigma_start=diag(rep(1,ncol(data[[1]]$z)), ncol=ncol(data[[1]]$z)), verbose=FALSE, nodes=1, tol=.0001) {
  P <- ncol(data[[1]]$z)
  K <- ncol(data[[1]]$x)
  n <- length(data)
  beta <- beta_start
  if(P == 1) Sigma <- matrix(Sigma_start, nrow=1)
  else Sigma <- Sigma_start
  
  psi_hats <- lapply(1:n, function(i) rep(0, 2 * P + P * (P-1)/2))
  iter <- 1
  while(TRUE) {
    # Get new psihats
    if(verbose) cat("Optimizing psi...\n")
    psi_hats <- optimize_psi(beta=beta, Sigma=Sigma, data=data, psi_start=psi_hats,nodes=nodes)
    if(verbose) cat("Optimizing Sigma...\n")
    Sigma <- matrix(0, nrow=P, ncol=P)
    for(i in 1:n) {
      mi <- psi_hats[[i]][1:P]
      Di <- psi_hats[[i]][(P+1):(2*P)]
      Li <- diag(rep(1,P))
      if(P > 1) Li[lower.tri(Li, diag=FALSE)] <- psi_hats[[i]][(2*P+1):(2*P+P*(P-1)/2)]
      Sigma <- Sigma + mi %*% t(mi) + Li %*% diag(exp(Di), nrow=P) %*% t(Li)
    }
    Sigma <- Sigma / n
    
    if(verbose) cat("Optimizing beta... \n")
    new_opt <- optim(beta,fn=elbo_psi_fixed, gr=elbo_grad_psi_fixed, Sigma=Sigma, data=data, psihats=psi_hats, nodes=nodes, verbose=FALSE, control=list(fnscale=-1, trace=0), method='L-BFGS-B')
    beta <- new_opt$par
    new_elbo <- new_opt$value
    if(iter > 1) {
      delta <- abs(new_elbo - elbo) / abs(elbo)
      if(verbose) cat("Old elbo = ", elbo, "; New elbo = ", new_elbo, "; Relative change = ", delta, "\n")
      if(delta < tol) break
    }
    iter <- iter + 1
    elbo <- new_elbo
  }
  return(list(beta=beta, Sigma=Sigma, psi=psi_hats, elbo=new_elbo))
}


sandwich_cov <- function(beta, Sigma_D, Sigma_L, m, D, L, data, verbose=TRUE, nodes=1) {
  require(numDeriv)
  K <- length(beta)
  P <- ncol(data[[1]]$z)
  theta_ln <- K + P * (P+1) / 2
  psi_ln <- 2*P + P * (P-1)/2
  N <- length(data)
  
  theta <- c(beta, Sigma_D, Sigma_L[lower.tri(Sigma_L, diag=FALSE)])
  
  require(parallel)
  abs <- mclapply(1:N, function(i) {
    if(verbose & ( i %% 100 == 0)) cat(i, " of ", N, " \n")
    
    grad <- grad(one_elbo_theta, theta, mi=m[i,], Di=D[i,], Li=L[[i]], Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)#one_elbo_thetagrad(beta, Sigma, m[i,], D[i,], L[[i]], Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
    hess <- hessian(one_elbo_thetapsi, c(theta, m[i,], D[i,], L[[i]][lower.tri(L[[i]], diag=FALSE)]), Yi=data[[i]]$y, Zi=data[[i]]$z, Xi=data[[i]]$x)
    
    Ai <- hess[1:theta_ln,1:theta_ln] - hess[1:theta_ln, -(1:theta_ln)] %*% solve(hess[(theta_ln + 1):(theta_ln + psi_ln), (theta_ln + 1):(theta_ln + psi_ln)]) %*% hess[-(1:theta_ln), 1:theta_ln]
    
    Bi <- outer(grad, grad)
    c(c(Ai), c(Bi))
  }, mc.cores=nodes)
  abs <- matrix(unlist(abs), nrow=length(abs), byrow=TRUE)
  abhat <- colMeans(abs)
  ahat <- matrix(abhat[1:theta_ln^2], nrow=theta_ln)
  bhat <- matrix(abhat[-(1:theta_ln^2)], nrow=theta_ln)
  return(list(ahat=ahat, bhat=bhat, sand_cov=solve(ahat) %*% bhat %*% solve(ahat) / N))
}

LDL <- function(mat) {
  U <- chol(mat)   ## start from factor given by chol 
  D <- diag(U)   ## extract sqrt(D) 
  L <- t(U/D)     ## get unit lower triangular factor 
  return(list(D=D^2, L=L))
}

score_info <- function(beta, Lambda, Theta, data, verbose=FALSE, nodes=1) {
  # Decompose Sigma = Theta %*% diag(exp(Lambda)) %*% t(Theta) (LDL decomp)
  K <- length(beta)
  P <- length(Lambda)
  theta_ln <- K + P * (P+1) / 2
  if(nodes > 1) {
    require(parallel)
    liks <- unlist(mclapply(1:length(data), function(i) {
      if(verbose & (i %% 100 == 0)) cat("Computing likelihood for obs ", i, " of ", length(data), ' \n')
      one_lik(beta, Lambda, Theta, one_data=data[[i]])
    }, mc.cores=nodes))
    grads <- mclapply(1:length(data), function(i) {
      if(verbose & (i %% 100 == 0)) cat(paste("Computing gradient for obs ", i, " of ", length(data), '\n'))
      one_lik_grad(beta, Lambda, Theta, one_data=data[[i]])
    }, mc.cores=nodes)
    grads <- matrix(unlist(grads), nrow=length(data), byrow=TRUE)
  } else{
    liks <- sapply(1:length(data), function(i) {
      if(verbose & (i %% 100 == 0)) print(paste("Computing likelihood for obs ", i, " of ", length(data)))
      one_lik(beta, Lambda, Theta, one_data=data[[i]])
    })
    grads <- t(sapply(1:length(data), function(i) {
      if(verbose & (i %% 100 == 0)) print(paste("Computing gradient for obs ", i, " of ", length(data)))
      one_lik_grad(beta, Lambda, Theta, one_data=data[[i]])
    }))
  }
  
  score <- colMeans(grads / liks)
  
  outers <- apply(grads, 1, function(row) outer(row, row))
  obs_info <- matrix(colMeans(t(outers) / liks^2), nrow=theta_ln)
  return(list(score=score, obs_info=obs_info))
}

one_lik <- function(beta, Lambda, Theta, one_data) {
  K <- length(beta)
  P <- length(Lambda)
  theta_ln <- K + P * (P+1) / 2
  
  Sigma <- Theta %*% diag(exp(Lambda), nrow=length(Lambda)) %*% t(Theta)
  det_Sigma <- det(Sigma)
  Sigma_inv <- solve(Sigma)
  sqrt_Sigma <- Theta %*% diag(exp(.5 * Lambda), nrow=length(Lambda))
  
  kernel <- function(gamma) {
    gammatilde <- gamma %*% t(sqrt_Sigma)
    mat <- sapply(1:nrow(one_data$x), function(j) {
      dbinom(one_data$y[j], size=1, prob=expit(sum(one_data$x[j,] * beta) + gammatilde %*% one_data$z[j,]))
    })
    apply(mat,1, prod)
  }
  
  return(quadrature(kernel, grid))
}

lik_fn_theta <- function(theta, one_data) {
  K <- ncol(one_data$x)
  P <- ncol(one_data$z)
  beta <- theta[1:K]
  Lambda <- theta[(K+1):(K+P)]
  Theta <- diag(rep(1,P), nrow=P)
  if(P > 1) Theta[lower.tri(Theta, diag=FALSE)] <- theta[-(1:9)]
  one_lik(beta, Lambda, Theta, one_data=one_data)
}

one_lik_grad <- function(beta, Lambda, Theta, one_data) {
  grad(lik_fn_theta, c(beta, Lambda, Theta[lower.tri(Theta, diag=FALSE)]), one_data=one_data)
}


log_lik_fn <- function(theta, data, verbose=FALSE) {
  sum(sapply(1:length(data), function(i) {
    if(verbose & (i %% 100 == 0)) print(paste("Computing likelihood for obs ", i, " of ", length(data)))
    lik_fn_theta(theta, one_data=data[[i]])
  }))
}

log_lik_grad <- function(theta, data, verbose=FALSE) {
  grads <- t(sapply(1:length(data), function(i) {
    if(verbose & (i %% 100 == 0)) print(paste("Computing gradient for obs ", i, " of ", length(data)))
    lik_grad_fn_theta(theta, one_data=data[[i]])
  }))
  colSums(grads)
}

lik_grad_fn_theta <- function(theta, one_data) {
  K <- ncol(one_data$x)
  P <- ncol(one_data$z)
  beta <- theta[1:K]
  Lambda <- theta[(K+1):(K+P)]
  Theta <- diag(rep(1,P), nrow=P)
  if(P > 1) Theta[lower.tri(Theta, diag=FALSE)] <- theta[-(1:9)]
  one_lik_grad(beta, Lambda, Theta, one_data=one_data)
}
