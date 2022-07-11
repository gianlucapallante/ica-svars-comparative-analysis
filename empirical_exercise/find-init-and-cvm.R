#source("fuctions-needed/testlik.R")

find_init_cvm <- function(res, n_lhs, ordering,C){
  ## n_lhs number of initial conditions
  ## res: reduced form residuals kxn
  ## C: Choleski decomposition of reduced form residuals variance-cov matrix
  u_chol <- res
  tu_chol <- t(res)
  set.seed(55) ## To make example reproducible
  lhs    <- 2*pi*improvedLHS(n = n_lhs, k = 3)
  
  
  
  ## Sample size
  #n <- nrow(resid(var.interactions))
  n <- ncol(res)
  
  pseudo.log.lik3D.super        <- function(c) {
    c11 <- c[1]
    c21 <- c[2]
    c31 <- c[3]
    c12 <- c[4]
    c22 <- c[5]
    c32 <- c[6]
    c13 <- c[7]
    c23 <- c[8]
    c33 <- c[9]
    
    c.1y <- c11*u_chol[,1] + c21*u_chol[,2] + c31*u_chol[,3]
    c.2y <- c12*u_chol[,1] + c22*u_chol[,2] + c32*u_chol[,3]
    c.3y <- c13*u_chol[,1] + c23*u_chol[,2] + c33*u_chol[,3]
    
    
    g <- rep(NA,1)
    g[1] <- -1*(n*log(1/sqrt(2*pi))-sum((c.1y)^2/2) + 
                  n*log(gamma((5+1)/2)/(sqrt(5*pi)*gamma(5/2)))-((5+1)/2)*sum(log(1+(c.2y)^2/5)) + 
                  n*log(gamma((12+1)/2)/(sqrt(12*pi)*gamma(12/2)))-((12+1)/2)*sum(log(1+(c.3y)^2/12)))
    g
  }
  
  ## Sub gaussian
  pseudo.log.lik3D.sub        <- function(c) {
    c11 <- c[1]
    c21 <- c[2]
    c31 <- c[3]
    c12 <- c[4]
    c22 <- c[5]
    c32 <- c[6]
    c13 <- c[7]
    c23 <- c[8]
    c33 <- c[9]
    
    c.1y <- c11*u_chol[,1] + c21*u_chol[,2] + c31*u_chol[,3]
    c.2y <- c12*u_chol[,1] + c22*u_chol[,2] + c32*u_chol[,3]
    c.3y <- c13*u_chol[,1] + c23*u_chol[,2] + c33*u_chol[,3]
    
    
    g <- rep(NA,1)
    g[1] <- -1*(n*log(1/sqrt(2*pi)) + sum(pi*(c.1y)^2+log(cosh((pi/2)*c.1y))) + 
                  n*log(1/sqrt(2*pi)) + sum(pi*(c.2y)^2+log(cosh((pi/2)*c.2y))) + 
                  n*log(1/sqrt(2*pi)) + sum(pi*(c.3y)^2+log(cosh((pi/2)*c.3y))))
    g
  }
  ## Constraint
  orth.constr3D    <- function(c){
    c11 <- c[1]
    c21 <- c[2]
    c31 <- c[3]
    c12 <- c[4]
    c22 <- c[5]
    c32 <- c[6]
    c13 <- c[7]
    c23 <- c[8]
    c33 <- c[9]
    g <- rep(NA,6)
    g[1] <- c11*c12 + c21*c22 + c31*c32 # c.1*c.2
    g[2] <- c11*c13 + c21*c23 + c31*c33 # c.1*c.3
    g[3] <- c12*c13 + c22*c23 + c32*c33 # c.2*c.3
    g[4] <- c11^2 + c21^2 + c31^2 - 1         # c.1*c.1
    g[5] <- c12^2 + c22^2 + c32^2 - 1         # c.2*c.2
    g[6] <- c13^2 + c23^2 + c33^2 - 1         # c.3*c.3
    g
  }
  
  Wscaled.dc <- vector("list", n_lhs)
  Wscaled.ica <- vector("list", n_lhs)
  w0          <- vector("list", n_lhs)
  A.ica       <- vector("list", n_lhs)
  A.dc        <- vector("list", n_lhs)
  A.ID.pml    <- vector("list", n_lhs)
  A.ID.dc     <- vector("list", n_lhs)
  A.ID.ica    <- vector("list", n_lhs)
  A.ID.pml    <- vector("list", n_lhs)
  
  ts_dim <- (1 + sqrt(1 + 8 * ncol(lhs))) / 2
  combns <- combn(ts_dim, 2, simplify = FALSE)
  
  for (kk in 1:nrow(lhs)){
    if (kk %% 100 == 0) {
      print(paste0("Initialization # ",kk))
    }
    #### Initialization ####
    we <- diag(ts_dim)
    pv <- lhs[kk,]
    for (i in seq_along(pv)) {
      
      tmp <- diag(ts_dim)
      tmp[combns[[i]][1], combns[[i]][1]] <- cos(pv[i])
      tmp[combns[[i]][2], combns[[i]][2]] <- cos(pv[i])
      tmp[combns[[i]][1], combns[[i]][2]] <- - sin(pv[i])
      tmp[combns[[i]][2], combns[[i]][1]] <- sin(pv[i])
      we <- we %*% tmp
    }
    w0[[kk]] <- we
    
    ## PML
    p.start           <- as.vector(w0[[kk]])
    
    please           <- auglag(par = p.start, fn = pseudo.log.lik3D.super, heq = orth.constr3D,
                               control.outer = list(trace=F))$par 
    Apml.raw         <- matrix(data = please,nrow = 3,byrow = F)
    Bpml             <- (C) %*% (Apml.raw)
    if (ordering == "maxfinder") {
      A.ID.pml[[kk]]   <- maxfinder(A = Bpml)$A.id
    }else if(ordering == "lanne-saik"){
      A.ID.pml[[kk]]  <- lanne_saik(original_M = Bpml)
    }
    
    Apml             <- solve(C) %*% A.ID.pml[[kk]]
    
    ## FastICA
    icares             <- fastICA(tu_chol, n.comp = ncol(sigg),tol=1e-14, maxit=3000, verbose=FALSE, w.init = w0[[kk]])
    W                  <- t((icares$K) %*% (icares$W)) 
    Wscal.ica          <- rescaleVar(W_hat = W,ut = u_chol)$Ws 
    Wscaled.ica[[kk]]  <- Wscal.ica %*% solve(C) 
    A.ica[[kk]]        <- solve(Wscaled.ica[[kk]]) # A is the mixing matrix
    
    ## Distance Covariance ##
    DC                 <- steadyICA(X = tu_chol,n.comp = ncol(sigg),w.init = w0[[kk]])
    Wscaled.dc[[kk]]   <- solve(DC$W) %*% solve(C)
    A.dc[[kk]]         <- solve(Wscaled.dc[[kk]])
    
    #Identification scheme
    if (ordering == "maxfinder") {
      A.ID.dc[[kk]]       <- maxfinder(A = A.dc[[kk]])$A.id
      A.ID.ica[[kk]]      <- maxfinder(A = A.ica[[kk]])$A.id
    }else if(ordering == "lanne-saik"){
      A.ID.dc[[kk]]       <- lanne_saik(original_M = A.dc[[kk]])
      A.ID.ica[[kk]]      <- lanne_saik(original_M = A.ica[[kk]]) 
    }
    
  }
  dd <- NULL
  ##CVM seems not to require any initialization
  ## the file dd contains the copula under the null hypothesis of statistical independence
  if (is.null(dd)) {
    dd <- indepTestSim(n, ncol(tu_chol), verbose = T)
  }
  #save(dd, file = 'indepNULL3D_n400.RData')
  lower <- rep(0, ncol(tu_chol) * (ncol(tu_chol) - 1)/2)
  upper <- rep(pi, ncol(tu_chol) * (ncol(tu_chol) - 1)/2)
  de_control <- list(itermax = 500, steptol = 100, 
                     trace = FALSE)
  ## First step of optimization with DEoptim
  de_res <- DEoptim(testlik, lower = lower, upper = upper, 
                    control = de_control, faklow = C, u = u.ica, dd = dd)
  
  
  k <- ncol(tu_chol)
  iter2  <- 150
  ## Second step of optimization. Creating randomized starting angles around the optimized angles
  ## here maybe you insert the LHS stuff
  theta_rot <- matrix(rnorm(n = (k * (k - 1)/2) * iter2, mean = de_res$optim$bestmem, 
                            sd = 0.3), (k * (k - 1)/2), iter2)
  theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
  # Start vectors for iterative optimization approach
  startvec_list <- as.list(as.data.frame(theta_rot))
  erg_list <- lapply(X = startvec_list, FUN = optim, fn = testlik, 
                     gr = NULL, faklow = C, u = u.ica, dd = dd, 
                     method = ifelse(k == 2, "Brent", "Nelder-Mead"), 
                     lower = ifelse(k == 2, -.Machine$double.xmax/2, -Inf), 
                     upper = ifelse(k == 2, .Machine$double.xmax/2, Inf), 
                     control = list(maxit = 1000), 
                     hessian = FALSE)
  # Print log-likelihood values from local maxima
  logliks <- sapply(erg_list, "[[", "value")
  if (min(logliks) < de_res$optim$bestval) {
    params <- sapply(erg_list, "[[", "par", simplify = FALSE)
    par_o <- params[[which.min(logliks)]]
    logs <- min(logliks)
    inc <- 1
  } else {
    par_o <- de_res$optim$bestmem
    logs <- de_res$optim$bestval
    inc <- 0
  }
  Acvm <- rotmat(par_o, C)   ## Estimated mixing matrix
  
  #c.ica <- vector("list", n_lhs)
  #c.dc  <- vector("list", n_lhs)
  #c.pml <- vector("list",n_lhs)
  
  
  results <- list(w0 =  w0,
                  Acvm  = Acvm, A.ID.dc = A.ID.dc, A.ID.ica = A.ID.ica,
                  A.ID.pml = A.ID.pml)
  
  return(results)
}
