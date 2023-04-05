library(BB)


#=====================#
# estimation of V
#=====================#

get.V <- function(n, rho, beta) {
  ind <- t(combn(1:n, 2))
  w <- rho/as.vector(1+exp(rho%*%beta))
  V <- colMeans(w[1:(n-1),]) %*% t(colMeans(w[1:(n-1),]))
  for (i in 2:n){
    set <- which(rowSums(ind==i)>0)
    V <- V + colMeans(w[set,]) %*% t(colMeans(w[set,]))
  }
  V <- V/(n-1)
  return(V)
}


#=====================#
# estimation of Q
#=====================#

get.Q <- function(rho, beta){
  m <- nrow(rho)
  exp.rho <- as.vector(exp(rho%*%beta))
  w <- sqrt(exp.rho)/(1+exp.rho) * rho
  Q <- t(w)%*%w/(2*m)
  return(Q)
}


#===============================#
# shift-adjusted estimation 
# under dataset shift (SAE1)
#===============================#


SAE1 <- function(y,x,x1,y1, rep=100){
  
  n1 <- nrow(x1) 
  out <- lm(y1~x1)
  beta.old <- out$coefficients[-1]
  sigma.old <- vcov(out)[-1,-1]*n1
  
  
  n <- nrow(x)
  ind <- combn(1:n, 2)
  y.diff <- y[ind[1,]] - y[ind[2,]]
  x.diff <- x[ind[1,],] - x[ind[2,],]
  rho <- y.diff * x.diff
  
  V.old <- get.V(n, rho, beta.old)
  
  sbeta <- function(beta){
    s <- rho/as.vector(1+exp(rho%*%beta)) 
    return (colMeans(s))
  }

  
  outb <- BBsolve(beta.old, sbeta, quiet=T, control = list(maxit=rep))
  
  beta.new <- outb$par
  Q <- get.Q(rho, beta.new) 
  sigma.new <- solve(Q) %*% V.old %*% solve(Q)/n
  
  list(beta.new=beta.new, sigma.new=sigma.new)
}


