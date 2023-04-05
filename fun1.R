library(BB)


#=====================#
# estimation of V
#=====================#

get.V <- function(nb, rho, beta) {
  ind <- t(combn(1:nb, 2))
  w <- rho/as.vector(1+exp(rho%*%beta))
  V <- colMeans(w[1:(nb-1),]) %*% t(colMeans(w[1:(nb-1),]))
  for (i in 2:nb){
    set <- which(rowSums(ind==i)>0)
    V <- V + colMeans(w[set,]) %*% t(colMeans(w[set,]))
  }
  V <- V/(nb-1)
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
# under dataset shift 
#===============================#


fun1 <- function(yb, xb, nb, beta.old, sigma.old, N, rep=100){
    
    ind <- combn(1:nb, 2)
    y.diff <- yb[ind[1,]] - yb[ind[2,]]
    x.diff <- xb[ind[1,],] - xb[ind[2,],]
    rho <- y.diff * x.diff
    
    sbeta <- function(beta){
      s <- rho/as.vector(1+exp(rho%*%beta)) 
      return ( colMeans(s))
    }
    
    # outb <- BBsolve(beta.old, sbeta, quiet=T, control = list(maxit=rep))
    outb <- BBsolve(beta.old, sbeta, quiet=T)
    beta1 <- outb$par
    Q <- get.Q(rho, beta1) 
    V <- get.V(nb, rho, beta1)
    s1 <- Q %*% solve(V) %*% Q*nb/(N+nb)
    s2 <- solve(sigma.old)*N/(N+nb)
    sigma.new <- solve(s1+s2)
    beta.new <- as.vector(solve(s1+s2)%*%(s1%*%beta1 + s2%*%beta.old))
  
    list(beta.new=beta.new, sigma.new=sigma.new)
}


