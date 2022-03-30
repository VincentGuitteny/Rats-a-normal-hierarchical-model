N = 5
T = 30
xbar = 22
 
  tau.c = dgamma(0.001,0.001)
  sigma <- 1 / sqrt(tau.c)
  alpha.c = dnorm(0.0,1.0E-6)
  alpha.tau = dgamma(0.001,0.001)
  beta.c = dnorm(0.0,1.0E-6)
  beta.tau = dgamma(0.001,0.001)
  alpha0 <- alpha.c - xbar * beta.c
  
  alpha = rep(1:N)
  beta = rep(1:N)
  Y = matrix(1:T*N, ncol=T, nrow=N, byrow=FALSE)
  mu = matrix(1:T*N, ncol=N, nrow=T, byrow=FALSE)
  
  for( i in 1 : N ) {
    alpha[i] = dnorm(alpha.c,alpha.tau)
    beta[i] = dnorm(beta.c,beta.tau)
    datatest = data[i,]
    for( j in 1 : T ) {
      mu[i , j] <- alpha[i] + beta[i] * (datatest[j] - xbar)
      Y[i , j] = dnorm(mu[i , j],tau.c)
      
    }
    
  }
  

alpha
mu






Y = matrix(1:T, ncol=2, nrow=5, byrow=FALSE)
Y[3,2]
mu = matrix(1:T*N, ncol=T, nrow=T, byrow=FALSE)
mu
