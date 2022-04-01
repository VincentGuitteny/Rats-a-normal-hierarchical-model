# Data
ni = 30 # nombres de lignes (rats)
nj = 5 # nombres de colonnes (jours)

y = structure(c(151, 145, 147, 155, 135, 159, 141, 159, 177, 134, 160, 143, 154, 171, 163, 160, 
                142, 156, 157, 152, 154, 139, 146, 157, 132, 160, 169, 157, 137, 153, 199, 199, 
                214, 200, 188, 210, 189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 
                212, 203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, 263, 237, 
                230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, 248, 234, 243, 259, 246, 
                253, 225, 229, 250, 237, 257, 261, 248, 219, 244, 283, 293, 312, 272, 280, 298, 
                275, 297, 350, 260, 313, 273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 
                272, 285, 286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 338, 
                376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 334, 302, 302, 323, 
                331, 345, 333, 316, 291, 324), 
              .Dim = c(ni, nj))

colnames(y) = c("Jour 8", "Jour 15", "Jour 22", "Jour 29", "Jour 36")

x = c(8.0, 15.0, 22.0, 29.0, 36.0)
xbar = 22.0

Gibbs = function(nchain, data){
  
  ni = nrow(data)
  nj = ncol(data)
  
  x = c(8.0, 15.0, 22.0, 29.0, 36.0)
  xbar = 22.0
  
  # État initial suggéré par les auteurs
  alpha = c(250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 
            250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250)
  beta = c(6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
           6, 6)
  alpha.c = 150
  beta.c = 10
  sigma.c = 1
  alpha.sigma = 1
  beta.sigma = 1
  
  alpha0 = alpha.c - xbar * beta.c
  
  init = c(alpha0, alpha.c, beta.c, alpha.sigma, beta.sigma, sigma.c, alpha, beta)
  
  # Début de la chaîne
  chain = matrix(NA, nchain + 1, 66)
  chain[1,] = init
  for (k in 1:nchain){
    # Mise à jour de alpha.c
    mean = (1000**2 * sum(alpha))/(alpha.sigma**2 + ni * 1000**2)
    sd = (1000**2 * alpha.sigma**2)/(alpha.sigma**2 + ni * 1000**2)
    
    alpha.c = rnorm(1, mean, sd)
    
    # Mise à jour de beta.c
    mean = (1000**2 * sum(beta))/(beta.sigma**2 + ni * 1000**2)
    sd = (1000**2 * beta.sigma**2)/(beta.sigma**2 + ni * 1000**2)
    
    beta.c = rnorm(1, mean, sd)
    
    # Mise à jour de alpha.sigma
    a = 0.001 + ni / 2
    b = (2 * 0.001 + sum((alpha - alpha.c)**2))/2
    
    alpha.sigma = 1/sqrt(rgamma(1, a, 1/b))
    
    # Mise à jour de beta.sigma
    a = 0.001 + ni / 2
    b = (2 * 0.001 + sum((beta - beta.c)**2))/2
      
    beta.sigma = 1/sqrt(rgamma(1, a, 1/b))
      
    # Mise à jour de sigma.c
    a = 0.001 + ni * nj / 2
    tot = 0
    for (i in 1:ni){
      for (j in 1:nj){
        s = (data[i,j] - alpha[i] - beta[i]*(x[j] - xbar))**2
        tot = tot + s
      }
    }
    b = (2 * 0.001 + tot)/2
    
    sigma.c = 1/sqrt(rgamma(1, a, 1/b))
    
    # Mise à jour de alpha
    for (i in 1:ni){
      mean = (alpha.c * sigma.c**2 + alpha.sigma**2 * sum(data[i,] - beta[i]*(x - xbar))) / 
        (sigma.c**2 + nj * alpha.sigma**2)
      sd = (alpha.sigma**2 * sigma.c**2) / (sigma.c**2 + nj * alpha.sigma**2)
        
      alpha[i] = rnorm(1, mean, sd)
    }
    
    # Mise à jour de beta
    for (i in 1:ni){
      mean = (beta.c * sigma.c**2 + beta.sigma**2 * sum((x - xbar) * (data[i,] - alpha[i]))) / 
        (sigma.c**2 + sum((x - xbar)**2) * beta.sigma**2)
      sd = (beta.sigma**2 * sigma.c**2) / (sigma.c**2 + sum((x - xbar)**2) * beta.sigma**2)
      
      beta[i] = rnorm(1, mean, sd)
    }
    
    # Mise à jour de alpha0
    alpha0 = alpha.c - xbar * beta.c
    
    # Mise à jour de la chaîne
    chain[k+1,] = c(alpha0, alpha.c, beta.c, alpha.sigma, beta.sigma, sigma.c, alpha, beta)
  }
  return(chain)
}

# 
chain = Gibbs(10**4, y)

library(coda)
burnin = 1:1000
plot(mcmc(chain[-burnin,])[,c(1,3,6)])

s = summary(mcmc(chain[-burnin,]))$statistics[c(1,3,6),]

y2 = matrix(NA, ni, nj)
for(i in 1:ni){
  for (j in 1:nj){
    y2[i,j] = rnorm(1, chain[10001,i+6] + chain[10001,i+36] * (x[j] - xbar), chain[10001,6]**2)
  }
}
