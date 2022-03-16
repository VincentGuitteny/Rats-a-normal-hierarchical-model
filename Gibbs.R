# Data
ni = 30 # nombres de lignes (rats)
nj = 5 # nombres de colonnes (jours)

y = structure(c(151, 145, 147, 155, 135, 159, 141, 159, 177, 134, 
                 160, 143, 154, 171, 163, 160, 142, 156, 157, 152, 154, 139, 146, 
                 157, 132, 160, 169, 157, 137, 153, 199, 199, 214, 200, 188, 210, 
                 189, 201, 236, 182, 208, 188, 200, 221, 216, 207, 187, 203, 212, 
                 203, 205, 190, 191, 211, 185, 207, 216, 205, 180, 200, 246, 249, 
                 263, 237, 230, 252, 231, 248, 285, 220, 261, 220, 244, 270, 242, 
                 248, 234, 243, 259, 246, 253, 225, 229, 250, 237, 257, 261, 248, 
                 219, 244, 283, 293, 312, 272, 280, 298, 275, 297, 350, 260, 313, 
                 273, 289, 326, 281, 288, 280, 283, 307, 286, 298, 267, 272, 285, 
                 286, 303, 295, 289, 258, 286, 320, 354, 328, 297, 323, 331, 305, 
                 338, 376, 296, 352, 314, 325, 358, 312, 324, 316, 317, 336, 321, 
                 334, 302, 302, 323, 331, 345, 333, 316, 291, 324), 
              .Dim = c(L, C))

colnames(y) = c("Jour 8", "Jour 15", "Jour 22", "Jour 29", "Jour 36")

x = c(8.0, 15.0, 22.0, 29.0, 36.0)
xbar = 22.0

Gibbs = function(nchain, data, prop_sd){
  
  x = c(8.0, 15.0, 22.0, 29.0, 36.0)
  xbar = 22.0
  
  # État initial
  alpha.c = data[1,1]
  beta.c = 0
  alpha.sigma = sd(data[,1])
  beta.sigma = 1
  sigma.c = sd(data)
  
  alpha = rep(mean(data[,1]), ni)
  beta = rep(mean(data[,2])/mean(data[,1]), ni)
  
  init = c(alpha.c, beta.c, alpha.sigma, beta.sigma, sigma.c, alpha, beta)
  
  # Début de la chaîne
  chain = matrix(NA, nchain +1, 7)
  chain[1,] = init
  for (i in 1:nchain){
    # Mise à jour de alpha.c
    prop = rnorm(1, alpha.c, prop_sd[1])
    
    mean = (prop_sd[1] * sum(alpha))/(alpha.sigma**2 + ni * prop_sd[1])
    sd = (prop_sd[1] * alpha.sigma**2)/(alpha.sigma**2 + ni * prop_sd[1])
    
    top = dnorm(prop, mean, sd)
    
    bottom = dnorm(alpha.c, mean, sd)
    
    acc_prob = exp(top - bottom)
    
    if (runif(1) < acc_prob){
      alpha.c = prop
    }
    
    # Mise à jour de beta.c
    prop = rnorm(1, beta.c, prop_sd[2])
    
    mean = (prop_sd[2] * sum(alpha))/(beta.sigma**2 + ni * prop_sd[2])
    sd = (prop_sd[2] * beta.sigma**2)/(beta.sigma**2 + ni * prop_sd[2])
    
    top = dnorm(prop, mean, sd)
    
    bottom = dnorm(beta.c, mean, sd)
    
    acc_prob = exp(top - bottom)
    
    if (runif(1) < acc_prob){
      beta.c = prop
    }
    
    # Mise à jour de alpha.sigma
    
    # Mise à jour de beta.sigma
    
    # Mise à jour de sigma.c
    
    # Mise à jour de alpha
    
    # Mise à jour de beta
    
    
    # Mise à jour de la chaîne
    chain[i+1,] = c(alpha.c, beta.c, alpha.sigma, beta.sigma, sigma.c, alpha, beta)
  }
}