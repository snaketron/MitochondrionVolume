require(parallel)
require(ggplot2)

S <- 30 # side dimension
Nmb <- 1000 # number of Mb in cell
Nb <- 100
masses <- c(1, 10, 100, 1000, 10000, 100000)

# Assumption:
# Mb density ~ 0.3*Mass^(-1/4)

mitdens.allometry <- function(m, n) {
  return(ceiling(x = 0.2*m^(-1/4)*n))
}


getMinDist1D <- function(mass, S, Nb, N, Nmb, mitdens.func) {
  R <- 1:(S)
  
  Nm <- mitdens.func(m = mass, n = N)
  
  mean.min.dist.mb <- numeric(length = Nb)
  for(b in 1:Nb) {
    dist.mb <- numeric(length = Nmb)
    cell <- numeric(length = S)
    
    # generate mitochondrion positions
    Im <- sample(x = R, size = Nm, replace = FALSE)
    cell[Im] <- 1
    i.Im <- which(cell == 1, arr.ind = TRUE)
    i.Im <- i.Im + 0.5
    
    # generate myoglobin positions
    Imb <- runif(n = Nmb, min = 1, max = S)
    
    # Mb-specific
    for(j in 1:Nmb) {
      dist.mb[j] <- min(sqrt( (Imb[j]-i.Im)^2 ))
    }
    mean.min.dist.mb[b] <- mean(dist.mb)
  }
  
  return(mean.min.dist.mb)
}

getMeanDist1D <- function(mass, S, Nb, N, Nmb, mitdens.func) {
  R <- 1:(S)
  
  Nm <- mitdens.func(m = mass, n = N)
  
  mean.min.dist.mb <- numeric(length = Nb)
  for(b in 1:Nb) {
    dist.mb <- numeric(length = Nmb)
    cell <- numeric(length = S)
    
    # generate mitochondrion positions
    Im <- sample(x = R, size = Nm, replace = FALSE)
    cell[Im] <- 1
    i.Im <- which(cell == 1, arr.ind = TRUE)
    i.Im <- i.Im + 0.5
    
    # generate myoglobin positions
    Imb <- runif(n = Nmb, min = 1, max = S)
    
    # Mb-specific
    for(j in 1:Nmb) {
      dist.mb[j] <- mean(sqrt( (Imb[j]-i.Im)^2 ))
    }
    mean.min.dist.mb[b] <- mean(dist.mb)
  }
  
  return(mean.min.dist.mb)
}


getMinDist2D <- function(mass, S,  Nb, N, Nmb, mitdens.func) {
  R <- 1:(S^2)
  
  Nm <- mitdens.func(m = mass, n = N)
  cat("Mass:", mass, "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
  
  mean.min.dist.mb <- numeric(length = Nb)
  for(b in 1:Nb) {
    dist.mb <- numeric(length = Nmb)
    
    cell <- matrix(data = 0, nrow = S, ncol = S)
    
    # generate mitochondrion positions
    Im <- sample(x = R, size = Nm, replace = FALSE)
    cell[Im] <- 1
    i.Im <- which(cell == 1, arr.ind = TRUE)
    i.Im <- i.Im - 0.5
    
    
    # generate myoglobin positions
    Imb <- replicate(n = 2, expr = runif(n = Nmb, min = 1, max = S))
    
    # Mb-specific
    for(j in 1:Nmb) {
      dist.mb[j] <- min(sqrt( (Imb[j, 1]-i.Im[, 1])^2 + (Imb[j, 2]-i.Im[, 2])^2 ))
    }
    mean.min.dist.mb[b] <- mean(dist.mb)
  }
  
  return (mean.min.dist.mb)
}

getMeanDist2D <- function(mass, S,  Nb, N, Nmb, mitdens.func) {
  R <- 1:(S^2)
  
  Nm <- mitdens.func(m = mass, n = N)
  cat("Mass:", mass, "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
  
  mean.min.dist.mb <- numeric(length = Nb)
  for(b in 1:Nb) {
    dist.mb <- numeric(length = Nmb)
    
    cell <- matrix(data = 0, nrow = S, ncol = S)
    
    # generate mitochondrion positions
    Im <- sample(x = R, size = Nm, replace = FALSE)
    cell[Im] <- 1
    i.Im <- which(cell == 1, arr.ind = TRUE)
    i.Im <- i.Im - 0.5
    
    
    # generate myoglobin positions
    Imb <- replicate(n = 2, expr = runif(n = Nmb, min = 1, max = S))
    
    # Mb-specific
    for(j in 1:Nmb) {
      dist.mb[j] <- mean(sqrt( (Imb[j, 1]-i.Im[, 1])^2 + (Imb[j, 2]-i.Im[, 2])^2 ))
    }
    mean.min.dist.mb[b] <- mean(dist.mb)
  }
  
  return (mean.min.dist.mb)
}



getMinDist3D <- function(mass, S, Nb, N, Nmb, mitdens.func) {
  R <- 1:(S^3)
  
  Nm <- mitdens.func(m = mass, n = N)
  cat("Mass:", mass, "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
  
  mean.min.dist.mb <- numeric(length = Nb)
  for(b in 1:Nb) {
    dist.mb <- numeric(length = Nmb)
    cell <- array(data = 0, dim = c(S, S, S))
    
    # generate mitochondrion positions
    Im <- sample(x = R, size = Nm, replace = FALSE)
    cell[Im] <- 1
    i.Im <- which(cell == 1, arr.ind = TRUE)
    i.Im <- i.Im - 0.5
    
    # generate myoglobin positions
    Imb <- replicate(n = 3, expr = runif(n = Nmb, min = 0, max = S))
    
    # Mb-specific
    for(j in 1:Nmb) {
      dist.mb[j] <- min(sqrt( (Imb[j, 1]-i.Im[, 1])^2 
                              + (Imb[j, 2]-i.Im[, 2])^2 
                              + (Imb[j, 3]-i.Im[, 3])^2))
    }
    
    # average min distance
    mean.min.dist.mb[b] <- mean(dist.mb)
  }
  return (mean.min.dist.mb)
}

getMeanDist3D <- function(mass, S, Nb, N, Nmb, mitdens.func) {
  R <- 1:(S^3)
  
  Nm <- mitdens.func(m = mass, n = N)
  cat("Mass:", mass, "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
  
  mean.min.dist.mb <- numeric(length = Nb)
  for(b in 1:Nb) {
    dist.mb <- numeric(length = Nmb)
    cell <- array(data = 0, dim = c(S, S, S))
    
    # generate mitochondrion positions
    Im <- sample(x = R, size = Nm, replace = FALSE)
    cell[Im] <- 1
    i.Im <- which(cell == 1, arr.ind = TRUE)
    i.Im <- i.Im - 0.5
    
    # generate myoglobin positions
    Imb <- replicate(n = 3, expr = runif(n = Nmb, min = 0, max = S))
    
    # Mb-specific
    for(j in 1:Nmb) {
      dist.mb[j] <- mean(sqrt( (Imb[j, 1]-i.Im[, 1])^2 
                                        + (Imb[j, 2]-i.Im[, 2])^2 
                                        + (Imb[j, 3]-i.Im[, 3])^2))
    }
    
    # average min distance
    mean.min.dist.mb[b] <- mean(dist.mb)
  }
  return (mean.min.dist.mb)
}


minD1 <- mclapply(X = masses, FUN = getMinDist1D, S = S, Nb = Nb, N = S, 
                  Nmb = Nmb, mitdens.func = mitdens.allometry, mc.cores = 4)
minD2 <- mclapply(X = masses, FUN = getMinDist2D, S = S, Nb = Nb, N = S^2, 
                  Nmb = Nmb, mitdens.func = mitdens.allometry, mc.cores = 4)
minD3 <- mclapply(X = masses, FUN = getMinDist3D, S = S, Nb = Nb, N = S^3, 
                  Nmb = Nmb, mitdens.func = mitdens.allometry, mc.cores = 4)

meanD1 <- mclapply(X = masses, FUN = getMeanDist1D, S = S, Nb = Nb, N = S, 
                   Nmb = Nmb, mitdens.func = mitdens.allometry, mc.cores = 4)
meanD2 <- mclapply(X = masses, FUN = getMeanDist2D, S = S, Nb = Nb, N = S^2, 
                   Nmb = Nmb, mitdens.func = mitdens.allometry, mc.cores = 4)
meanD3 <- mclapply(X = masses, FUN = getMeanDist3D, S = S, Nb = Nb, N = S^3, 
                   Nmb = Nmb, mitdens.func = mitdens.allometry, mc.cores = 4)

minD1 <- unlist(lapply(X = minD1, FUN = mean))
minD2 <- unlist(lapply(X = minD2, FUN = mean))
minD3 <- unlist(lapply(X = minD3, FUN = mean))

meanD1 <- unlist(lapply(X = meanD1, FUN = mean))
meanD2 <- unlist(lapply(X = meanD2, FUN = mean))
meanD3 <- unlist(lapply(X = meanD3, FUN = mean))




summary <- rbind(data.frame(mean.min.dist = minD1, mean.dist = meanD1, dimension = "D1", mass = masses),
                 data.frame(mean.min.dist = minD2, mean.dist = meanD2, dimension = "D2", mass = masses),
                 data.frame(mean.min.dist = minD3, mean.dist = meanD3,  dimension = "D3", mass = masses))


ggplot(data = summary)+
  geom_point(aes(x = mass, y = mean.min.dist, col = dimension))+
  theme_bw()+
  scale_x_log10()+
  theme(legend.position = "top")+
  annotation_logticks(base = 10, sides = "b")


ggplot(data = summary)+
  geom_point(aes(x = mass, y = mean.dist, col = dimension))+
  theme_bw()+
  scale_x_log10()+
  theme(legend.position = "top")+
  annotation_logticks(base = 10, sides = "b")


# coefficients
summary(lm(log10(minD1) ~ log10(masses)))
summary(lm(log10(minD2) ~ log10(masses)))
summary(lm(log10(minD3) ~ log10(masses)))

