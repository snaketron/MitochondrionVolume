N <- 8000 # number of max mitochondria in cell 
Nmb <- 100 # number of Mb in cell
Nb <- 100
masses <- c(1, 10, 50, 100, 500, 1000, 2000, 3000)

# Assumption:
# Mb density ~ 0.3*Mass^(-1/4)

mitdens.allometry <- function(m, n) {
  return(ceiling(x = 0.3*m^(-1/4)*n))
}


getDist1D <- function(masses, S, Nb, N, Nmb, mitdens.func) {
  R <- 1:(S)
  
  # result.mt <- vector(mode = "list", length = length(masses))
  result.mb <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    Nm <- mitdens.func(m = masses[i], n = N)
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    
    mean.dist.mb <- numeric(length = Nb)
    
    for(b in 1:Nb) {
      # dist.mt <- numeric(length = Nm)
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
      mean.dist.mb[b] <- mean(dist.mb)
    }
    result.mb[[i]] <- mean.dist.mb
  }
  return (list(result.mb = result.mb))
}


getDist2D <- function(masses, S,  Nb, N, Nmb, mitdens.func) {
  R <- 1:(S^2)
  
  result.mb <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    Nm <- mitdens.func(m = masses[i], n = N)
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    
    mean.dist.mb <- numeric(length = Nb)
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
        dist.mb[j] <-min(sqrt( (Imb[j, 1]-i.Im[, 1])^2 + (Imb[j, 2]-i.Im[, 2])^2 ))
      }
      mean.dist.mb[b] <- mean(dist.mb)
    }
    
    result.mb[[i]] <- mean.dist.mb
  }
  
  return (list(result.mb = result.mb))
}


getDist3D <- function(masses, S,  Nb, N, Nmb, mitdens.func) {
  R <- 1:(S^3) # 3D space
  
  result.mb <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    Nm <- mitdens.func(m = masses[i], n = N)
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    
    mean.dist.mb <- numeric(length = Nb)
    for(b in 1:Nb) {
      dist.mb <- numeric(length = Nmb)
      cell <- array(data = 0, dim = c(S, S, S))
      
      # generate mitochondrion positions
      Im <- sample(x = R, size = Nm, replace = FALSE)
      cell[Im] <- 1
      i.Im <- which(cell == 1, arr.ind = TRUE)
      i.Im <- i.Im - 0.5
      
      # generate myoglobin positions
      Imb <- replicate(n = 3, expr = runif(n = Nmb, min = 1, max = S))
      
      # Mb-specific
      for(j in 1:Nmb) {
        dist.mb[j] <- min(sqrt( (Imb[j, 1]-i.Im[, 1])^2 
                                + (Imb[j, 2]-i.Im[, 2])^2 
                                + (Imb[j, 3]-i.Im[, 3])^2))
      }
      mean.dist.mb[b] <- mean(dist.mb)
    }
    result.mb[[i]] <- mean.dist.mb
  }
  return (list(result.mb = result.mb))
}



D1 <- getDist1D(masses = masses, S = N, Nb = Nb, N = N, 
                Nmb = Nmb, mitdens.func = mitdens.allometry)
D2 <- getDist2D(masses = masses, S = floor(N^(1/2)), 
                Nb = Nb, N = floor(N^(1/2))^2, 
                Nmb = Nmb, mitdens.func = mitdens.allometry)
D3 <- getDist3D(masses = masses, S = floor(N^(1/3)), 
                Nb = Nb, N = floor(N^(1/3))^3, 
                Nmb = Nmb, mitdens.func = mitdens.allometry)

summary(lm(log10(unlist(lapply(D1$result.mb, mean))) ~ log10(masses)))
summary(lm(log10(unlist(lapply(D2$result.mb, mean))) ~ log10(masses)))
summary(lm(log10(unlist(lapply(D3$result.mb, mean))) ~ log10(masses)))

summary <- c()
for(i in 1:length(masses)) {
  temp <- rbind(data.frame(mean.dist = mean(D1$result.mb[[i]]),
                           dimension = "D1", mass = masses[i]),
                data.frame(mean.dist = mean(D2$result.mb[[i]]),
                           dimension = "D2", mass = masses[i]),
                data.frame(mean.dist = mean(D3$result.mb[[i]]),
                           dimension = "D3", mass = masses[i]))
  summary <- rbind(summary, temp)
}
rm(i, temp)






require(ggplot2)


ggplot(data = summary)+
  facet_wrap(facets = ~dimension, nrow = 1)+
  geom_point(aes(x = mass, y = mean.dist))+
  theme_bw()



ggplot(data = summary)+
  facet_wrap(facets = ~dimension, nrow = 1)+
  geom_point(aes(x = mass, y = mean.dist))+
  theme_bw()


ggplot(data = summary)+
  geom_point(aes(x = mass, y = mean.dist, col = dimension))+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  theme(legend.position = "top")+
  annotation_logticks(base = 10, sides = "bl")

