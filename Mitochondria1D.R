N <- 10000 # number of max mitochondria in cell 
# Nm <- d_mit * N # number of mitochondria
Nm <- 100 # number of mitochondria in cell
Nmb <- 1000 # number of Mb in cell
Nb <- 100
masses <- c(1, 10, 50, 100, 500, 1000, 2000, 3000)



getDist1D <- function(masses, S,  Nb, N) {
  R <- 1:(S) # 2D space
  
  result.mt <- vector(mode = "list", length = length(masses))
  result.mb <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    # Nmb <- ceiling(x = masses[i]^(1.21))
    # Nm <- 100
    Nm <- ceiling(x = 0.3*masses[i]^(-1/4)*N)
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    
    mean.dist.mt <- numeric(length = Nb)
    mean.dist.mb <- numeric(length = Nb)
    
    for(b in 1:Nb) {
      dist.mt <- numeric(length = Nm)
      dist.mb <- numeric(length = Nmb)
      cell <- numeric(length = S)
      
      # generate mitochondrion positions
      Im <- sample(x = R, size = Nm, replace = FALSE)
      cell[Im] <- 1
      i.Im <- which(cell == 1, arr.ind = TRUE)
      i.Im <- i.Im + 0.5
      
      
      # generate myoglobin positions
      Imb <- runif(n = Nmb, min = 1, max = S)
      
      
      # mitochondrion-specific
      for(j in 1:Nm) {
        dist.mt[j] <- min(sqrt( (Imb-i.Im[j])^2 ))
      }
      mean.dist.mt[b] <- mean(dist.mt)
      
      
      # Mb-specific
      for(j in 1:Nmb) {
        dist.mb[j] <- min(sqrt( (Imb[j]-i.Im)^2 ))
      }
      mean.dist.mb[b] <- mean(dist.mb)
    }
    
    result.mt[[i]] <- mean.dist.mt
    result.mb[[i]] <- mean.dist.mb
  }
  
  return (list(result.mt = result.mt, 
               result.mb = result.mb))
}


getDist2D <- function(masses, S,  Nb, N) {
  R <- 1:S^2 # 2D space
  
  result.mt <- vector(mode = "list", length = length(masses))
  result.mb <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    # Nmb <- ceiling(x = masses[i]^(1.21))
    # Nm <- 100
    Nm <- ceiling(x = 0.3*masses[i]^(-1/4)*N)
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    
    mean.dist.mt <- numeric(length = Nb)
    mean.dist.mb <- numeric(length = Nb)
    
    for(b in 1:Nb) {
      dist.mt <- numeric(length = Nm)
      dist.mb <- numeric(length = Nmb)
      cell <- matrix(data = 0, nrow = S, ncol = S)
      
      # generate mitochondrion positions
      Im <- sample(x = R, size = Nm, replace = FALSE)
      cell[Im] <- 1
      i.Im <- which(cell == 1, arr.ind = TRUE)
      i.Im <- i.Im - 0.5
      
      
      # generate myoglobin positions
      Imb <- replicate(n = 2, expr = runif(n = Nmb, min = 1, max = S))
      
      
      # mitochondrion-specific
      for(j in 1:Nm) {
        dist.mt[j] <- min(sqrt( (Imb[, 1]-i.Im[j, 1])^2 + (Imb[, 2]-i.Im[j, 2])^2 ))
      }
      mean.dist.mt[b] <- mean(dist.mt)
      
      
      # Mb-specific
      for(j in 1:Nmb) {
        dist.mb[j] <-min(sqrt( (Imb[j, 1]-i.Im[, 1])^2 + (Imb[j, 2]-i.Im[, 2])^2 ))
      }
      mean.dist.mb[b] <- mean(dist.mb)
    }
    
    result.mt[[i]] <- mean.dist.mt
    result.mb[[i]] <- mean.dist.mb
  }
  
  return (list(result.mt = result.mt, 
               result.mb = result.mb))
}


getDist3D <- function(masses, S,  Nb, N) {
  R <- 1:S^3 # 2D space
  
  result.mt <- vector(mode = "list", length = length(masses))
  result.mb <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    # Nmb <- ceiling(x = masses[i]^(1.21))
    # Nm <- 100
    Nm <- ceiling(x = 0.3*masses[i]^(-1/4)*N)
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    
    mean.dist.mt <- numeric(length = Nb)
    mean.dist.mb <- numeric(length = Nb)
    
    for(b in 1:Nb) {
      dist.mt <- numeric(length = Nm)
      dist.mb <- numeric(length = Nmb)
      cell <- array(data = 0, dim = c(S, S, S))
      
      # generate mitochondrion positions
      Im <- sample(x = R, size = Nm, replace = FALSE)
      cell[Im] <- 1
      i.Im <- which(cell == 1, arr.ind = TRUE)
      i.Im <- i.Im - 0.5
      
      
      # generate myoglobin positions
      Imb <- replicate(n = 3, expr = runif(n = Nmb, min = 1, max = S))
      
      
      # mitochondrion-specific
      for(j in 1:Nm) {
        dist.mt[j] <- min(sqrt( (Imb[, 1]-i.Im[j, 1])^2 
                                + (Imb[, 2]-i.Im[j, 2])^2 
                                + (Imb[, 3]-i.Im[j, 3])^2))
      }
      mean.dist.mt[b] <- mean(dist.mt)
      
      
      # Mb-specific
      for(j in 1:Nmb) {
        dist.mb[j] <- min(sqrt( (Imb[j, 1]-i.Im[, 1])^2 
                                + (Imb[j, 2]-i.Im[, 2])^2 
                                + (Imb[j, 3]-i.Im[, 3])^2))
      }
      mean.dist.mb[b] <- mean(dist.mb)
    }
    
    result.mt[[i]] <- mean.dist.mt
    result.mb[[i]] <- mean.dist.mb
  }
  
  return (list(result.mt = result.mt, 
               result.mb = result.mb))
}



D1 <- getDist1D(masses = masses, S = N, Nb = Nb, N = N)
D2 <- getDist2D(masses = masses, S = floor(N^(1/2)), Nb = Nb, N = N)
D3 <- getDist3D(masses = masses, S = floor(N^(1/3)), Nb = Nb, N = N)

summary(lm(log10(unlist(lapply(D1$result.mb, mean))) ~ log10(masses)))
summary(lm(log10(unlist(lapply(D2$result.mb, mean))) ~ log10(masses)))
summary(lm(log10(unlist(lapply(D3$result.mb, mean))) ~ log10(masses)))


Mb.conc <- ceiling(x = masses^(1.21))
# Nm.conc <- ceiling(x = 100*masses^(-1/4))


meanD <- numeric(length = length(D2$mean.dist))
for(s in 1:length(D2$mean.dist)) {
  meanD[s] <- mean(D2$mean.dist[[s]])
}

summary(lm(log10(meanD)~log10(masses)))
summary(lm(log10(meanD)~log10(Mb.conc)))
summary(lm(log10(Mb.conc)~log10(masses)))

cols <- c("#000000", "#ff0000", "#ffb3b3", "#99b3ff", "#0039e6", "#00cc33", "darkgray")
plot(density(d$mean.dist[[1]]), xlim = c(0, 6))
for(s in 1:length(d$mean.dist)) {
  if(s > 1) {
    points(density(d$all.dist[[s]]), type = "l", col = cols[s])
  }
}

# require(ggplot2)
# ggplot(data = data.frame(d = meanD, n = d$Nms))+
#   geom_point(aes(x = n, y = d))+
#   scale_y_log10()+
#   scale_x_log10()

# ggplot(data = data.frame(d = meanD, n = d$Nms))+
#   geom_point(aes(x = n, y = d))

masses <- (1:10)^10

summary(lm(log10(d$Nms) ~ log10(meanD)))
summary(lm(log10(d$Nms) ~ log10(masses)))
summary(lm(log10(meanD) ~ log10(d$Nms)))
summary(lm(d$Nms ~ meanD))
summary(lm(meanD ~ d$Nms))


# make a grid Nm x Nmb and compute ratio for different grids