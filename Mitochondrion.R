N <- 10000 # number of max mitochondria in cell 
# Nm <- d_mit * N # number of mitochondria
Nm <- 100 # number of mitochondria in cell
Nmb <- 1000 # number of Mb in cell
S <- floor(N^(1/2))
Nb <- 100
masses <- c(1, 10, 50, 100, 500, 1000, 2000, 3000)


getDist3D <- function(masses, Nm, S,  Nb) {
  R <- 1:(S^3) # 3D space
  
  as <- vector(mode = "list", length = length(masses))
  bs <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    Nmb <- ceiling(x = masses[i]^(1.21))
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    # browser()
    
    mean.dist <- numeric(length = Nb)
    all.dist <- matrix(data = NA, nrow = Nb, ncol = Nmb)
    for(b in 1:Nb) {
      d <- numeric(length = Nmb)
      cell <- array(data = 0, dim = c(S, S, S))
      
      Im <- sample(x = R, size = Nm, replace = FALSE)
      Imb <- replicate(n = 3, expr = runif(n = Nmb, min = 1, max = S))
      if(is.vector(x = Imb) == TRUE) {
        Imb <- matrix(data = Imb, nrow = 1)
      }
      
      cell[Im] <- 1
      i.Im <- which(cell == 1, arr.ind = TRUE)
      i.Im <- i.Im + 0.5
      
      for(j in 1:nrow(Imb)) {
        d[j] <- min(sqrt( (Imb[j, 1]-i.Im[, 1])^2 
                             + (Imb[j, 2]-i.Im[, 2])^2 
                             + (Imb[j, 3]-i.Im[, 3])^2 ))
        all.dist[b, j] <- d[j]
      }
      mean.dist[b] <- mean(d)
    }
    
    # Nm <- ceiling((3/4)*Nm)
    as[[i]] <- mean.dist
    bs[[i]] <- all.dist
  }
  
  return (list(mean.dist = as, 
               all.dist = bs))
}


getDist2D <- function(masses, S,  Nb) {
  R <- 1:(S^2) # 2D space
  
  as <- vector(mode = "list", length = length(masses))
  bs <- vector(mode = "list", length = length(masses))
  for(i in 1:length(masses)) {
    Nmb <- ceiling(x = masses[i]^(1.21))
    Nm <- 100
    # Nm <- ceiling(x = 0.3*masses[i]^(-1/4))
    cat("Mass:", masses[i], "; Nm:", Nm, "; Nmb:", Nmb, "\n", sep = ' ')
    
    # browser()
    mean.dist <- numeric(length = Nb)
    all.dist <- matrix(data = NA, nrow = Nb, ncol = Nm)
    for(b in 1:Nb) {
      d <- numeric(length = Nm)
      cell <- array(data = 0, dim = c(S, S))
      
      Im <- sample(x = R, size = Nm, replace = FALSE)
      Imb <- replicate(n = 2, expr = runif(n = Nmb, min = 1, max = S))
      if(is.vector(x = Imb) == TRUE) {
        Imb <- matrix(data = Imb, nrow = 1)
      }
      
      cell[Im] <- 1
      i.Im <- which(cell == 1, arr.ind = TRUE)
      i.Im <- i.Im + 0.5
      
      
      for(j in 1:nrow(i.Im)) {
        d[j] <- min(sqrt( (Imb[, 1]-i.Im[j, 1])^2 
                          + (Imb[, 2]-i.Im[j, 2])^2 ))
        all.dist[b, j] <- d[j]
      }
      
      mean.dist[b] <- mean(d)
    }
    
    as[[i]] <- mean.dist
    bs[[i]] <- all.dist
  }
  
  return (list(mean.dist = as, 
               all.dist = bs))
}



D3 <- getDist3D(masses = masses, Nm = Nm, S = S, Nb = Nb)
D2 <- getDist2D(masses = masses, S = S, Nb = Nb)


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