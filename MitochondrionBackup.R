Nm <- 2000 # number of mitochondria (20 % of cell)
Nmb <- 500 # number of Mb
N <- 10000 # number of spaces in cell 
S <- floor(N^(1/3))
Nb <- 200
Ns <- 10


getD <- function(Nm, S,  Nb, Ns) {
  R <- 1:(S^3) # 3D space
  
  Nms <- numeric(length = Ns)
  a <- vector(mode = "list", length = Ns)
  for(s in 1:Ns) {
    cat("Ns:", s, "; Nm:", Nm, "\n", sep = ' ')
    
    d <- matrix(data = 0, nrow = Nb, ncol = Nmb)
    for(i in 1:Nb) {
      cell <- array(data = 0, dim = c(S, S, S))
      
      Im <- sample(x = R, size = Nm, replace = FALSE)
      Imb <- sample(x = R[!R %in% Im], size = Nmb, replace = FALSE)
      
      cell[Im] <- 1
      cell[Imb] <- -1
      
      i.Im <- which(cell == 1, arr.ind = TRUE)
      i.Imb <- which(cell == -1, arr.ind = TRUE)
      
      # browser()
      for(j in 1:nrow(i.Imb)) {
        d[i, j] <- min(sqrt( (i.Imb[j, 1]-i.Im[, 1])^2 
                             + (i.Imb[j, 2]-i.Im[, 2])^2 
                             + (i.Imb[j, 3]-i.Im[, 3])^2 ))
      }
    }
    
    Nms[s] <- Nm
    Nm <- ceiling((3/4)*Nm)
    a[[s]] <- d
  }
  
  return (list(a = a, 
               Nms = Nms))
}


d <- getD(Nm = Nm, S = S, Nb = Nb, Ns = Ns)


plot(density(rowMeans(d$a[[1]])), xlim = c(1, 10))
meanD <- numeric(length = 10)
for(s in 1:10) {
  meanD[s] <- mean(rowMeans(d$a[[s]]))
  if(s > 1) {
    points(density(rowMeans(d$a[[s]])), type = "l")
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