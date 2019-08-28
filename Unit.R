A <- replicate(n = 1, expr = runif(n = 10^4, min = 0, max = 1))
B <- replicate(n = 1, expr = runif(n = 10^4, min = 0, max = 1))
d1 <- (sqrt((A[, 1]-B[, 1])^2))
mean(d1)
hist(d1) # 0.33

A <- replicate(n = 2, expr = runif(n = 10^4, min = 0, max = 1))
B <- replicate(n = 2, expr = runif(n = 10^4, min = 0, max = 1))
d2 <- (sqrt((A[, 1]-B[, 1])^2 + (A[, 2]-B[, 2])^2))
mean(d2)
hist(d2) # 0.52

A <- replicate(n = 3, expr = runif(n = 10^4, min = 0, max = 1))
B <- replicate(n = 3, expr = runif(n = 10^4, min = 0, max = 1))
d3 <- (sqrt((A[, 1]-B[, 1])^2 + (A[, 2]-B[, 2])^2 + (A[, 3]-B[, 3])^2))
hist(d3)
mean(d3)  # 0.66
