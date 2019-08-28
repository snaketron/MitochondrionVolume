# Given are two points on a unit sphere uniformly at random. 
# What is the expected distance between them?
n <- 3

2^(n-1)/sqrt(pi)*gamma(n/2)^2/gamma(n-1/2) 

