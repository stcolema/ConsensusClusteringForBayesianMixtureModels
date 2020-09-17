
library(magrittr)
library(mdiHelpR)
library(pheatmap)

set.seed(1)

R <- 200000
N <- 9
K <- 3
N_k <- N / K

cl <- rep(c(1:K), each = N_k)



# cm <- matrix(0, nrow = N, ncol = N)
# 
# for( r in 1:R ){
#   c_star <- 1:N
#   i <- sample(N, size = 1)
#   c_i <- cl[i]
#   J <- which(cl == c_i)
#   J <- J[J != i]
#   j <- sample(J, size = 1)
#   c_star[c(i, j)] <- c_star[i]
#   cc <- createSimilarityMat(matrix(c_star, nrow = 1))
#   cm <- cm + cc
# }
# cm <- cm / R
# 
# pheatmap::pheatmap(cm)
# 
# cm <- matrix(0, nrow = N, ncol = N)
# for( r in 1:R ){
#   c_star <- sample(1:K, size = N, replace = T)
#   i <- sample(N, 1)
#   c_i <- cl[i]
#   J <- which(cl == c_i)
#   J <- J[J != i]
#   j <- sample(J, 1)
#   c_star[c(i, j)] <- c_star[i]
#   cc <- createSimilarityMat(matrix(c_star, nrow = 1))
#   cm <- cm + cc
# }
# cm <- cm / R
# pheatmap::pheatmap(cm)

K_max <- min(round(N/2), 50)

cm <- matrix(0, nrow = N, ncol = N)
w <- matrix(0, nrow = R, ncol = K_max)

for( r in 1:R ){
  set.seed(r)
  # v0 <- rbeta(K_max, 1, 1)
  # stick <- 1
  # w0 <- rep(0, K_max)
  # 
  # for(i in 1:K_max){
  #   w0[i] <- v0[i] * stick
  #   stick <- stick - w0[i]
  # }
  
  c_star <- sample(1:K_max, size = N, replace = T) # , prob = w0)
  i <- sample(N, 1)
  c_i <- cl[i]
  J <- which(cl == c_i)
  J <- J[J != i]
  j <- sample(J, 1)
  c_star[c(i, j)] <- c_star[i]
  cc <- createSimilarityMat(matrix(c_star, nrow = 1))
  cm <- cm + cc
  # w[r, ] = w0
}
cm <- cm / R

pheatmap::pheatmap(cm)

cm[1, 2]
cm[1, N]
cm[1, 2] - cm[1, N]
2/(N*(N_k - 1))
(1/K_max) + (1 - 1/K_max)*(2/(N*(N_k - 1)))
