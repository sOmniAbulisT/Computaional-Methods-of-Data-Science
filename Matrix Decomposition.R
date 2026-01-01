rm(list = ls())

set.seed(20250903)
library(Matrix)
test_mat <- matrix(rnorm(5*5), 5, 5)
test_mat <- test_mat + diag(1e-8, 5)
LU(test_mat)

A <- matrix(c(3, -1, 2, 1, 2, 3, 2, -2, 1), nrow = 3, byrow = T)
LU(A)

n <- nrow(A)
L <- diag(1, n)
U <- matrix(0, n , n)

###----- first row of U -----### 
U[1, ] <- A[1, ]

###----- first column of l -----###
# l21
L[2, 1] <- A[2, 1] / U[1 ,1]

# l31
L[3, 1] <- A[3, 1] / U[1, 1]

for(i in 2:3){
  L[i, 1] <- A[i, 1] / U[1, 1]
}

###----- second row of U -----###
# u22
U[2, 2] <- A[2, 2] - L[2, 1] * U[1, 2]

# u23
U[2, 3] <- A[2, 3] - L[2, 1] * U[1, 3]

for(j in 2:3){
  U[2, j] <- A[2, j] - L[2, 1] * U[1, j]
}

###----- second column of L -----###
# l32
L[3, 2] <- (A[3, 2] - L[3, 1] * U[1, 2]) / U[2, 2]

###----- third row of U -----###
# u33
U[3, 3] <- A[3, 3] - L[3, 1] * U[1, 3] - L[3, 2] * U[2, 3]

###----- QR Decomposition -----###
A <- matrix(c(-1, -1, 1, 1, 3, 3, -1, -1, 5, 1, 3, 7), nrow = 4, ncol = 3, byrow = TRUE)

## Gram-Schmidt Algorithm
Q <- matrix(0, nrow = 4, ncol = 3)
R <- matrix(0, nrow = 3, ncol = 3)

# First column of Q and R
q1_tilde <- A[, 1]
R[1, 1] <- q1_tilde^2 |> sum() |> sqrt()
Q[, 1] <- q1_tilde / R[1, 1]

# Second column of Q and R
R[1, 2] <- Q[, 1] |> t() %*% A[, 2]
q2_tilde <- A[, 2] - R[1, 2] %*% Q[, 1]
R[2, 2] <- q2_tilde^2 |> sum() |> sqrt()
Q[, 2] <- q2_tilde / R[2, 2]

# Third column of Q and R
R[1, 3] <- Q[, 1] |> t() %*% A[, 3]
R[2, 3] <- Q[, 2] |> t() %*% A[, 3]
q3_tilde <- A[, 3] - R[1, 3] %*% Q[, 1] - R[2, 3] %*% Q[, 2]
R[3, 3] <- q3_tilde^2 |> sum() |> sqrt()
Q[, 3] <- q3_tilde / R[3, 3]
all.equal(Q%*%R, A)

## Modified Gram-Schmidt Algorithm
A <- matrix(c(-1, -1, 1, 1, 3, 3, -1, -1, 5, 1, 3, 7), nrow = 4, ncol = 3, byrow = TRUE)
B <- A
Q <- matrix(0, nrow = 4, ncol = 3); R <- matrix(0, nrow = 3, ncol = 3)

# K = 1 (first row of R and first column of Q)
R[1, 1] <- B[, 1]^2 |> sum() |> sqrt()
Q[, 1] <- B[, 1] / R[1, 1]

R[1, 2] <- Q[, 1] %*% B[, 2] |> sum()
B[, 2] <- B[, 2] - R[1, 2] %*% Q[, 1]

R[1, 3] <- Q[, 1] %*% B[, 3] |> sum()
B[, 3] <- B[, 3] - R[1, 3] %*% Q[, 1]

# K = 2 (second row of R and second column of Q)
R[2, 2] <- B[, 2]^2 |> sum() |> sqrt()
Q[, 2] <- B[, 2] / R[2, 2]

R[2, 3] <- Q[, 2] %*% B[, 3] |> sum()
B[, 3] <- B[, 3] - R[2, 3] %*% Q[, 2]

# K = 3 (third row of R and third column of Q)
R[3, 3] <- B[, 3]^2 |> sum() |> sqrt()
Q[, 3] <- B[, 3] / R[3, 3]

a <- QR(A)

### power iteration ###
A <- matrix(c(4,1,0,1,3,1,0,1,2), 3, 3, byrow = 1)
PowerIter(A)

#### SVD ####
A <- matrix(c(1,0,0,0,2,0,0,3,0,0,0,0,0,0,0,0,2,0,0,0), nrow = 4, ncol = 5, byrow = TRUE)

### Step 1
m <- nrow(X); n <- ncol(X)
M <- t(X)%*%X
d <- PowerIter(M)

### Step 2
S <- matrix(0, m, n)

for(i in seq_len(m)){
  S[i, i] <- d$EigenValue[i]
}

V <- d$vectors

### Step 3
U <- matrix(0, m, m)
for(j in seq_len(m)){
  U[, j] <- (X%*%V[, j])/S[j, j]
}

U%*%S%*%t(V)
SVD(X)
### ICA ###
## Centering
A <- matrix(c(1,3,1,2,2,3,0,3,5,4,4,5,5,5,3,4), nrow = 8, ncol = 2, byrow = T)
A <- t(A)
mu <- rowMeans(A)

D <- A - mu

## Decorreltion
m <- ncol(D)
COV <- (D %*% t(D))/(m - 1)

eig <- PowerIter(COV)
lambda <- eig$EigenValue
vecs <- eig$EigenVector

D <- t(D)
U <-  D%*%vecs

## Whitening
Z <- diag(1/sqrt(lambda)) %*% t(U)
Z <- round(-Z, 2)
