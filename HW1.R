######################################################################
####----- This programming assignment is Homework 1 for CMDS -----####
######################################################################

###----- library packages -----###
library(readr)
library(dplyr)
library(tidyr)
library(stringr)

###----- Data loading -----###
lines <- read_lines("CAmaxTemp.txt")

dat <- tibble(raw = lines) |>
  extract(
    raw,
    into = c("ID", "Location", "YearRange", "values"),
    regex = "^(\\d{5})([A-Z ,_]+?)\\s{2,}(\\d{6}-\\d{6})\\s+(.*)$"
  ) |>
  mutate(
    ID = as.integer(ID),
    Location = str_squish(Location),
    values = str_squish(values)       
  ) |>
  separate_wider_delim(
    cols = values, delim = " ",
    names = c(paste0("M", 1:12), "Max")
  ) |>
  mutate(across(M1:Max, as.integer))

X <- as.matrix(dat[1:12, 4:15])

###----- Source required function -----###
source("LU.R")
source("QR.R")
source("eigen.R")

###----- Question 1 -----###
## X: the data matrix
### (a)
n <- nrow(X); L <- diag(1, n); U <- matrix(0, n, n)

# first row of U
U[1, ] <- X[1, ]

# first column of L
for(i in 2:12){
  L[i, 1] <- X[i, 1] / U[1, 1]
}

# second row of U
for(j in 2:12){
  U[2, j] <- X[2, j] - L[2, 1] * U[1, j]
}

# second column of L
for(i in 3:12){
  L[i, 2] <- (X[i, 2] - L[i, 1] * U[1, 2]) / U[2, 2]
}

# third row of U
for(j in 3:12){
  U[3, j] <- X[3, j] - L[3, 1]*U[1, j] - L[3, 2]*U[2, j]
}

# third column of L
for(i in 4:12){
  L[i, 3] <- (X[i, 3] - L[i, 1]*U[1, 3] - L[i, 2]*U[2, 3]) / U[3, 3]
}

result <- LU(X)
U <- result$UpperTriangular; L <- result$LowerTriangular
all.equal(L%*%U, X)
L%*%U
X
LU(X)$L %*% LU(X)$U
### (b) 
m <- nrow(X); n <- ncol(X)
Q <- matrix(0, m, n); R <- matrix(0, n, n); B <- X

# First row of R and first column of Q
for(j in seq_len(n)){
  if(j == 1){
    R[1, j] <- B[, j]^2 |> sum() |> sqrt()
    Q[, j] <- B[, j] / R[j, j]
  } else {
    R[1, j] <- Q[, 1] %*% B[, j] |> sum()
    B[, j] <- B[, j] - R[1, j] %*% Q[, 1]
  }
}

# Second row of R and second column of Q
for(j in 2:n){
  if(j == 2){
    R[2, j] <- B[, j]^2 |> sum() |> sqrt()
    Q[, j] <- B[, j] / R[j, j]
  } else {
    R[2, j] <- Q[, 2] %*% B[, j] |> sum()
    B[, j] <- B[, j] - R[2, j] %*% Q[, 2]
  }
}

# Third row of R and third column of Q
for(j in 3:n){
  if(j == 3){
    R[3, j] <- B[, j]^2 |> sum() |> sqrt()
    Q[, j] <- B[, j] / R[j, j]
  } else {
    R[3, j] <- Q[, 3] %*% B[, j] |> sum()
    B[, j] <- B[, j] - R[3, j] %*% Q[, 3]
  }
}

QR(X)

# Inverse of X
n <- nrow(QR(X)[[1]]); m <- ncol(t(QR(X)[[2]])); x <- matrix(0, n, m)
for(i in n:1){
  RHS <- B[i, ]
  if(i < n) RHS - R[i, (i+1):n] %*% x[(i+1):n, ]
  x[i, ] <- RHS / R[i, i]
}

### (c)
x0 <- rep(1, nrow(X))

# step 1
y1 <- X %*% x0; lam1 <- max(abs(y1)); x1 <- y1 / lam1

# step 2
y2 <- X %*% x1; lam2 <- max(abs(y2)); x2 <- y2 / lam2

# step 3
y3 <- X %*% x2; lam3 <- max(abs(y3)); x3 <- y3 / lam3

for(i in seq_len(1000)){
  y  <- X %*% x0
  lambda <- max(abs(y))
  x  <- y / lambda
  lambda
}

### (d)
n <- nrow(X); Lambda <- diag(n)
S1 <- QR(X); Q1 <- S1$Q; R1 <- S1$R
A1 <- R1 %*% Q1
Lambda <- Lambda %*% Q1

S2 <- QR(A1); Q2 <- S2$Q; R2 <- S2$R
A2 <- R2 %*% Q2
Lambda <- Lambda %*% Q2

S3 <- QR(A2); Q3 <- S3$Q; R3 <- S3$R
A3 <- R3 %*% Q3
Lambda <- Lambda %*% Q3

for(k in seq_len(2000)){
  S <- QR(X)
  Qk <- S$Q; Rk <- S$R
  Ak <- Rk %*% Qk
  Lambda <- Lambda %*% Qk
  
  off <- sqrt(sum(Ak[lower.tri(Ak)]^2))
  if (off < 1e-12) break 
}
OthoIter(X)
###----- Question 2 -----###
### (a)
n <- nrow(X); p <- ncol(X)
mu <- colSums(X) / n
Xc <- sweep(X, MARGIN = 2, STATS = mu, FUN = "-")
COV <- (t(Xc) %*% Xc)/(n - 1)

### (b)
B <- COV
n <- nrow(B)
V <- matrix(0, n, 3)      
Lambda <- numeric(3)

for(i in seq_len(3)){
  x0 <- rnorm(12)
  if (i > 1) {
    for (h in 1:(i-1)) x0 <- x0 - sum(x0 * V[,h]) * V[,h]   
  }
  x0 <- x0 / max(abs(x0))
  
  iter <- 0
  while(TRUE){
    iter <- iter + 1
    y <- B %*% x0
    normy <- max(abs(y))
    xk <- as.vector(y / normy)
    
    if (max(abs(xk - x0)) <= 1e-12 || iter >= 2000) { x0 <- xk; break }
    x0 <- xk
  }
  
  v <- x0 / sqrt(sum(x0^2))
  lambda <- as.numeric(t(v) %*% B %*% v)
  
  V[, i] <- v
  Lambda[i] <- lambda
  
  B <- B - lambda * (v %*% t(v))
}
cumulative <- sum(Lambda) / sum(diag(COV))
V
cumulative

### (c)
score <- Xc %*% V
scatterplot3d(score[, 1], score[, 2], score[, 3], 
              xlab = "PC1", ylab = "PC2", zlab = "PC3", 
              pch = 16)

### (d)
SVD(X)

### (f)
lowrank <- SVD(X)
U3 <- lowrank$U[, seq_len(3)]
S3 <- lowrank$S[seq_len(3), seq_len(3)]
V3 <- lowrank$V[, seq_len(3)]
X3 <- U3 %*% S3 %*% t(V3)

###----- Question 3 -----###
X

Y <- cbind(X[,2], X[,6], X[,10])
## Pre-processing 1 Centering
mu <- rowMeans(t(Y))
D <- Y - mu

## Pre-processing 2 Decorrelation
m <- nrow(D)
COV <- (t(D) %*% D)/(m - 1)
eig <- PowerIter(COV)
lambda <- diag(eig$EigenValue)
V <- eig$EigenVector

U <-  D%*%V

## Pro-processing 3 Whitening Scaling
Z <- D %*% V %*% sqrt(lambda)
Z

## (d)
ICA(Z)
