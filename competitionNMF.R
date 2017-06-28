#####
## OPTIONAL TASK: Competition
#Submitted by
#Harshita Jhavar
#Matriculation Number - 2566267
#####
require(Matrix)
#
# Use this stub if you want to take part to the competition
# The following commands should work in a fresh R instance
#source('competitionNMF.R')
#A <- matrix(runif(10000), ncol=100)
#res <- nmf(A, 12)
#norm(A - res$W %*% res$H, "F")

nmf <- function(A, rank)
{
  ## YOUR CODE HERE
  W <- matrix(runif(nrow(A)*rank), nrow(A), rank) # Initialize W and H
  H <- matrix(runif(rank*ncol(A)), rank, ncol(A))
  K <- nnzero(A) #Computing number of non zero entries in A
  m <- nrow(W) 
  n <- ncol(H)
  ep <- 2.220e-14
  rho.H <- 1 + ((K + m*rank)/(n*rank + n))
  rho.W <- 1 + ((K + n*rank)/(m*rank + m))
  alpha <- 1
  max.iter <- 1+alpha*rho.W
  H<-H*((t(W)%*%A)/((t(W)%*%W%*%H) + ep))
  W<-W*((A%*%t(H))/((W%*%H%*%t(H)) + ep))
  W.previous <- W
  H.previous <- H
  epsilon <- 0.5
  default.maxiter <- 500
  #Main Loop
  for  (iter in 1:  default.maxiter){ 
    for (i in 1: max.iter){ 
      if(i == 1){
        temp <- epsilon * norm(W - W.previous)}
      if(norm(W-W.previous,'F') <= temp) break #Stopping Criterion
    }
  W.previous = W
  H<-H*((t(W)%*%A)/((t(W)%*%W%*%H) + ep))
  W<-W*((A%*%t(H))/((W%*%H%*%t(H)) + ep))
}
    updated.result <- list(W=W, H=H)
    return(updated.result)
}

