## Data Mining and Matrices
## Summer semester 2017
##
## Assignment 2
## Due 18 June 2017 at 23:59

## Student's information
## Name: Harshita Jhavar
## Matriculation #:2566267

##
## This file contains stubs of code to help you to do your 
## assignment. You can fill your parts at the indicated positions
## and return this file as a part of your solution. Remember that
## if you define any functions etc. in other files, you must include 
## those files, as well. You do not have to include standard libraries
## or any files that we provide.
##
## Remember to fill your name and matriculation number above.
##

## Preamble
###########

## utils.R contains some helper functions
## This might require
#install.packages("Matrix")
source("utils.R")

## Package I used for calculating pseudo inverse
##install.packages("corpcor")
library(corpcor)
## Task 1
##########

## utils.R provides a boilerplate code for Euclidean NMF
## The function prototype is 
## nmf.euclid <- function(A, k, func, nreps=1, maxiter=300, minerr=2.220e-14)
## where A is data, k is rank, func is a pointer to a function that gets A, W, and H and returns 
## updated W and H, nreps is the number of repetitions, maxiter is the maximum number of iterative
## updates and minerr is the smallest error after which we consider we have converged.
##
## Here we define two helper functions that actually implement the updates.

## Boilerplate for ALS updates
nmf.als <- function(A, W, H) {
  ## YOUR UPDATE CODE HERE
  H<-pseudoinverse(W)%*%A
  H[H < 0]<-0 #Truncating negative values to zero
  W<-A%*%pseudoinverse(H)
  W[W < 0]<-0
  ## END YOUR UPDATE CODE
  
  ## Output the updated result
  updated.result <- list(W=W, H=H)
  return(updated.result)
}

## Boilerplate for multiplicative updates
nmf.mult <- function(A, W, H) {
  ## YOUR UPDATE CODE HERE
  ep<-2.220e-14
  H<-H*((t(W)%*%A)/((t(W)%*%W%*%H) + ep))
  W<-W*((A%*%t(H))/((W%*%H%*%t(H)) + ep))
  ## END YOUR UPDATE CODE
  
  ## Output the updated result
  updated.result <- list(W=W, H=H)
  return(updated.result)
}

## Boilerplate for gradient descent
nmf.grad <- function(A, W, H) {
  ## YOUR UPDATE CODE HERE
  neta <- diag(1 / rowSums(t(W)%*%W))
  G<-t(W)%*%W%*%H - t(W)%*%A
  H<- H - (neta%*%G)
  H[H < 0]<-0 #Truncating negative values to zero
  ## END YOUR UPDATE CODE
  
  ## Output the updated result
  updated.result <- list(W=W, H=H)
  return(updated.result)
}

## Load the news data
news <- as.matrix(read.csv("news.csv"))

## Sample use of nmf.euclid for ALS updates:
#res <- nmf.euclid(news, 10, nmf.als) # Computes NMF with ALS updates


## Compare the three algorithms.
## Notice that nmf.euclid returns the best error and the errors per iteration for the run that yield 
## the best error, e.g. if "res" is as above,
#plot(res, type="l")
## shows how the error updates

## YOUR PART STARTS HERE
## plot for ALS
ALS.error = nmf.euclid(news, k=20,maxiter = 30, nreps = 1, nmf.als)
plot(unlist(ALS.error$err.per.iter),type = "l", main="Error Plot for Alternating least squares", sub=NULL, xlab="#iterations", ylab="Error value", col = "red")
print(ALS.error$err) #Best Error 

## plot for multiplicative updates
mult.err = nmf.euclid(news, k=20, maxiter = 30, nreps = 1, nmf.mult)
plot(unlist(mult.err$err.per.iter),type="l", main="Lee&Sung's Multiplicative NMF Algorithm", sub=NULL, xlab="#iterations", ylab="Error value", col = "blue")
print(mult.err$err)

## plot for NMF via gradient descent using Oblique Projected Landweber updates
grad.err = nmf.euclid(news, k=20, maxiter =30, nreps = 1, nmf.grad)
plot(unlist(grad.err$err.per.iter),type="l", main="NMF via gradient descent using Oblique Projected Landweber updates", sub=NULL, xlab="#iterations", ylab="Error value", col = "green")
print(grad.err$err)

## YOUR PART ENDS HERE

#############
## Task 2
#############

## Remember to normalize the data before applying nmf.euclid.
## If res <- nmf.euclid(news, 20, ...), you can print the top-10 terms by
#for (k in 1:nrow(res$H)) print(sort(res$H[k,], decreasing=TRUE)[1:10])

## YOUR PART STARTS HERE
#Normalizing
news<-news/sum(news)
#Finding the NMF for k=20 using ALS method
ALS.error <- nmf.euclid(news, k=20, nmf.als)
#Printing the top 10 terms
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#Finding the NMF for k=5 using ALS method
ALS.error <- nmf.euclid(news, k=5, nmf.als)
#Printing the top 10 terms
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#Finding the NMF for k=14 using ALS method
ALS.error <- nmf.euclid(news, k=14, nmf.als)
#Printing the top 10 terms
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#Finding the NMF for k=32 using ALS method
ALS.error <- nmf.euclid(news, k=32, nmf.als)
#Printing the top 10 terms
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#Finding the NMF for k=40using ALS method
ALS.error <- nmf.euclid(news, k=40, nmf.als)
#Printing the top 10 terms
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#K-L divergence
#For k=20
ALS.error <- nmf.gkl(news, 20)
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#For k=5
ALS.error <- nmf.gkl(news, 5)
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#For k=14
ALS.error <- nmf.gkl(news, 14)
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#For k=32
ALS.error <- nmf.gkl(news, 32)
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])
#For k=40
ALS.error <- nmf.gkl(news, 40)
for (k in 1:nrow(ALS.error$H)) print(sort(ALS.error$H[k,], decreasing=TRUE)[1:10])

## YOUR PART ENDS HERE

## You can compute the KL-divergence minimizing NMF using nmf.gkl found in utils.R
#res <- nmf.gkl(news, 20)


#############
## Task 3
#############

## Remember to use the normalized data!

## Compute the KL-divergence minimizing NMF for k=20
k <- 20
nmf.res <- nmf.gkl(news, k)

# Normalize to W'SH'
nmf.wsh.res <- nmf.wsh(nmf.res)
## and to SW'H'
nmf.swh.res <- nmf.swh(nmf.res)

## Do some clustering
clust <- kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster

## How good is this? The smaller the better
nmi.news(clust)
# 0.8705277
## Now do the the other clusterings and other paramter choices.

## YOUR PART STARTS HERE
#Plotting the clustering by kmeans done above for k=20
plot(clust, main="kmeans Clustering for k=20", col=clust)

#k-means on the first k principal components (Karhunen L´oeve transform) for k=20
SVD<-svd(news)
V = SVD$v
KLT<-news%*%V[,1:20]
clust<- kmeans(KLT, 20, nstart=20,iter.max=15)$cluster
nmi.news(clust)
#0.8794345
plot(clust, main="kmeans on the first 20 principal components (Karhunen L´oeve transform) for k=20", col=clust)

#k-means on the W matrix of the NMF (using KL divergence),
clust<-kmeans(nmf.res$W, 20, nstart=20)$cluster
nmi.news(clust)
#0.8604751
plot(clust, main="k-means on the W matrix of the NMF (using KL divergence) for k=20", col=clust)

#k-means on the W' matrix of factorization W'ΣH' obtained from the NMF, and
clust<-kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
#0.8559922
plot(clust, main="k-means on the W' matrix of factorization W'ΣH' obtained from the NMF for k=20", col=clust)

#k-means on the W' matrix of factorization ΣW'H' obtained from the NMF.
clust<-kmeans(nmf.swh.res$W, 20, nstart=20, iter.max = 15)$cluster
nmi.news(clust)
#0.7969381
plot(clust, main="k-means on the W' matrix of factorization ΣW'H' obtained from the NMF for k=20", col=clust)
###################################################################################################### K=5
#Trying different ranks of the matrix factorization
#k=5
k <- 5
nmf.res <- nmf.gkl(news, k)

# Normalize to W'SH'
nmf.wsh.res <- nmf.wsh(nmf.res)
## and to SW'H'
nmf.swh.res <- nmf.swh(nmf.res)

##kmeans for k=5
clust <- kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
# 0.8663
plot(clust, main="kmeans Clustering for k=5", col=clust)

#k-means on the first k principal components (Karhunen L´oeve transform) for k=5
KLT<-news%*%V[,1:5]
clust<- kmeans(KLT, 20, nstart=20,iter.max=15)$cluster
nmi.news(clust)
#0.8641288
plot(clust, main="kmeans on the first 5 principal components (Karhunen L´oeve transform) for k=5", col=clust)

#k-means on the W matrix of the NMF (using KL divergence),
clust<-kmeans(nmf.res$W, 20, nstart=20)$cluster
nmi.news(clust)
# 0.9104881
plot(clust, main="k-means on the W matrix of the NMF (using KL divergence) for k=5", col=clust)

#k-means on the W' matrix of factorization W'ΣH' obtained from the NMF, and
clust<-kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
# 0.8952437
plot(clust, main="k-means on the W' matrix of factorization W'ΣH' obtained from the NMF for k=5", col=clust)

#k-means on the W' matrix of factorization ΣW'H' obtained from the NMF.
clust<-kmeans(nmf.swh.res$W, 20, nstart=20, iter.max = 15)$cluster
nmi.news(clust)
#0.8081977
plot(clust, main="k-means on the W' matrix of factorization ΣW'H' obtained from the NMF for k=5", col=clust)
###################################################################################### K = 14
#k=14
k <- 14
nmf.res <- nmf.gkl(news, k)

# Normalize to W'SH'
nmf.wsh.res <- nmf.wsh(nmf.res)
## and to SW'H'
nmf.swh.res <- nmf.swh(nmf.res)

##kmeans for k=14
clust <- kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
#0.836479
plot(clust, main="kmeans Clustering for k=14", col=clust)

#k-means on the first k principal components (Karhunen L´oeve transform) for k=14
KLT<-news%*%V[,1:14]
clust<- kmeans(KLT, 20, nstart=20,iter.max=15)$cluster
nmi.news(clust)
#0.8775555
plot(clust, main="kmeans on the first 14 principal components (Karhunen L´oeve transform) for k=14", col=clust)

#k-means on the W matrix of the NMF (using KL divergence),
clust<-kmeans(nmf.res$W, 20, nstart=20)$cluster
nmi.news(clust)
# 0.7964861
plot(clust, main="k-means on the W matrix of the NMF (using KL divergence) for k=14", col=clust)

#k-means on the W' matrix of factorization W'ΣH' obtained from the NMF, and
clust<-kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
#0.8195208
plot(clust, main="k-means on the W' matrix of factorization W'ΣH' obtained from the NMF for k=14", col=clust)

#k-means on the W' matrix of factorization ΣW'H' obtained from the NMF.
clust<-kmeans(nmf.swh.res$W, 20, nstart=20, iter.max = 15)$cluster
nmi.news(clust)
# 0.7697037
plot(clust, main="k-means on the W' matrix of factorization ΣW'H' obtained from the NMF for k=14", col=clust)
################################################################################################################ K=32
#k=32
k <-32
nmf.res <- nmf.gkl(news, k)

# Normalize to W'SH'
nmf.wsh.res <- nmf.wsh(nmf.res)
## and to SW'H'
nmf.swh.res <- nmf.swh(nmf.res)
clust <- kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
#0.8767417
plot(clust, main="kmeans Clustering for k=32", col=clust)

#k-means on the first k principal components (Karhunen L´oeve transform) for k=32
KLT<-news%*%V[,1:32]
clust<- kmeans(KLT, 20, nstart=20,iter.max=15)$cluster
nmi.news(clust)
#0.8673437
plot(clust, main="kmeans on the first 32 principal components (Karhunen L´oeve transform) for k=32", col=clust)

#k-means on the W matrix of the NMF (using KL divergence),
clust<-kmeans(nmf.res$W, 20, nstart=20)$cluster
nmi.news(clust)
#0.8578099
plot(clust, main="k-means on the W matrix of the NMF (using KL divergence) for k=32", col=clust)

#k-means on the W' matrix of factorization W'ΣH' obtained from the NMF, and
clust<-kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
#0.9159513
plot(clust, main="k-means on the W' matrix of factorization W'ΣH' obtained from the NMF for k=32", col=clust)

#k-means on the W' matrix of factorization ΣW'H' obtained from the NMF.
clust<-kmeans(nmf.swh.res$W, 20, nstart=20, iter.max = 15)$cluster
nmi.news(clust)
#0.8276406
plot(clust, main="k-means on the W' matrix of factorization ΣW'H' obtained from the NMF for k=32", col=clust)
#k=40
k <-40
nmf.res <- nmf.gkl(news, k)

# Normalize to W'SH'
nmf.wsh.res <- nmf.wsh(nmf.res)
## and to SW'H'
nmf.swh.res <- nmf.swh(nmf.res)
clust <- kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
# 0.9717254
plot(clust, main="kmeans Clustering for k=40", col=clust)

#k-means on the first k principal components (Karhunen L´oeve transform) for k=40
KLT<-news%*%V[,1:40]
clust<- kmeans(KLT, 20, nstart=20,iter.max=15)$cluster
nmi.news(clust)
# 0.8553645
plot(clust, main="kmeans on the first 40 principal components (Karhunen L´oeve transform) for k=40", col=clust)

#k-means on the W matrix of the NMF (using KL divergence),
clust<-kmeans(nmf.res$W, 20, nstart=20)$cluster
nmi.news(clust)
# 0.9360338
plot(clust, main="k-means on the W matrix of the NMF (using KL divergence) for k=40", col=clust)

#k-means on the W' matrix of factorization W'ΣH' obtained from the NMF, and
clust<-kmeans(nmf.wsh.res$W, 20, nstart=20)$cluster
nmi.news(clust)
# 0.9247936
plot(clust, main="k-means on the W' matrix of factorization W'ΣH' obtained from the NMF for k=40", col=clust)

#k-means on the W' matrix of factorization ΣW'H' obtained from the NMF.
clust<-kmeans(nmf.swh.res$W, 20, nstart=20, iter.max = 15)$cluster
nmi.news(clust)
#0.8291794
plot(clust, main="k-means on the W' matrix of factorization ΣW'H' obtained from the NMF for k=40", col=clust)
## YOUR PART ENDS HERE

