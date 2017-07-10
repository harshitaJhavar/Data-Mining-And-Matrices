## Data Mining and Matrices
## Summer semester 2017
##
## Assignment 3
## Due 9 July 2017 at 23:59

## Student's information
## Name: Harshita Jhavar
## Matriculation #: 2566267

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

#setwd("/home/harshita/Desktop/DMM/Data Analysis Assignments/Assignment3/assignment_3")

## Preamble
###########

## We use library rworldmap to create pretty maps
## If you don't have it, you have to uncomment the following line
#install.packages("rworldmap")
## If you have issues with rworldmap library, you can omit it (see below)
library(rworldmap)
library(corpcor) #For caculating pseudoinverse
## Same thing for plotrix
#install.packages("plotrix")
library(plotrix)

## We use R pacakge fastICA to provide us a version of FastICA.
## You might have to install it:
#install.packages("fastICA")
library(fastICA)

## utils.R contains some helper functions
## This might require
#install.packages("Matrix")
source("utils.R")

## Task 1
#########

## Implement the two-phase CX algorithm in the following boilerplate
twophase.cx <- function(A, k, fudge=5) {
	# A is the input data
	# k is the final number of columns
	# fudge determines how many columns are sampled in the 
	#       first phase (fudge*k, default: fudge = 5)
	
	## You can use function qr for the second part
	## If qr.out <- qr(M), then qr.out$pivot gives a
	## permutation of the columns of M, in particular
	## qr$pivot[1:k] gives the indices of the first
	## k columns
	
	## YOUR PART STARTS HERE
  SVD<-svd(A)
  V<-SVD$v
  Vt<-t(V)
  p<-array(data = NA, dim = ncol(Vt), dimnames = NULL)
  
  for(i in 1:ncol(Vt)){
    p[i]=(norm(Vt[,i],"2")*norm(Vt[,i],"2"))/k  
  }
  c=k*fudge
  S1<-matrix(data = 0, nrow=ncol(A), ncol=c, dimnames = NULL)
  draw<-array(data = NA, dim = c, dimnames = NULL)
  scale<-array(data = NA, dim = c, dimnames = NULL)
  Vkt<-Vt[1:k,]
  ##Step 1: Random sampling
  for(i in 1:c){
    draw[i]<-sample(1:ncol(Vkt),1,replace=TRUE) 
    S1[draw[i],i]<-1
    scale[i]<-1/sqrt(c*p[draw[i]])  
  }
  D1<-diag(scale)
  R1<-matrix(data = NA, nrow=k, ncol=c, dimnames = NULL)
  R1<-Vkt%*%S1%*%D1
  S2<-matrix(data = 0, nrow=c, ncol=k, dimnames = NULL)
  ##Step 2: QR Factorization
  T<-qr(R1)
  F<-T$pivot
  Fk<-F[1:k]
  draw2<-array(data = NA, dim = k, dimnames = NULL)
  for(i in 1:k) 
  {
    S2[Fk[i],i]=1
  }
  b<-array(data=0,dim=k,dimnames=NULL)
  R2<-A%*%S1%*%S2
  for(i in 1:ncol(A)){
    for(j in 1:k){
      if(R2[,j]==A[,i]){
        b[j]<-i
      }
    }
  }
  X<-pseudoinverse(R2)%*%A
  err<-matrix(data=NA, nrow=nrow(A),ncol=ncol(A),dimnames=NULL)
  err<-norm(A-R2%*%X, 'F')
  return(list(cols=b,X=X,err=err))
	## YOUR PART ENDS HERE
	
	# The function should return a list of three elements:
	# cols : the *indices* of the columns in C (i.e. C <- A[,cols])
	# X    : the matrix X
	# err  : reconstruction error, err <- norm(A - C %*% X, 'F')
}

## This function calls twophase.cx repeatedly and returns the 
## decomposition with least error
cx <- function(A, k, nreps=20, ...) {
	best.res = list(cols=NA, X=NA, err=Inf)
	for (rep in 1:nreps) {
		res <- twophase.cx(A, k, ...)
		if (res$err < best.res$err) {
			best.res <- res
		}
	}
	return(best.res)
}

## Load climate data & coordinates
data <- as.matrix( read.csv("worldclim.csv") )
coord <- read.csv("coordinates.csv")

## Normalize the data to z-scores

## YOUR PART START HERE

data.normal <-scale(data, scale=TRUE, center=TRUE)

## YOUR PART ENDS HERE

## Compute the CX factorization of the *transpose* of the climate
## data for 5 factors
cx.res <- cx(t(data.normal), 5)

## Set up the plotting
xLim <- c(min(coord["lon"]), max(coord["lon"]))
yLim <- c(min(coord["lat"]), max(coord["lat"]))
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)

## Plot the locations of the selected columns
## First column is marked with A, second with B, and so on
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=2, pch=19, cex=1.5)
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=4, pch=65:127, cex=.6)

## Plot the row of X corresponding to the first column
x <- cx.res$X[1,]
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")


## Do similar plots to other rows of X and try to interpret the
## results. Try also CX decompositions with different k and different fudge factors (incl. fudge=1).
## Implement also deterministic.cx that just selects the top-k columns based on their probabilities.

## YOUR PART STARTS HERE
## Plot the row of X corresponding to the second column
x <- cx.res$X[2,]
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")
## Plot the row of X corresponding to the third column
x <- cx.res$X[3,]
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## With rank = 4 and fudge = 8
cx.res <- cx(t(data.normal), 4,8)
xLim <- c(min(coord["lon"]), max(coord["lon"]))
yLim <- c(min(coord["lat"]), max(coord["lat"]))
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=2, pch=19, cex=1.5)
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=4, pch=65:127, cex=.6)
#Plot corresponsing to the first column
x <- cx.res$X[1,]
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## With rank = 5 and fudge = 1
cx.res <- cx(t(data.normal), 5, 1)
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=2, pch=19, cex=1.5)
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=4, pch=65:127, cex=.6)
#Plot corresponsing to the first column
x <- cx.res$X[1,]
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## With rank = 6 and fudge = 6
cx.res <- cx(t(data.normal), 6, 6)
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=2, pch=19, cex=1.5)
points(coord[cx.res$cols,1], coord[cx.res$cols,2], col=4, pch=65:127, cex=.6)
#Plot corresponsing to the first column
x <- cx.res$X[1,]
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## YOUR PART ENDS HERE


## Task 2
#########

## Implement the convex_cone algorithm for NNCX inside this boilerplate
convex.cone <- function(R, k) {
  ## R is the input matrix (assumed to be nonnegative)
  ## k is the rank (i.e. number of columns to choose)
  b<-array(data=0,dim=k,dimnames=NULL)
  c<-matrix(data=NA, nrow=nrow(R),ncol=1,dimnames=NULL) 
  #Finding the column with the highest norm
  norms_columns_R <- list()
  for (i in 1:ncol(R)){
    norms_columns_R[i]<-norm(R[i],'F')}
  #Assigning the column c with the highest norm
  c <- as.matrix(R[,c(which.max(norms_columns_R))])
  #Normalizing c to a unit norm
  c <- c / sqrt(sum(c^2))
  SVD<-svd(R)
  V<-SVD$v
  Vkt <- V[1:k,]
  #Solving for nonnegative x that minimizes ||R-cx(t)||
  print(dim(pseudoinverse(c)))
  print(dim(t(Vkt)))
  x<-pseudoinverse(c)%*%t(Vkt)
  err<-matrix(data=NA, nrow=nrow(R),ncol=ncol(R),dimnames=NULL) 
  err<-norm(t(Vkt)-c%*%x,'F')
  return(list(cols=abs(c),X=x,err=err))
}

#This function calls cx_cone repeatedly and returns the decomposition with least error
cx_cone <- function(A, k, nreps=20, ...){
  R <- A
  best.res = list(cols=NA, x=NA, err=Inf)
  for (rep in 1:nreps) {
    res <- convex.cone(R, k, ...)
    if (res$err < best.res$err) {
      best.res <- res
    }
  }
  return(best.res)
}
## data.normal is the normalized data from the above task
## Shift it to nonnegative values

## YOUR PART STARTS HERE
#Normalizing the world climate data with non-negative values
data.nn <- data.normal
for(j in 1:ncol(data.nn)){
  minimum <- min(data.nn[,c(j)])
  for(i in 1:nrow(data.nn)){
    if(data.nn[i,j]<0){
      data.nn[i, j] <- data.nn[i, j]-minimum 
    }
  }
}  

## YOUR PART ENDS HERE

## Run the convex.cone algorithm with data.nn and different values of k and analyze the
## results based on the reconstruction error and visualizations (use the visualization
## methods from the previous task).

## YOUR PART STARTS HERE
#For k=5
nncx.res <- cx_cone(t(data.nn), k=5)
#Setting the plotting
xLim <- c(min(coord["lon"]), max(coord["lon"]))
yLim <- c(min(coord["lat"]), max(coord["lat"]))
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[nncx.res$cols,1], coord[nncx.res$cols,2], col=2, pch=19, cex=1.5)
points(coord[nncx.res$cols,1], coord[nncx.res$cols,2], col=4, pch=65:127, cex=.6)
#Plot corresponsing to the first column for k=5
x <- nncx.res$X[1,]
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

#For k=4
nncx <- convex.cone(data.nn, k=4)
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[nncx.res$cols,1], coord[nncx.res$cols,2], col=2, pch=19, cex=1.5)
points(coord[nncx.res$cols,1], coord[nncx.res$cols,2], col=4, pch=65:127, cex=.6)
#Plot corresponsing to the first column for k=5
x <- nncx.res$X[1,]
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

#For k=6
nncx <- convex.cone(data.nn, k=6)
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[nncx.res$cols,1], coord[nncx.res$cols,2], col=2, pch=19, cex=1.5)
points(coord[nncx.res$cols,1], coord[nncx.res$cols,2], col=4, pch=65:127, cex=.6)
#Plot corresponsing to the first column for k=5
x <- nncx.res$X[1,]
points(coord[,1], coord[,2], col=color.scale(x, c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(x), 4), round(mean(x), 4), round(max(x), 4)), color.scale(sort(x), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")


## YOUR PART ENDS HERE.

## Task 3
##########

## Load data
house.price <- read.csv("us_housing_prices.csv", row.names=1)

## The places are rows
rownames(house.price)
## And columns are months
colnames(house.price)[1:18]
colnames(house.price)[seq(1, by=12, to=ncol(house.price))]

## Let's plot all the lines
## The plots stay more readable when you plot less than 10 lines at a time
## Notice that this function expects the time series to be the *rows* of the input
plot.time.series(house.price[1:5,])
plot.time.series(house.price[6:10,])
plot.time.series(house.price[11:15,])
plot.time.series(house.price[16:20,])
## You can also try
#plot.time.series(house.price[1:10,], 2)

## Next we need to handle the not-available numbers
## First we'll replace them with 0:
hp <- house.price
hp[is.na(house.price)] <- 0

## Then let's transpose and scale the data
hp <- scale(t(hp))

## Then let's compute the full ICA
## verbose=TRUE will print out the iterations and how well do we converge
## Replacing the "20" with some smaller number k will give k independent components
## by taking only the top-k singular vectors of U as the whitened data.
hp.ica <- fastICA(hp, 20, fun="logcosh", verbose=TRUE)

## YOUR PART: explain what the different matrices in hp.ica mean and show how you can
## reconstruct the (scaled) housing prices in Los Angeles, California.
a<-hp.ica$A
rownames(hp.ica$S) <- colnames(house.price)
s<-hp.ica$S
cala<-s%*%a

## YOUR PART ENDS HERE


## We can plot some of the independent components
## Remember that we have to transpose them for plotting
## We could also change the signs if needed
rownames(hp.ica$S) <- colnames(house.price)
with(hp.ica, plot.time.series(t(S[,1:5])) )

## We can also plot the scatterplot
with(hp.ica, plot(S[,1:2], asp=1) )
## And we can get the labels of the outliers with "identify"
## The command "identify" allows you to click any point in the scatter plot
## and it will annotate it; use this to annotate some outlier points.
## You can end the identification process by pressing any other than the 
## left mouse key. For more information, see
#?identify
with(hp.ica, identify(S[,1:2], labels=rownames(S)) )


## YOUR PART STARTS HERE
##timeseries plots
## CA Los Angeles
with(hp.ica, plot.time.series(t(s[,1:2])) )
s[212:330,2]<-abs(s[212:330,2])
with(hp.ica, plot.time.series(t(s[,1:2])) )
s<-hp.ica$S
## CA San Diego
with(hp.ica, plot.time.series(t(s[,2:3])) )
s[200:330,3]<-abs(s[200:330,3])
with(hp.ica, plot.time.series(t(s[,2:3])) )
s<-hp.ica$S
##Scatter plots
with(hp.ica, plot(s[,1:2], xlab="AZ Phoenix", ylab="CA Los Angeles" ,asp=1) )
with(hp.ica, identify(s[,1:2], labels=rownames(s)) )
#[1] 244 261
with(hp.ica, plot(s[,7:8], xlab="FL Miami", ylab="FL Tampa" ,asp=1) )
with(hp.ica, identify(s[,7:8], labels=rownames(s)) )
#[1] 271 299 329
with(hp.ica, plot(s[,4:5], xlab="CA San Francisco", ylab="CO Denver" ,asp=1) )
with(hp.ica, identify(s[,4:5], labels=rownames(s)) )
#[1]  22 266
with(hp.ica, plot(s[,6:7], xlab="Washington DC", ylab="Miami FL" ,asp=1) )
with(hp.ica, identify(s[,6:7], labels=rownames(s)) )
#[1]  1 271 302
##Depict the economic crisis from 1988-2014
with(hp.ica, plot.time.series(t(s[,2:3])) )
s[212:330,2]<-abs(s[212:330,2])
with(hp.ica, plot.time.series(t(s[,2:3])) )
s<-hp.ica$S
with(hp.ica, plot.time.series(t(s[,2:3])) )
s[200:330,3]<-abs(s[200:330,3])
with(hp.ica, plot.time.series(t(s[,2:3])) )
##AVG_CAL
avgrow<-function(x){
  SUM<-array(0,dim=length(nrow(x)))
  for(i in 1:nrow(x)){
    count=0
    sum=0
    for(j in 1:ncol(x)){
      if(is.na(x[i,j])){
        count=count+1
      }
      else{
        sum=sum+x[i,j]
      }
    }
    SUM[i]<-sum/(ncol(x)-count)
  }
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if(is.na(x[i,j])){
        x[i,j]<-SUM[i]
      }
    }
  }
  x
}
## Observations:
hpt<-avgrow(house.price)
hpts<-scale(t(hpt))
hpt.ica <- fastICA(hpts, 20, fun="logcosh", verbose=TRUE)
rownames(hpt.ica$S) <- colnames(house.price)
st<-hpt.ica$S
with(hpt.ica, plot.time.series(t(st[,1:2])) )
st[250:330,2]<-abs(st[250:330,2])
with(hpt.ica, plot.time.series(t(st[,1:2])) )
st[150:230,2]<-(-1)*(st[150:230,2])
with(hpt.ica, plot.time.series(t(st[,1:2])) )
## Replace 1st
replace1st<-function(x){
  A = list()
  for(i in 1: nrow(x)){
    for(j in 1:ncol(x)){
      if(!is.na(x[i,j])){
        A[i]<-x[i,j]
        break
      }
      
    } 
  }
  
  for(i in 1: nrow(x)){
    for(j in 1:ncol(x)){
      if(is.na(x[i,j])){
        x[i,j]<-A[i]
      }
    } 
  }
  return(x)}
  ## Observations
  hps<-replace1st(house.price)
  hps<-scale(t(hps))
  hps.ica <- fastICA(hps, 20, fun="logcosh", verbose=TRUE)
  rownames(hps.ica$S) <- colnames(house.price)
  sts<-hps.ica$S
  with(hps.ica, plot.time.series(t(sts[,1:2])) )
  ##Whitened data:
  hp <- house.price
  hp[is.na(house.price)] <- 0
  hp <- scale(t(hp))
  SVD<-svd(hp)
  U<-SVD$u
  plot(U[,1:2])
  ## For ICA over U
  rownames(U) <- colnames(house.price)
  hp.w.ica <- fastICA(U, 20, fun="logcosh", verbose=TRUE)
  s.w<-hp.w.ica$S
  with(hp.w.ica, plot.time.series(t(s.w[,1:2])) )
  ##Scree test
  plot(SVD$d)
  #Result comes as 2 from observation
  ##Guttman Kaiser
  GuttmanK<-function(x){
    count<-0
    for(i in 1:length(x)){
      if(x[i]>1){
        count<-count+1
      }
    }
    count
  }
  ##
  GuttmanK(SVD$d)
 # [1] 13
  ##Result taken as 13
  hp.w.ica <- fastICA(U, 13, fun="logcosh", verbose=TRUE)
  rownames(hp.w.ica$S) <- colnames(house.price)
  s.w<-hp.w.ica$S
  with(hp.ica, plot.time.series(t(s.w[,1:2])) )

## YOUR PART ENDS HERE
# Under test.clust.alg(Clustering) value = .9732194
  