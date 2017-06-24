## Data Mining and Matrices
## Summer semester 2017
##
## Assignment 1
## Due 28 May 2017 at 23:59

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

## Preamble
###########

## We use library rworldmap to create pretty maps
## If you don't have it, you have to uncomment the following line
#install.packages("rworldmap")
## If you have issues with rworldmap library, you can omit it (see below)
library(rworldmap)

## Same thing for plotrix
#install.packages("plotrix")
library(plotrix)

## utils.R contains some helper functions
## This might require
#install.packages("clue")
source("utils.R")

## Task 1
##########

## Load the data matrix and coordinates
data <- as.matrix( read.csv("worldclim.csv") )
coord <- read.csv("coordinates.csv")

## Familiarize yourself with the data here e.g.
#summary(data)
#colnames(data)
## You should not include this part to your report. 

## Compute the SVD and extract the matrices
SVD = svd(data)
U = SVD$u
S = SVD$d
V = SVD$v

## Prepare for plotting: the extremes of the map
xLim <- c(min(coord["lon"]), max(coord["lon"]))
yLim <- c(min(coord["lat"]), max(coord["lat"]))

## Get and plot a map
## If you couldn't get rworldmap working, ignore these lines
map <- getMap(resolution="low")
plot(map, xlim=xLim, ylim=yLim, asp=1)

## Plot the first column of U. The color indicates the value, with
## red being low, green being middle, and blue being high
## Try
#?color.scale
## for more information about the color scale.
## If you don't have working rworldmap, replace 'points' with 'plot'
points(coord[,1], coord[,2], col=color.scale(U[,1], c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), color.spec="rgb"))

## Alternative plot with different color scheme and filled circles
points(coord[,1], coord[,2], col=color.scale(U[,1], c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)

## A color legend to explain the colors
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(U[,1]), 4), round(mean(U[,1]), 4), round(max(U[,1]), 4)), color.scale(sort(U[,1]), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## Plot the second column
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=color.scale(U[,2], c(0,1), 0.8, 1, color.spec="hsv"))
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(U[,2]), 4), round(mean(U[,2]), 4), round(max(U[,2]), 4)), color.scale(sort(U[,2]), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## YOUR PART STARTS HERE
#Normalizing the data to Z-Scores
#Z-Score is (x-mean(x))/sd(x)) where x is any value in a particular column of the matrix and mean and sd are defined with respect to every column in a matrix
#scale command directly computes the Z-Score center: when center = TRUE, the objects’ column means are 
#subtracted from the values in those columns (ignoring NAs)
#if scale = TRUE, the centered column values are divided by the column’s standard deviation (when center is also TRUE)
normalized_data <- scale(data, center = TRUE, scale = TRUE)
#Computing the SVD of the normalized data and extracting the matrices
SVD_of_normalized_data = svd(normalized_data)
U_of_SVD_from_normalized_data = SVD_of_normalized_data$u
S_of_SVD_from_normalized_data = SVD_of_normalized_data$d
V_of_SVD_from_normalized_data = SVD_of_normalized_data$v

#Plotting the first column of U for the normalized data
map_normalized <- getMap(resolution="low")
plot(map_normalized, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=color.scale(U_of_SVD_from_normalized_data[,1], c(0,1), 0.8, 1, color.spec="hsv"), cex=.6, pch=19)
## A color legend to explain the colors
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(U_of_SVD_from_normalized_data[,1]), 4), round(mean(U_of_SVD_from_normalized_data[,1]), 4), round(max(U_of_SVD_from_normalized_data[,1]), 4)), color.scale(sort(U_of_SVD_from_normalized_data[,1]), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## Plot the second column
plot(map_normalized, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=color.scale(U_of_SVD_from_normalized_data[,2], c(0,1), 0.8, 1, color.spec="hsv"))
color.legend(xLim[1]+1, yLim[1]-5, xLim[2]-1, yLim[1]-3, c(round(min(U_of_SVD_from_normalized_data[,2]), 4), round(mean(U_of_SVD_from_normalized_data[,2]), 4), round(max(U_of_SVD_from_normalized_data[,2]), 4)), color.scale(sort(U_of_SVD_from_normalized_data[,2]), c(0,1), 0.8, 1, color.spec="hsv"), gradient="x")

## YOUR PART ENDS HERE

#############
## Task 2
#############

## You have to implement the rank selection techniques yourself. 
## Remember to use the normalized data.

## If S contains the singular values, you can plot them by
plot(S, type="l")

## YOUR PART STARTS HERE
#Computing the SVD of the normalized data
SVD_of_normalized_data = svd(normalized_data)
U_of_SVD_from_normalized_data = SVD_of_normalized_data$u
S_of_SVD_from_normalized_data = SVD_of_normalized_data$d
V_of_SVD_from_normalized_data = SVD_of_normalized_data$v

#Selecting k with Guttman-Kaiser Criterion
k_GKC = 0 #Rank variable
for(i in 1:length(S_of_SVD_from_normalized_data)){
  if(S_of_SVD_from_normalized_data[i] < 1){
    break #Followed slide11 from Chapter 1
  }
  else{
    k_GKC = k_GKC + 1
  }
}
print(paste0(k_GKC, " is rank with Guttman-Kaiser Criterion."))
#90% of explained Frobenius norm
k_90_Frob = 0 #Rank variable
Square_of_S_of_SVD_from_normalized_data = S_of_SVD_from_normalized_data ^2 #Storing the square of all singular values
threshold = 0.9 * sum(Square_of_S_of_SVD_from_normalized_data) #90% of the sum of the squares of all the singular values
#Defining a function to find sum of first n elements in a vector
sumfun<-function(x,start,end){
  return(sum(x[start:end]))
}
temp=0 #temporary variable to hold sum
for(i in 1:length(Square_of_S_of_SVD_from_normalized_data)){
  temp = sumfun(Square_of_S_of_SVD_from_normalized_data,1,i) #Sum of square of singular values from 1 to i
  if(temp >= threshold){
    k_90_Frob = i
    break
  }
}
print(paste0(k_90_Frob, " is rank with 90% of explained Frobenius Norm"))
#Scree test
#Plotting the singular values in decreasing order
plot(sort(S_of_SVD_from_normalized_data, decreasing = TRUE), type="l", xlab = "Index of SVD Values", ylab = "SVD values of normalized dataset in decreasing order",main = "Scree plot")
#This scree plot shows that 6 of those factors explain most of the variability because the line starts to straighten after factor 6.
#The remaining factors explain a very small proportion of the variability and are likely unimportant.
print("6 is rank with Scree Test")
#Entropy Based Method
relative_contribution = NULL
sum_of_product_relative_contribution = NULL
for(i in 1:length(Square_of_S_of_SVD_from_normalized_data)){
  relative_contribution[i] = Square_of_S_of_SVD_from_normalized_data[i] / sum(Square_of_S_of_SVD_from_normalized_data)
  #Calculating sum of ri*log (ri)
  if(relative_contribution[i] == 0){
    sum_of_product_relative_contribution[i] = 0
  }
  else{
    sum_of_product_relative_contribution[i] = (relative_contribution[i] * log(relative_contribution[i]))
  }
}
entropy = -1 * (1/log(48)) * sum(sum_of_product_relative_contribution)
k_entropy = 0 #Rank Variable
for(i in 1:length(Square_of_S_of_SVD_from_normalized_data)){
  if(sumfun(relative_contribution,1,i) < entropy)
    break
  else{
    k_entropy = i
  }
}
print(paste0(k_entropy, " is rank with entropy based method"))
#Random flipping of signs
row<-nrow(normalized_data)
col<-ncol(normalized_data)
SVD<-svd(normalized_data)
U<-SVD$u
D<-diag(SVD$d)
V<-SVD$v
average_over_different_norm_values<-array(data=0,dim=col-1,dimnames=NULL)
#Flip over 10 random matrices
for(i in 1:10){
  rnd<-sample(c(1,-1),row*col,TRUE,NULL)
  rnd_matrix<-matrix(rnd,nrow=row,ncol=col)
  rnd_matrix<-rnd_matrix*normalized_data
  svd_rnd<-svd(rnd_matrix)
  rnd_U<-svd_rnd$u
  rnd_D<-diag(svd_rnd$d)
  rnd_V<-svd_rnd$v
  different_norm_values<-array(data=NA,dim=col-1,dimnames=NULL)
  for(k in 1:(col-1)){
    residual_normalized_data<-U[,k:col]%*%(D[k:col,k:col])%*%V[k:col,]
    residual_rnd<-rnd_U[,k:col]%*%rnd_D[k:col,k:col]%*%rnd_V[k:col,]
    different_norm_values[k]<-(norm(residual_normalized_data,"2")-norm(residual_rnd,"2"))/norm(residual_normalized_data,"F")
    }
  average_over_different_norm_values<-(average_over_different_norm_values + different_norm_values)
}
average_different_norm_values<-(average_over_different_norm_values)/10
k_randomflip = which.min(abs(average_different_norm_values)) #Taking the minimum value for k
print(paste0(k_randomflip," is the rank by random flips"))

## YOUR PART ENDS HERE

#############
## Task 3
#############

## We again use the normalized data. If that is contained in matrix D
## we can compute the k-means to 5 clusters with 10 re-starts as
climate.clustering <- kmeans(normalized_data, 5, iter.max=100, nstart=10)$cluster

## The plotting is as in Task 1 (we assume map and xLim and yLim are 
## ready)
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=climate.clustering)

## If we want to compare two clusterings, we should match their labels.
## The utils.R file contains match() function to align the labels of
## two clusterings. Here's how to use it

## Open a new device window, e.g. (in Linux)
#X11()

## Compute a new clustering with just one repetition
#other.clustering <- kmeans(D, 5, iter.max=100, nstart=1)$cluster

## Plot the clustering
#plot(map, xlim=xLim, ylim=yLim, asp=1)
#points(coord[,1], coord[,2], col=other.clustering)

## Match the labels
#oc.matched = match(other.clustering, to=climate.clustering)

## Plot the matched labels
#plot(map, xlim=xLim, ylim=yLim, asp=1)
#points(coord[,1], coord[,2], col=oc.matched)


## YOUR PART STARTS HERE
#A new plot in which the x-axis position is the first left singular vector, 
#the y-axis position is the second left singular vector and the color of marker is defined by clustering.
plot(U_of_SVD_from_normalized_data[,1],U_of_SVD_from_normalized_data[,2],type = "p", col = climate.clustering, xlab = "First Left Singular Vector", ylab = "Second left Singular Vector")

#Karhumen Loeve transformation or PCA transformation
#k=2 as we have to project the data into two dimensional subspace
truncated_V_of_SVD_from_normalized_data = matrix(c(V_of_SVD_from_normalized_data[,1],V_of_SVD_from_normalized_data[,2]),nrow=48,ncol=2)
normalized_data_PCA_transformed = normalized_data %*% truncated_V_of_SVD_from_normalized_data
#Applying kMeans on the PCA transformed data
climate_clustering_PCA <- kmeans(normalized_data_PCA_transformed, 5, iter.max=100, nstart=10)$cluster
#Plotting this clustering on world map
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=climate_clustering_PCA)
#Plotting of its first two left singular vectors on x and y axis respectively
SVD_of_normalized_data_PCA_transformed = svd(normalized_data_PCA_transformed)
U_of_SVD_from_normalized_data_PCA_transformed = SVD_of_normalized_data_PCA_transformed$u
S_of_SVD_from_normalized_data_PCA_transformed = SVD_of_normalized_data_PCA_transformed$d
V_of_SVD_from_normalized_data_PCA_transformed = SVD_of_normalized_data_PCA_transformed$v
plot(U_of_SVD_from_normalized_data_PCA_transformed[,1],U_of_SVD_from_normalized_data_PCA_transformed[,2],type = "p", col = climate_clustering_PCA, xlab = "First Left Singular Vector after PCA transformation", ylab = "Second left Singular Vector after PCA transformation")
#Comparing the two clusterings
X11()
## Compute a new clustering with just one repetition
climate_clustering_PCA_1_repetition <- kmeans(normalized_data_PCA_transformed, 5, iter.max=100, nstart=1)$cluster

## Plot the clustering
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=climate_clustering_PCA_1_repetition)

## Match the labels
oc.matched = match(climate_clustering_PCA_1_repetition, to=climate.clustering)

## Plot the matched labels
plot(map, xlim=xLim, ylim=yLim, asp=1)
points(coord[,1], coord[,2], col=oc.matched)

## YOUR PART ENDS HERE
