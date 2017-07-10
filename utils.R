require(Matrix)

## plot rows of a matrix as separate series
plot.time.series <- function(A, cols=1) {
	k <- nrow(A)
	x <- 1:ncol(A)
	mfrow.save <- par("mfrow")
	mar.save <- par("mar")
	par(mfrow=c(k,cols))
	labseq <- seq(1, by=12, ncol(A))
	# Reduce margin from top
	par(mar=c(5,4,1,2) + 0.1)
	for (i in 1:(k-cols)) {
		plot(x, t(A[i,]), type="l", xaxt="n", xlab="", ylab=rownames(A)[i])
		axis(1, at=labseq, labels=NA)
	}
	for (i in (k - cols + 1):k) {
		plot(x, t(A[i,]), type="l", xaxt="n", xlab="", ylab=rownames(A)[i])
		axis(1, at=labseq, labels=colnames(A)[labseq], las=2)
	}
	# Set parameters back
	par(mfrow=mfrow.save)
	par(mar=mar.save)
}

