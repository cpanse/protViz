#R

context("centroid") 


# splits by 0 intensities
.groundtruth.peakgroups<- function(x){
	V <- x != 0.0; 
	n <- length(V); 
	pg <- rep(0, n)

	count <- 0; 
	for (i in 1:length(V)){ 
		if (V[i]){
			if(!V[i-1]){count <- count + 1}; 
			pg[i] <- count} 
			}
	pg
}

test_that("compare to ground truth centroids.", {

	p <- .getProfileMS2()
	mZ <- p$mZ
	intensity <- p$intensity

	peakgrps <- split(1:n, .groundtruth.peakgroups(intensity))
	peakgrps <- peakgrps[names(peakgrps) != "0"]

	rv <- lapply(peakgrps, FUN = function(i) {
		intensity.auc <- .trapez(mZ[i], intensity[i])
		mZ.centroid <- weighted.mean(x = mZ[i], w = intensity[i])
		data.frame(mZ = mZ.centroid, intensity = intensity.auc, n=length(i))
	})

	centroid.groundtruth <- do.call("rbind", rv)

	# compare ten most intense peaks
	topN <- 10
	centroid.groundtruth.topN <- centroid.groundtruth[order(centroid.groundtruth$intensity, decreasing=TRUE)[1:topN],]

	rv <- centroid(mZ, intensity, 20)
	rv.topN <- rv[order(rv$intensity, decreasing=TRUE)[1:topN],]


	expect_equal(
	  c(summary(fit <- lm(rv.topN$mZ ~ centroid.groundtruth.topN$mZ))$r.squared,
	  summary(fit <- lm(rv.topN$intensity ~ centroid.groundtruth.topN$intensity))$r.squared),
	  c(1,1), tolerance = 0.001)

})
