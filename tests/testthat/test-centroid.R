#R

context("profile mode to centroid") 


test_that("compare to ground truth (Orbitrap) centroids.", {
	p <- .getProfileMS2()
	mZ <- p$mZ
	intensity <- p$intensity
	n <- length(intensity)

	peakgrps.groundtruth <- split(1:n, .determine.peakgroups.orbitrap(intensity))
	peakgrps.groundtruth <- peakgrps.groundtruth[names(peakgrps.groundtruth) != "0"]

	rv <- lapply(peakgrps.groundtruth, FUN = function(i) {
		intensity.auc <- .trapez(mZ[i], intensity[i])
		mZ.centroid <- weighted.mean(x = mZ[i], w = intensity[i])
		data.frame(mZ = mZ.centroid, intensity = intensity.auc, n=length(i))
	})

	centroid.groundtruth <- do.call("rbind", rv)

	# compare ten most intense peaks
	topN <- 10
	centroid.groundtruth.topN <- centroid.groundtruth[order(centroid.groundtruth$intensity, decreasing=TRUE)[1:topN],]

	rv <- centroid(mZ, intensity, tolppm=20)
	rv.topN <- rv[order(rv$intensity, decreasing=TRUE)[1:topN],]


	expect_equal(
	  c(summary(fit <- lm(rv.topN$mZ ~ centroid.groundtruth.topN$mZ))$r.squared,
	  summary(fit <- lm(rv.topN$intensity ~ centroid.groundtruth.topN$intensity))$r.squared),
	  c(1,1), tolerance = 0.01)

})
