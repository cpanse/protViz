#R

# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/R/peakplot.R $
# $Id: peakplot.R 6683 2014-09-18 06:52:41Z cpanse $
# $Date: 2014-09-18 08:52:41 +0200 (Thu, 18 Sep 2014) $

# TODO: make a class peakplot
# peakplot.print
# peakplot.plot

.peakplot.putlabel <- function(MASS, INTENSITY, LABEL, 
    l.col="green", 
    delta=0, 
    yMin, 
    maxIntensity=max(INTENSITY)) {

	noise<-seq(0, 2 * delta, length=length(MASS))

	segments(MASS, 
		INTENSITY+0.03*maxIntensity,
		MASS, 
		1.1*maxIntensity,
        lty=2,
		#yMin+noise,
		col=l.col, 
		pch=5, 
		cex=0.5,
		lwd=0.5)

	segments(MASS, INTENSITY, MASS, rep(0, length(MASS)), lwd=1.5, col=l.col)


	text(MASS, yMin+noise, round(MASS,2), 
		cex=0.50,
		pos=1,
		#offset=0.0,
		col="red",
		srt=90)
}

.peakplot.label <- function(spec, match, itol=0.6, ...){
    # filtering the ions
    # TODO(cp): 1. assign highest peak within itol range; 2. col/pch setting

    LABEL.abc<-(abs(match$mZ.Da.error) < itol) & (regexpr("[abc].*", match$label) > 0)
    LABEL.xyz<-(abs(match$mZ.Da.error) < itol) & (regexpr("[xyz].*", match$label) > 0)

    points(spec$mZ[match$idx[LABEL.abc]], spec$intensity[match$idx[LABEL.abc]], col="black", ...)
    points(spec$mZ[match$idx[LABEL.abc]], spec$intensity[match$idx[LABEL.abc]], col="black", type='h')
    points(spec$mZ[match$idx[LABEL.xyz]], spec$intensity[match$idx[LABEL.xyz]], col="blue", ...)
    points(spec$mZ[match$idx[LABEL.xyz]], spec$intensity[match$idx[LABEL.xyz]], col="blue", type = 'h' )
}

.peakplot.pie <- function(spec, match){ 

    LABEL.abc<-abs(match$mZ.Da.error < 0.6) & (regexpr("[abc].*", match$label) > 0)
    LABEL.xyz<-abs(match$mZ.Da.error < 0.6) & (regexpr("[xyz].*", match$label) > 0)

    i.abc<-spec$intensity[match$idx[LABEL.abc]]
    i.xyz<-spec$intensity[match$idx[LABEL.xyz]]

    l.abc<-match$label[LABEL.abc]
    l.xyz<-match$label[LABEL.xyz]

    i.rest<-sum(spec$intensity)-sum(i.abc)-sum(i.xyz)

    pie(c(i.abc,i.xyz,i.rest), c(l.abc, l.xyz, "rest"), col=c(rep("blue",length(i.abc)), rep("grey",length(i.abc)), "white"))
}                      

peakplot <- function(peptideSequence,
    spec, 
    FUN = defaultIon, 
    fi = fragmentIon(peptideSequence, FUN=FUN)[[1]],
    sub = paste(peptideSequence, spec$title, sep=" / "),
    itol = 0.6,
    pattern.abc = "[abc].*",
    pattern.xyz = "[xyz].*",
    ion.axes = TRUE, ...){ 

    # n <- nchar(peptideSequence)
    m <- psm(peptideSequence, spec, FUN, fi = fi, plot = FALSE)
    max.intensity <- max(spec$intensity, na.rm=TRUE)

    plot(spec$mZ,
    	spec$intensity,
        xlab = 'm/z',
        ylab = 'intensity',
        type = 'h',
        sub = sub,
        axes='F', ...
    ) 

    LABEL.abc <- (abs(m$mZ.Da.error) < itol) & (regexpr(pattern.abc, m$label) > 0)
    LABEL.xyz <- (abs(m$mZ.Da.error) < itol) & (regexpr(pattern.xyz, m$label) > 0)

    if (ion.axes){
        if (length(m$idx[LABEL.abc]) > 0){
            axis(1, spec$mZ[m$idx[LABEL.abc]], m$label[LABEL.abc],las=2)
        }
        axis(2)
        if (length(m$idx[LABEL.xyz]) > 0){
            axis(3, spec$mZ[m$idx[LABEL.xyz]], m$label[LABEL.xyz], col.axis='blue', las=2)
        }
    }else{
        axis(1)
        axis(2)
        a.at <- spec$mZ[m$idx[LABEL.abc | LABEL.xyz]]
        a.label <- m$label[LABEL.abc | LABEL.xyz]

        if (length(a.at) > 0) {
            axis(3,a.at, a.label, col.axis='black', las=2)
        } else {
            print ("WARNING")
            print (a.at)
            print (a.label)
        }
    }
    box()

    axis(4, seq(0, max.intensity, length = 6), seq(0, 100, length = 6))

    .peakplot.label(spec = spec, match = m, itol = itol, pch = 22)

    sortedFragmentIonsTable <- data.frame(
      label = c(m$label[LABEL.abc], m$label[LABEL.xyz]),
      mass = c(spec$mZ[m$idx[LABEL.abc]], spec$mZ[m$idx[LABEL.xyz]]))

    if (nrow(sortedFragmentIonsTable) > 0){
	    sortedFragmentIonsTable <- sortedFragmentIonsTable[order(sortedFragmentIonsTable$mass), ]
	    legend("right", sprintf("% 10.3f   %s", sortedFragmentIonsTable$mass, sortedFragmentIonsTable$label),
	      title = "Fragment Ions",
	      bty = 'n',
	      cex = 0.65)
    }

    return(m)
}    


.featureDensityPlot <- function(data, n=ncol(data), nbins = 30){
  my.col<-rainbow(n);
  mids<-numeric()
  density<-numeric()
  for (i in 1:n) { 
    h<-hist(data[,i],nbins, plot=F)
    mids<-c(mids, h$mids)
    density<-c(density, h$density)
  }
  plot(mids,density, type='n')
  for (i in 1:n) { 
    h<-hist(data[,i],nbins, plot=F)
    lines(h$mids,h$density, col=my.col[i])
  }
  legend("topleft", names(data), cex=0.5,
         text.col=my.col
  )
}


