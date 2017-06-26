# protViz - Visualizing and Analyzing Mass Spectrometry Related Data in Proteomics

## Documentation

The package ships with two pdf vignettes.

```
vignette('protViz')
vignette('PTM_MarkerFinder')
```

## Installation

### CRAN

```
install.packages('protViz')
```


### Docker

```
docker pull cpanse/protviz \
&& docker run -d -p 8792:8787 cpanse/protviz     
```

connect to http://yourdockerhost:8791  using a web browser

* username: rstudio
* password: rstudio


## Documentation

The package ships with a package vignette (browseVignettes('protViz') and a reference manual (just type ?protViz on the R shell).

Both documents are also available on the [package's CRAN](https://CRAN.R-project.org/package=protViz) page.


## Related approaches

* [RforProteomics](http://bioconductor.org/packages/RforProteomics/)

[![Research software impact](http://depsy.org/api/package/cran/protViz/badge.svg)](http://depsy.org/package/r/protViz)
