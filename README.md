[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/protViz)](https://cran.r-project.org/package=protViz)
[![](https://cranlogs.r-pkg.org/badges/grand-total/protViz)](https://cran.r-project.org/package=protViz)
[![](https://cranlogs.r-pkg.org/badges/protViz)](https://cran.r-project.org/package=protViz)
[![codecov](https://codecov.io/github/cpanse/protViz/branch/master/graph/badge.svg?token=lYkh9C6v5O)](https://app.codecov.io/github/cpanse/protViz)
# protViz - Visualizing and Analyzing Mass Spectrometry Related Data in Proteomics

## Documentation

The package ships with two pdf vignettes.

```
browseVignettes('protViz')

vignette('protViz')
vignette('PTM_MarkerFinder')
```

## Installation

### CRAN

```
install.packages('protViz')
```


### from [github](https://github.com/cpanse/protViz)

install the latest development version

```{r}
install.packages('devtools')
library(devtools)
install_git('https://github.com/cpanse/protViz', build_vignettes = FALSE, quiet = FALSE)
library(protViz)
```

or

```{r}
BiocManager::install('cpanse/protViz')
```

### R CMD build hints


```{r}
Rcpp::compileAttributes()

tools::package_native_routine_registration_skeleton(".", character_only = FALSE)


RcppExport SEXP _rcpp_module_boot_MyModule();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_MyModule", (DL_FUNC) &_rcpp_module_boot_MyModule, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_minModuleEx(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

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

* [RforProteomics](https://bioconductor.org/packages/RforProteomics/)
* [Spectra](https://github.com/rformassspectrometry/Spectra)
