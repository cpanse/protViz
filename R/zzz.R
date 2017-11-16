#R

.onAttach <- function(lib, pkg){
  if(interactive()){
    version <- packageVersion('protViz')
    packageStartupMessage("Package 'protViz' version ", version)
    
    invisible()
  }
}

loadModule("FastaMod", TRUE)
