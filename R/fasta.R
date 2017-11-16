#R

setMethod('summary', "Fasta",
          function(object, ...) {
            list("filename" = NA,
                 "number_of_IDs" = object$getDesc(),
                 "number_of_REVs" = NA,
                 "number_of_Contaminats" = NA,
                 "number_of_AAs" = NA,
                 "number_of_tryptics_peptides" = object$getNumberOfTrypticPeptides(),
                 "object_size" = object.size(object)
            )
          }          
)

#summary.FASTA <- function(object, revpattern = "^>REV.*", conpattern = "^>.*FGCZCont.*", ...){
  
  
#}

