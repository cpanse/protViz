\name{pgImporter}

\alias{pgImporter}

\title{Progenesis Importer Function}

\description{

    The function reads a CSV file which is exported using ProgenesisLCMS. The object is flexible when it comes to exported fields, ideally, everything is in the CSV. 
    The delimiter is a semicolon (this setting is ProgenesisLCMS version AND system language setting sensitive). All the information is kept and accessible.
    
    There is an option which allows you to switch the separator.

}

\details{


    This importer is dedicated to the commercial software ProgenesisLCMS. Nevertheless it is has a simple structure and can be adapted and used also for other maps.

    todo: Pass an additonal argument to the pgImporter function which 
    allows to define the separator.

}

\usage{

    pgImporter(file, sep=';')

}

\arguments{

  \item{file}{A csv file.}

  \item{sep}{The field separator character. default vaule is set to a semicolon.}

}

\author{

   Christian Panse, Hubert Rehauer, Jonas Grossmann 2012

}
