\name{aa_mw_table}
\alias{aa_mw_table}
\docType{data}
\title{
Table of amino acids molecular weight
}
\description{
This data.frame contains the table of amino acid residues molecular weight (mw).
It  is used to calculate mw of modified residues. The first two columns are
fixed and the rest columns are applied to store modification types. In practice,
this table is copied as template and modified accordingly for each plot.
}
\usage{data("aa_mw_table")}
\format{
    A data frame with 25 observations on the following 7 variables.
    \describe{
        \item{\code{aa}}{a character vector}
        \item{\code{weight}}{a numeric vector}
        \item{\code{aa_varmod}}{a character vector}
        \item{\code{labelmod}}{a logical vector}
        \item{\code{reporterion}}{a logical vector}
        \item{\code{reporterion_group}}{a logical vector}
        \item{\code{labelmod_group}}{a logical vector}
    }
}

\source{
http://www.matrixscience.com/help/aa_help.html
}
\examples{
data(aa_mw_table)
}
\keyword{datasets}
