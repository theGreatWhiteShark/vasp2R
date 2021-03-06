% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vasp-utils.R
\name{vasp.diff}
\alias{vasp.diff}
\title{Calculates the difference of two charge densities.}
\usage{
vasp.diff(x, y)
}
\arguments{
\item{x}{First object of class "vasp". It is used as the reference.}

\item{y}{Second object of class "vasp".}
}
\value{
Input x with different values in the 'charge' element.
}
\details{
The charge of the second input is subtracted from the one of the first input.
}
\author{
Philipp Mueller
}
\seealso{
Other vasp: \code{\link{vasp.bonds}},
  \code{\link{vasp.import}}, \code{\link{vasp.reproduce}},
  \code{\link{vasp.rotate.cell}}
}

