% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vasp-import.R
\name{vasp.import}
\alias{vasp.import}
\title{Import of an CHGCAR or CHG file.}
\usage{
vasp.import(path)
}
\arguments{
\item{path}{Character vector containing the path to the CHGCAR file. (As for now you MUST to include in the name of the file at the end of 'path' and it MUST contain the letters 'CHG' if the format is according to the CHG file produced by VASP or 'CHGCAR' if according to the CHGCAR file)}
}
\value{
List of three named data.frames:
\itemize{
  \item{ charge: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to a grid point and "charge" containing the charge density on this specific grid point.}
  \item{ atoms: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to the atoms within the cell and "type" providing a factor to distinguish the individual types of atoms.}
  \item{ lattice:  Consists of three columns: "x", "y", "z" containing the components of the three lattice vectors of the unit cell.}
 }
}
\details{
The CHGCAR file got a very odd format which has to be reformatted to ensure the its compatibility with R. In it the number of grid points of the axes are given in the 12th row of the file. In the actual charge density following from the 13th row onward the x coordinate is the fastest and the z the slowest varying index. If you want to assign a scale to your results, keep in mind that the charge density is divided by the volume of the unit cell. The individual lattice vectors in the CHGCAR file are given in different rows and their x, y and z are provided in the different columns.
}
\author{
Philipp Mueller
}
\seealso{
Other vasp: \code{\link{vasp.bonds}},
  \code{\link{vasp.diff}}, \code{\link{vasp.reproduce}},
  \code{\link{vasp.rotate.cell}}
}

