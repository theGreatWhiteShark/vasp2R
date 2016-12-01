##' Output of the \code{\link{vasp.import}} function applied on the CHGCAR file of a single silicium atom in a fcc cell.
##'
##' @format Data frame containing four columns
##' \itemize{
##'   \item{ charge: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to a grid point and "charge" containing the charge density on this specific grid point.}
##'   \item{ atoms: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to the atoms within the cell and "type" providing a factor to distinguish the individual types of atoms.}
##'   \item{ lattice:  Consists of three columns: "x", "y", "z" containing the components of the three lattice vectors of the unit cell.}
##' }
##'
##' @name SiChgcar
NULL

##' Output of the \code{\link{vasp.import}} function applied on the CHGCAR file of a carbon monooxide molecule attached to a Ni surface. Example from \url{http://vasp.at/vasp-workshop/slides/handsonIII.pdf}
##'
##' @format Data frame containing four columns
##' \itemize{
##'   \item{ charge: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to a grid point and "charge" containing the charge density on this specific grid point.}
##'   \item{ atoms: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to the atoms within the cell and "type" providing a factor to distinguish the individual types of atoms.}
##'   \item{ lattice:  Consists of three columns: "x", "y", "z" containing the components of the three lattice vectors of the unit cell.}
##' }
##'
##' @name NiChgcar
NULL
