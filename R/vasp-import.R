##' @title Import of an CHGCAR or CHG file.
##'
##' @details The CHGCAR file got a very odd format which has to be reformatted to ensure the its compatibility with R. In it the number of grid points of the axes are given in the 12th row of the file. In the actual charge density following from the 13th row onward the x coordinate is the fastest and the z the slowest varying index. If you want to assign a scale to your results, keep in mind that the charge density is divided by the volume of the unit cell. The individual lattice vectors in the CHGCAR file are given in different rows and their x, y and z are provided in the different columns.
##'
##' @param path Character vector containing the path to the CHGCAR file. (As for now you MUST to include in the name of the file at the end of 'path' and it MUST contain the letters 'CHG' if the format is according to the CHG file produced by VASP or 'CHGCAR' if according to the CHGCAR file)
##'
##' @family vasp
##' @export
##' @useDynLib vasp2R
##' 
##' @return List of three named data.frames:
##' \itemize{
##'   \item{ charge: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to a grid point and "charge" containing the charge density on this specific grid point.}
##'   \item{ atoms: Consists of four columns: "x", "y", "z" containing the Cartesian coordinates to the atoms within the cell and "type" providing a factor to distinguish the individual types of atoms.}
##'   \item{ lattice:  Consists of three columns: "x", "y", "z" containing the components of the three lattice vectors of the unit cell.}
##'  }
##' @author Philipp Mueller
vasp.import <- function( path ){   
    if ( !grep( "CHG", path ) )
        stop( "The CHGCAR or CHG file has to be provided for!" )    
    lattice.constant <- read.table( path, skip = 1, nrows = 1 )
    lattice.vectors <- read.table( path, skip = 2, nrows = 3 )* lattice.constant[[ 1 ]]
    names( lattice.vectors ) <- c( "x", "y", "z" )
    lattice.types <- read.table( path, skip = 5, nrows = 1, as.is = TRUE )
    lattice.numbers <- read.table( path, skip = 6, nrows = 1 )
    lattice.coordinatesystem <- read.table( path, skip = 7, nrows = 1, as.is = TRUE )
    if ( lattice.coordinatesystem[[ 1 ]] != "Direct" && lattice.coordinatesystem[[ 1 ]] != "D" )
        stop( "vasp.import: Wait. The charge density and the coordinates in the CHGCAR are not in direct form? This wasn't accounted in the import function. Please check!" )
    lattice.positions <- read.table( path, skip = 8, nrows = sum( lattice.numbers ) )
    lattice.grid <- read.table( path, skip = ( 9 + sum( lattice.numbers ) ), nrows = 1 )
    lattice.volume <- ( lattice.vectors[ 1, 1 ]^2 + lattice.vectors[ 1, 2 ]^2 + lattice.vectors[ 1, 3 ]^2 )*
        ( lattice.vectors[ 2, 1 ]^2 + lattice.vectors[ 2, 2 ]^2 + lattice.vectors[ 2, 3 ]^2 )*
        ( lattice.vectors[ 3, 1 ]^2 + lattice.vectors[ 3, 2 ]^2 + lattice.vectors[ 3, 3 ]^2 )

    ## Determining the number of lines which span the entries of the charge density
    charge.entries <- lattice.grid[[ 1 ]]* lattice.grid[[ 2 ]]* lattice.grid[[ 3 ]]
    if ( grep( "CHG", path ) ){
        ## CHG was supplied
        charge.lines <- charge.entries/ 10
    } else {
        ## CHGCAR was supplied
        charge.lines <- charge.entries/ 5
    }

    ## Import the raw format of the charge grid.
    charge.raw <- readr::read_table( path, skip = ( 10 + sum( lattice.numbers ) ),
                                    n_max= charge.lines, col_names = FALSE,
                                    progress = FALSE )/ lattice.volume
    if ( grep( "CHG", path ) ){
        ## This ensures the functionality of the function even when the CHG file is provided.
        names( charge.raw ) <- c( "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10" )
    } else {
        names( charge.raw ) <- c( "V1", "V2", "V3", "V4", "V5" )
    }
    ## The main difference between the CHGCAR and the CHG file is the format. The charge density is stored with higher precision in the CHGCAR and with 5 columns per row. The CHG contains 10 columns per row.
    ## After this step there is just one format regardless of the input file.
    ## The rows of the charge density will be concatenated one by one. So A[a,b;c,d] becomes V[a,b,c,d]. 
    if ( grep( "CHG", path ) ){
        charge.vector <- .Call( "importCHG", PACKAGE = "vasp2R",
                               charge.raw )
    } else
        charge.vector <- .Call( "importCHGCAR", PACKAGE = "vasp2R",
                               charge.raw )

    ## Reformatting the charge density
    ## This assumes that the FFT grid is spread homogeneously on the axes 
    x.number <- lattice.grid[[ 1 ]]
    x.range <- seq( 0, 1, , x.number )
    y.number <- lattice.grid[[ 2 ]]
    y.range <- seq( 0, 1, , y.number )
    z.number <- lattice.grid[[ 3 ]]
    z.range <- seq( 0, 1, , z.number )
    ## The elements in the *.range vectors are multiplied with the lattice vectors to determine the actual position of the elements.

    ## Converts the charge entries of charge.vector into a data frame with the Cartesian coordinates in the first three columns and the corresponding charge density of the grid point in the fourth column.
    charge.df <- .Call( "formatCharge", PACKAGE = "vasp2R",
                       charge.vector, x.range, y.range, z.range,
                       lattice.vectors )
    
    ## So why not extracting the positions of the individual atoms from the CHGCAR too?
    atoms.positions <- data.frame( x = rep( 0, nrow( lattice.positions ) ),
                                     y = rep( 0, nrow( lattice.positions ) ),
                                     z = rep( 0, nrow( lattice.positions ) ) )
    ## For dummies (since the day just began)
    lattX <- lattice.vectors$x
    lattY <- lattice.vectors$y
    lattZ <- lattice.vectors$z
    xGrid <- lattice.positions[[ 1 ]]
    yGrid <- lattice.positions[[ 2 ]]
    zGrid <- lattice.positions[[ 3 ]]
    for ( ii in 1 : nrow( lattice.positions ) ){
        atoms.positions$x[[ ii ]] <- lattX[[ 1 ]]* xGrid[ ii ] + 
            lattX[[ 2 ]]* yGrid[ ii ] + lattX[[ 3 ]]* zGrid[ ii ]
        atoms.positions$y[[ ii ]] <- lattY[[ 1 ]]* xGrid[ ii ] + 
            lattY[[ 2 ]]* yGrid[ ii ] + lattY[[ 3 ]]* zGrid[ ii ]
        atoms.positions$z[[ ii ]] <- lattZ[[ 1 ]]* xGrid[ ii ] + 
            lattZ[[ 2 ]]* yGrid[ ii ] + lattZ[[ 3 ]]* zGrid[ ii ]
        
        ## To assure that the atom lies within the supercell
        while( atoms.positions$y[[ ii ]] < 0 ){
            atoms.positions$y[[ ii ]] <- atoms.positions$y[[ ii ]] + ( lattY[[ 1 ]] + lattY[[ 2 ]] )
            atoms.positions$x[[ ii ]] <- atoms.positions$x[[ ii ]] + lattX[[ 2 ]]
        }
        while ( atoms.positions$y[[ ii ]]* 1.1 > ( lattY[[ 1 ]] + lattY[[ 2 ]] ) ) {
            atoms.positions$y[[ ii ]] <- atoms.positions$y[[ ii ]] - ( lattY[[ 1 ]] + lattY[[ 2 ]] )
            atoms.positions$x[[ ii ]] <- atoms.positions$x[[ ii ]] - lattX[[ 2 ]]
        }
        while( atoms.positions$x[[ ii ]] < 0 )
            atoms.positions$x[[ ii ]] <- atoms.positions$x[[ ii ]] + lattX[[ 1 ]]
        while ( atoms.positions$x[[ ii ]]* 1.1 > ( lattX[[ 1 ]] + abs( lattX[[ 2 ]]* atoms.positions$y[[ ii ]]/ lattY[[ 2 ]] ) ) )
            atoms.positions$x[[ ii ]] <- atoms.positions$x[[ ii ]] - lattX[[ 1 ]]
    }
    
    ## The position of the kinds of atoms will be stored in one single data frame with a column named "type" indicating the kind of the atom
    atoms.positions$type <- "place holder"
    lattice.numbers.count <- 0
    for ( tt in 1 : length( lattice.numbers ) ){
        if ( tt == 1 ){
            atoms.positions$type[ 1 : lattice.numbers[[ tt ]] ] <- as.character( lattice.types[[ tt ]] )
            lattice.numbers.count <- lattice.numbers[[ tt ]]
        } else {
            atoms.positions$type[ ( lattice.numbers.count + 1 ) :
                                     ( lattice.numbers.count + lattice.numbers[[ tt ]] ) ] <-  as.character( lattice.types[[ tt ]] )
            lattice.numbers.count <- lattice.numbers.count + lattice.numbers[[ tt ]]
        }
    }
    output <- list( charge = charge.df, atoms = atoms.positions, lattice = lattice.vectors )
    class( output ) <- c( "vasp", "list" )
    
    return( output )
}
