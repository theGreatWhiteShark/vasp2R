##' @title vasp.rotate.cell
##' @description Rotates the content of a cell with a specific angle.
##'
##' @details In this version of the function only a rotation in the x-y plane is possible. The supplied content of the cell has to consist of named columns "x" and "y" for the Cartesian coordinates. If no angle is provided the angle between the two lattice vectors is determined and the content will be shifted by half of this angle in the positive direction. A rotation matrix is applied to all the coordinates in the cell using a cpp loop. Generic function working with both the specific content of the cell, like the atom positions, and an object of class "vasp".
##'
##' @param cell Data frame containing at least two named columns "x" and "y" containing Cartesian coordinates or an object of class "vasp".
##' @param lattice Data frame of dimension 3x3 containing the three lattice vectors in separate rows and their coordinates in different columns. See \code{\link{vasp.import}}. If cell is an object of class "vasp" this argument can is omitted.
##' @param angle Specifies the angle with which the cell is going to be rotated. If no value is supplied the angle between the two lattice vectors is determined and the content will be shifted by half of this angle in the positive direction.
##' 
##' @family vasp
##' @seealso \code{\link{vasp.import}}
##' @return Same object like the input "cell" but with modified coordinates.
##' @author Philipp Mueller
vasp.rotate.cell <- function( cell, lattice = NULL, angle = NULL ){
    UseMethod( "vasp.rotate.cell" )
}
vasp.rotate.cell.vasp <- function( cell, lattice = NULL, angle = NULL ){
    atoms.new <- vasp.rotate.cell.default( cell$atoms, cell$lattice, angle )
    charge.new <- vasp.rotate.cell.default( cell$charge, cell$lattice, angle )

    cell.new <- list( atoms = atoms.new, charge = charge.new, lattice = cell$lattice )
    class( cell.new ) <- c( "vasp", "data.frame" )
    return( cell.new )
}
vasp.rotate.cell.default <- function( cell, lattice = NULL, angle = NULL ){
    if ( is.null( cell$x ) || is.null( cell$y ) )
        stop( "vasp.rotate.cell.default: supplied cell must have at least two named columns containing the x and y coordinates" )
    if ( is.null( lattice ) )
        stop( "vasp.rotate.cell.default: A data frame containing the lattice vectors has to be supplied!" )

    if ( is.null( angle ) ){
        angle <- acos( ( lattice[ 1, 1 ]* lattice[ 2, 1 ] + lattice[ 1, 2 ]* lattice[ 2, 2 ] )/
                          ( sqrt( lattice[ 1, 1 ]^2 + lattice[ 1, 2 ]^2 )*
                               sqrt( lattice[ 2, 1 ]^2 + lattice[ 2, 2 ]^2 ) ) )
        angle <- -angle/ 2
    }
    
    rotate.cell <- inline::cxxfunction( signature( cell = "numeric",
                                     angle = "double" ),
                          plugin = "Rcpp", body = '
Rcpp::DataFrame cellDataFrame(cell);
Rcpp::NumericVector xOld = cellDataFrame[ "x" ];
Rcpp::NumericVector yOld = cellDataFrame[ "y" ];
Rcpp::NumericVector xNew(xOld.size());
Rcpp::NumericVector yNew(yOld.size());
double phi = Rcpp::as<double>(angle);

for (int ii=0; ii<xOld.size(); ii++){
   xNew[ii] = xOld[ii]* cos(phi) - yOld[ii]* sin(phi);
   yNew[ii] = xOld[ii]* sin(phi) + yOld[ii]* cos(phi);
}

return Rcpp::DataFrame::create(Rcpp::Named("x")=xNew,
                               Rcpp::Named("y")=yNew);')
                              
    plane.new <- rotate.cell( cell, angle )
    cell$x <- plane.new$x
    cell$y <- plane.new$y
    return( cell )
}


##' @title vasp.bonds
##' @description Generates a data frame containing the start and the endpoints of all bounds in the system.
##'
##' @details Two atoms share a bond if they a closer than a specified distance. Since this function was originally written for a layered system the bonds are only created within one type of atoms. But its quite easily extendible to bonds between several kinds of atoms. Also compatible with object of class "vasp".
##'
##' @param atoms Data frame containing the x, y and z position of the atoms in cartesian coordinates as named columns and a fourth column names "type" containing a character specifying the atom species. See \code{\link{vasp.import}}$atoms.
##' @param distance Vector of length unique( atoms$type ) containing the bonding distances between the individual atoms. To establish a bond two atoms must be closer than this value. The order in this distances vector has to be the same as in unique( atoms$type ).
##'
##' @return Data frame with seven named columns "x.begin", "y.begin", "z.begin", "x.end", "y.end", "z.end" and "type" containing the start and end points of the bonds in the system as well as the type of atoms between which the bond is established. The later is necessary for the correct coloring and can be easily extended to bonds of many different atom types.
##' @family vasp
##' @author Philipp Mueller
vasp.bonds <- function( atoms, distance ){
    if ( any( class( atoms ) == "vasp" ) )
        atoms <- atoms$atoms
    if ( is.null( atoms$x ) || is.null( atoms$y ) || is.null( atoms$z ) )
        stop( "vasp.bonds: A data frame containing named columns 'x', 'y' and 'z' has to be provided in bonds." )
    if ( length( unique( atoms$type ) ) != length( distance ) )
        stop( "vasp.bonds: The number of bonding distance thresholds in bonds has to be of the same length than the number of types in the supplied data frame." )

    bond <- data.frame( x.begin = rep( -999, nrow( atoms )* 5 ),
                       y.begin = rep( -999, nrow( atoms )* 5 ),
                       z.begin = rep( -999, nrow( atoms )* 5 ),
                       x.end = rep( -999, nrow( atoms )* 5 ),
                       y.end = rep( -999, nrow( atoms )* 5 ),
                       z.end = rep( -999, nrow( atoms )* 5 ),
                       type = rep( "bla", nrow( atoms )* 5 ),
                       stringsAsFactors = FALSE )
    bond.count <- 1
    bond.type.count <- 0
    for ( tt in unique( atoms$type ) ){
        pos <- atoms[ atoms$type == tt, ]

        for ( ii in 1 : ( nrow( pos ) - 1 ) ){
            for ( jj in ( ii + 1 ) : nrow( pos ) ){
                if ( ( ( pos$x[[ ii ]] - pos$x[[ jj ]] )^2 +
                       ( pos$y[[ ii ]] - pos$y[[ jj ]] )^2 +
                       ( pos$z[[ ii ]] - pos$z[[ jj ]] )^2  ) <=
                    distance[ which( unique( atoms$type == tt ) ) ] ){
                    bond$x.begin[[ bond.count ]] <- pos$x[[ ii ]]
                    bond$y.begin[[ bond.count ]] <- pos$y[[ ii ]]
                    bond$z.begin[[ bond.count ]] <- pos$z[[ ii ]]
                    bond$x.end[[ bond.count ]] <- pos$x[[ jj ]]
                    bond$y.end[[ bond.count ]] <- pos$y[[ jj ]]
                    bond$z.end[[ bond.count ]] <- pos$z[[ jj ]]
                    bond.count <- bond.count + 1
                }
            }            
        }
        if ( bond.type.count !=  bond.count ){
            bond$type[ bond.type.count : ( bond.count - 1 ) ] <- tt
            bond.type.count <- bond.count
        }
    }
    bond.final <- bond[ bond$x.begin != -999, ]
    return( bond.final )
}

##' @title reproduce
##' @description Reproduce the content of a unit cell in multiple direction of its lattice vectors.
##'
##' @details The number of reproductions has to be provided as a sequence. E.g. seq(-3,3). Generic function working with both the specific content of the x, like the atom positions, and an object of class "vasp". 
##'
##' @param x Data frame containing at least two named columns "x" and "y" containing Cartesian coordinates or an object of class "vasp" to reproduce its atom positions and charge in one step.
##' @param lattice Data frame of dimension 3x3 containing the three lattice vectors in separate rows and their coordinates in different columns. See \code{\link{vasp.import}}. If an object of class "vasp" was supplied in x this argument is omitted.
##' @param x.rep Number of reproductions of the first lattice vector of the unit cell given in a sequence. 
##' @param y.rep Number of reproductions of the second lattice vector of the unit cell given in a sequence.
##' @param z.rep Number of reproductions of the third lattice vector of the unit cell given in a sequence.
##'
##' @family vasp
##' @return Input 'x' extended by the additional content of the reproduced unit cell.
##' @author Philipp Mueller
reproduce <- function( x, lattice = NULL, x.rep = NULL, y.rep = NULL, z.rep = NULL, x.window = NULL, y.window = NULL, z.window = NULL ){
    UseMethod( "reproduce" )
}
reproduce.vasp <- function( x, lattice = NULL, x.rep = NULL, y.rep = NULL, z.rep = NULL, x.window = NULL, y.window = NULL, z.window = NULL ){
    x$atoms <- reproduce.default( x$atoms, x$lattice, x.rep, y.rep, z.rep, x.window, y.window, z.window )
    x$charge <- reproduce.default( x$charge, x$lattice, x.rep, y.rep, z.rep, x.window, y.window, z.window )
    class( x ) <- c( "vasp", "data.frame" )
    return( x )
}    
reproduce.default <- function( x, lattice = NULL, x.rep = NULL, y.rep = NULL, z.rep = NULL, x.window = NULL, y.window = NULL, z.window = NULL ){
    if ( is.null( x$x ) || is.null( x$y ) || is.null( x$z ) )
        stop( "vasp.reproduce.default: A data frame containing named columns 'x', 'y' and 'z' has to be provided in bonds." )
    if ( is.null( lattice ) )
        stop( "vasp.reproduce.default: A data frame containing the lattice vectors has to be supplied!" )

    reproductions <- length( x.rep ) + length( y.rep ) + length( z.rep )
    if ( is.null( x.rep ) )
        x.rep <- 0
    if ( is.null( y.rep ) )
        y.rep <- 0
    if ( is.null( z.rep ) )
        z.rep <- 0
    
    output <- data.frame( x= rep( x$x, reproductions ), y = rep( x$y, reproductions ),
                         z = rep( x$z, reproductions ), something = rep( x[[ 4 ]], reproductions ) )
    names( output )[[ 4 ]] <- names( x )[[ 4 ]]

    ## Reproducing the x
    rc <- 1
    cl <- nrow( x )
    for ( ii in x.rep ){
        for ( jj in y.rep ){
            for ( kk in z.rep ){
                ## The original x is already included
                if ( ii == 0 && jj == 0 && kk == 0 )
                    next
                
                output[ ( cl* rc + 1 ) : ( cl*( rc + 1 ) ), 1 ] <- x$x + lattice[ 1, 1 ]* ii +
                    lattice[ 2, 1 ]* jj + lattice[ 3, 1 ]* kk
                output[ ( cl* rc + 1 ) : ( cl*( rc + 1 ) ), 2 ] <- x$y + lattice[ 1, 2 ]* ii +
                    lattice[ 2, 2 ]* jj + lattice[ 3, 2 ]* kk
                output[ ( cl* rc + 1 ) : ( cl*( rc + 1 ) ), 3 ] <- x$z + lattice[ 1, 3 ]* ii +
                    lattice[ 2, 3 ]* jj + lattice[ 3, 3 ]* kk
                output[ ( cl* rc + 1 ) : ( cl*( rc + 1 ) ), 4 ] <- x[[ 4 ]]
                rc <- rc + 1
            }
        }
    }

    ## Selecting a specific window of the reproduced x
    if ( ! is.null( x.window ) ){
        if ( length( x.window ) != 2 )
            error( "vasp.reproduce.default: provided window of the x axis seems not to be in the correct format" )
        if ( x.window[[ 1 ]] > x.window[[ 2 ]] )
            x.window <- rev( x.window )
        output <- output[ output$x > x.window[[ 1 ]] & output$x < x.window[[ 2 ]], ]
    }
    if ( ! is.null( y.window ) ){
        if ( length( y.window ) != 2 )
            error( "vasp.reproduce.default: provided window of the y ayis seems not to be in the correct format" )
        if ( y.window[[ 1 ]] > y.window[[ 2 ]] )
            y.window <- rev( y.window )
        output <- output[ output$y > y.window[[ 1 ]] & output$y < y.window[[ 2 ]], ]
    }
    if ( ! is.null( z.window ) ){
        if ( length( z.window ) != 2 )
            error( "vasp.reproduce.default: provided window of the z azis seems not to be in the correct format" )
        if ( z.window[[ 1 ]] > z.window[[ 2 ]] )
            z.window <- rev( z.window )
        output <- output[ output$z > z.window[[ 1 ]] & output$z < z.window[[ 2 ]], ]
    }    
    return( output )
}

##' @title diff.vasp
##' @description Calculates the charge difference of two charge densities.
##'
##' @details Additional method for the S3 generic base::diff. Only the charge of the second input is subtracted from the charge of the first input.
##'
##' @param x First object of class "vasp". Its used as the reference.
##' @param y Second object of class "vasp".
##'
##' @family vasp
##'
##' @return Input x with different values in the 'charge' element.
##' @author Philipp Mueller
diff.vasp <- function( x, y ){
    ## The restriction in the z grid of plot.charge.calc is removed.
    diff <- x
    diff$charge$charge <- x$charge$charge - y$charge$charge
    output <- list( charge = diff$charge, atoms = x$atoms, lattice = x$lattice )
    class( output ) <- c( "vasp", "list" )
    return( output )
}

plane <- function( x, height = mean( x$atoms[ x$atoms$type == unique( x$atoms$type )[ 1 ], ]$z, ... ) ){
    UseMethod( "plane" )
}
##' @title plane.vasp
##' @description Extracts a plane of charge rectangle to the z axis.
##'
##' @details Method of a S3 generic. The element containing all details about the atom positions is just copied to the output. Only the charge density is sliced.
##'
##' @param x Provided VASP result of class "vasp".
##' @param height Position on the z axis at which the slice through the charge density is going to be extracted. The default value is the mean z value of all atoms of the first type listed in x$atoms$type.
##'
##' @family vasp
##'
##' @return Object of class vasp containing just one z plane in its $charge element.
##' @author Philipp Mueller
plane.vasp <- function( x, height = mean( x$atoms[ x$atoms$type == unique( x$atoms$type )[ 1 ], ]$z ) ){
    ## Plotting the electron density at the z grid nearest to the Graphene layer.
    ## Where x is the latter one and y is the reference
    plot.height <- which( abs( unique( x$charge$z ) - height ) == 
                             min( abs( unique( x$charge$z ) - height ) ) )
    plot.data.diff <-  x$charge[ x$charge[ "z" ] == unique( x$charge$z )[ plot.height ], ]
    output <- list( charge = plot.data.diff, atoms = x$atoms )
    class( output ) <- c( "vasp", "data.frame" )
    return( output )
}
