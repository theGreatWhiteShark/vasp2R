##' @title vasp.import
##' @description Imports the CHGCAR file of a VASP calculation into R. Should also run with the CHG file.
##'
##' @details The CHGCAR file got a very odd format which has to be reformatted to ensure the its compatibility with R. In the CHGCAR file the number of grid points of the specific axes are given in the 12 column. In the charge data the x coordinate is the fastest and the z the slowest varying index. The charge density is also divided by the volume of the cell inside of this script (VASP internal stuff). The individual lattice vectors in the are given in different rows and their x, y and z are provided in the different columns.
##'
##' @param path Character vector containing the path to the CHGCAR file.
##'
##' @family vasp
##' 
##' @return List of three named data frames:
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
    if ( grep( "CHGCAR", path ) ){
        ## CHGCAR was supplied
        charge.lines <- charge.entries/ 5
    } else {
        charge.lines <- charge.entries/ 10
    }

    ## Import the raw format of the charge grid.
    charge.raw <- readr::read_table( path, skip = ( 10 + sum( lattice.numbers ) ), n_max= charge.lines, col_names = FALSE )/ lattice.volume
    if ( grep( "CHGCAR", path ) ){
        names( charge.raw ) <- c( "V1", "V2", "V3", "V4", "V5" )
    } else {
        ## This ensures the functionality of the function even when the CHG file is provided.
        names( charge.raw ) <- c( "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10" )
    }
    ## The main difference between the CHGCAR and the CHG file is the format. The charge density is stored with higher precision in the CHGCAR and with 5 columns per row. The CHG contains 10 columns per row.
    if ( grep( "CHGCAR", path ) ){
        concatC <- '
            Rcpp::DataFrame dfin( input );
            Rcpp::NumericVector v1 = dfin["V1"];
            Rcpp::NumericVector v2 = dfin["V2"];
            Rcpp::NumericVector v3 = dfin["V3"];
            Rcpp::NumericVector v4 = dfin["V4"];
            Rcpp::NumericVector v5 = dfin["V5"];
        
            int dfoutLength = v1.size()* 5;
            Rcpp::NumericVector dfout( dfoutLength );

            for ( int ii = 0; ii < v1.size(); ii++ ){
                dfout[ ii* 5 ] = v1[ ii ];
                dfout[ ( ii* 5 ) + 1 ] = v2[ ii ];
                dfout[ ( ii* 5 ) + 2 ] = v3[ ii ];
                dfout[ ( ii* 5 ) + 3 ] = v4[ ii ];
                dfout[ ( ii* 5 ) + 4 ] = v5[ ii ];
            }
            return dfout;
        '
    } else {
        concatC <- '
            Rcpp::DataFrame dfin( input );
            Rcpp::NumericVector v1 = dfin["V1"];
            Rcpp::NumericVector v2 = dfin["V2"];
            Rcpp::NumericVector v3 = dfin["V3"];
            Rcpp::NumericVector v4 = dfin["V4"];
            Rcpp::NumericVector v5 = dfin["V5"];
            Rcpp::NumericVector v6 = dfin["V6"];
            Rcpp::NumericVector v7 = dfin["V7"];
            Rcpp::NumericVector v8 = dfin["V8"];
            Rcpp::NumericVector v9 = dfin["V9"];
            Rcpp::NumericVector v10 = dfin["V10"];
        
            int dfoutLength = v1.size()* 10;
            Rcpp::NumericVector dfout( dfoutLength );

            for ( int ii = 0; ii < v1.size(); ii++ ){
                dfout[ ii* 10 ] = v1[ ii ];
                dfout[ ( ii* 10 ) + 1 ] = v2[ ii ];
                dfout[ ( ii* 10 ) + 2 ] = v3[ ii ];
                dfout[ ( ii* 10 ) + 3 ] = v4[ ii ];
                dfout[ ( ii* 10 ) + 4 ] = v5[ ii ];
                dfout[ ( ii* 10 ) + 5 ] = v6[ ii ];
                dfout[ ( ii* 10 ) + 6 ] = v7[ ii ];
                dfout[ ( ii* 10 ) + 7 ] = v8[ ii ];
                dfout[ ( ii* 10 ) + 8 ] = v9[ ii ];
                dfout[ ( ii* 10 ) + 9 ] = v10[ ii ];
            }
            return dfout;
        '
    }
    ## After this step there is just one format regardless of the input file.
    ## The rows of the charge density will be concatenated one by one. So A[a,b;c,d] becomes V[a,b,c,d]. 
    concatCharge <- inline::cxxfunction( signature( input =  "data.frame" ), concatC, plugin = "Rcpp" )
    charge.vector <- concatCharge( input = ( charge.raw ) )

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
    formatCharge <- inline::cxxfunction( signature( chargeVector = "numeric",
                                        xGridVector = "numeric",
                                        yGridVector = "numeric",
                                        zGridVector = "numeric",
                                        latticeVector = "numeric" ),
                             plugin = "Rcpp",
                             body = '    
    Rcpp::NumericVector charge(chargeVector);
    Rcpp::NumericVector xPos(charge.size());
    Rcpp::NumericVector yPos(charge.size());
    Rcpp::NumericVector zPos(charge.size());
    Rcpp::NumericVector xGrid(xGridVector);
    Rcpp::NumericVector yGrid(yGridVector);
    Rcpp::NumericVector zGrid(zGridVector);
    Rcpp::DataFrame latt(latticeVector);
    // The x, y and z components of the lattice vectors are given in vectors and not the lattice vectors themselves!
    Rcpp::NumericVector lattX = latt["x"];
    Rcpp::NumericVector lattY = latt["y"];
    Rcpp::NumericVector lattZ = latt["z"];
    int xCount = 0;
    int yCount = 0;
    int zCount = 0;

    for (int ii=0; ii<charge.size(); ii++) {
        xPos[ii] = lattX[0]*xGrid[xCount] + lattX[1]*yGrid[yCount] + lattX[2]*zGrid[zCount];
        yPos[ii] = lattY[0]*xGrid[xCount] + lattY[1]*yGrid[yCount] + lattY[2]*zGrid[zCount];
        zPos[ii] = lattZ[0]*xGrid[xCount] + lattZ[1]*yGrid[yCount] + lattZ[2]*zGrid[zCount];

        xCount++;
        if ( xCount == xGrid.size() ) {
            xCount = 0;
            yCount++;
        }
        if ( yCount == yGrid.size() ) {
            yCount = 0;
            zCount++;
        }
    }

    return Rcpp::DataFrame::create(Rcpp::Named("x")=xPos,
                                   Rcpp::Named("y")=yPos,
                                   Rcpp::Named("z")=zPos,
                                   Rcpp::Named("charge")=charge);'
                             )
    charge.df <- formatCharge( charge.vector, x.range, y.range, z.range, lattice.vectors)
    
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
