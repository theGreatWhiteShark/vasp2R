
// includes from the plugin

#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;


// user includes


// declarations
extern "C" {
SEXP formatCharge( SEXP chargeVector, SEXP xGridVector, SEXP yGridVector, SEXP zGridVector, SEXP latticeVector) ;
}

// definition

SEXP formatCharge( SEXP chargeVector, SEXP xGridVector, SEXP yGridVector, SEXP zGridVector, SEXP latticeVector ){
BEGIN_RCPP
    
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
                                   Rcpp::Named("charge")=charge);
END_RCPP
}



