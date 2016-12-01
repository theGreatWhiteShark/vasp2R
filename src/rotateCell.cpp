
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
SEXP rotateCell( SEXP cell, SEXP angle) ;
}

// definition

SEXP rotateCell( SEXP cell, SEXP angle ){
BEGIN_RCPP

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
                               Rcpp::Named("y")=yNew);
END_RCPP
}



