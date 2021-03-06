
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
SEXP importCHGCAR( SEXP input) ;
}

// definition

SEXP importCHGCAR( SEXP input ){
BEGIN_RCPP

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
        
END_RCPP
}



