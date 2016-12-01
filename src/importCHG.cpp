
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
SEXP importCHG( SEXP input) ;
}

// definition

SEXP importCHG( SEXP input ){
BEGIN_RCPP

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
        
END_RCPP
}



