#ifndef __MACRO_H__
#define __MACRO_H__



// ****************************************************************************
// ** This header defines the symbolic constants and macros.                 **
// ** For clarity, useless options defined in the makefile will be "undef"   **
// ** in the end of this file.                                               **
// ****************************************************************************


// ########################
// ## Symbolic Constants ##
// ########################

// option == NONE --> the option is turned off
#define NONE      0


// GPU architecture
#define FERMI        1
#define KEPLER       2
#define MAXWELL      3
#define PASCAL       4
#define VOLTA        5
#define TURING       6


// extreme values
#ifndef __INT_MAX__
#  define __INT_MAX__      2147483647
#endif

#ifndef __LONG_MAX__
#  define __LONG_MAX__     9223372036854775807L
#endif

#ifndef __UINT_MAX__
#  define __UINT_MAX__     ( __INT_MAX__*2U + 1U )
#endif

#ifndef __FLT_MAX__
#  define __FLT_MAX__      3.40282347e+38F
#endif

#ifndef __FLT_MIN__
#  define __FLT_MIN__      1.17549435e-38F
#endif


// NULL values
#ifndef NULL
#  define NULL             0
#endif

#ifndef NULL_INT
#  define NULL_INT         __INT_MAX__
#endif

#ifndef NULL_REAL
#  define NULL_REAL        __FLT_MAX__
#endif

#ifndef NULL_BOOL
#  define NULL_BOOL        false
#endif


// macro for the function "Aux_Error"
#define ERROR_INFO         __FILE__, __LINE__, __FUNCTION__


// floating-point type for MPI
#ifdef FLOAT8
#  define GAMER_MPI_REAL   MPI_DOUBLE
#else
#  define GAMER_MPI_REAL   MPI_FLOAT
#endif



// ############
// ## Macros ##
// ############

// single/double-precision mathematic functions
#ifdef FLOAT8
#  define  SQRT( a )        sqrt( a )
#  define   POW( a, b )      pow( a, b )
#  define   SIN( a )         sin( a )
#  define  ATAN( a )        atan( a )
/*
#  define  FABS( a )        fabs( a )
#  define   COS( a )         cos( a )
#  define   LOG( a )         log( a )
#  define  FMAX( a, b )     fmax( a, b )
#  define  FMIN( a, b )     fmin( a, b )
#  define  FMOD( a, b )     fmod( a, b )
#  define ATAN2( a, b )    atan2( a, b )
*/
#else
#  define  SQRT( a )        sqrtf( a )
#  define   POW( a, b )      powf( a, b )
#  define   SIN( a )         sinf( a )
#  define  ATAN( a )        atanf( a )
/*
#  define  FABS( a )        fabsf( a )
#  define   COS( a )         cosf( a )
#  define   LOG( a )         logf( a )
#  define  FMAX( a, b )     fmaxf( a, b )
#  define  FMIN( a, b )     fminf( a, b )
#  define  FMOD( a, b )     fmodf( a, b )
#  define ATAN2( a, b )    atan2f( a, b )
*/
#endif


// sign function
#define SIGN( a )       (  ( (a) < (real)0.0 ) ? (real)-1.0 : (real)+1.0  )


// max/min functions
#define MAX( a, b )     (  ( (a) > (b) ) ? (a) : (b)  )
#define MIN( a, b )     (  ( (a) < (b) ) ? (a) : (b)  )


// square/cube function
#define SQR(  a )       ( (a)*(a)     )
#define CUBE( a )       ( (a)*(a)*(a) )



// ################################
// ## Remove useless definitions ##
// ################################



#endif  // #ifndef __MACRO_H__
