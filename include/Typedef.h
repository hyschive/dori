#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__



// ****************************************************************************
// ** This header defines the "typedef and enum"                             **
// ** Please DO NOT modify the number assigned to any enumerator constant !! **
// ****************************************************************************


// single/double precision for the accumulation arrays
#ifdef FLOAT8_ACC
typedef double real_acc;
#else
typedef float  real_acc;
#endif


// single/double precision for all floating-point variables
#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif


// program initialization
enum OptInitMethod_t { INIT_FILE=1, INIT_FUNC=2 };

// gravity types
enum OptGravityType_t { GRAVITY_SELF=1, GRAVITY_EXTERNAL=2, GRAVITY_BOTH=3 };

// external acceleration initialization
enum OptExtMethod_t { EXT_FILE=1, EXT_FUNC=2 };

// scheme for interpolating external acceleration/jerk to the particles' position
enum OptExtParInt_t { EXT_PAR_INT_CIC = 1, EXT_PAR_INT_TSC = 2 };

// scheme for calculating external acceleration derivatives
enum OptExtAccDer_t { EXT_ACC_DER_QUAD = 1, EXT_ACC_DER_QUAR = 2 };



#endif  // #ifndef __TYPEDEF_H__
