#include "Dori.h"


real M_ext = 1.0;
real R_ext = 1.0;



//----------------------------------------------------------------------
// Function    :  Ext_AddAccFromFunc
// Description :  Add external acceleration from the analytical function
//
// Note        :  1. External acceleration (and jerk) will be ADDED to the input array
//                2. External jerk is calculated only when EXT_JERK == true
//
// Parameter   :  NPar   : Number of particles to be calculated
//                MyPos  : Position array
//                MyVel  : Velocity array
//                MyAcc  : Acceleration array
//                MyJerk : Jerk array
//                Time   : Target physical time
//----------------------------------------------------------------------
void Ext_AddAccFromFunc( const int NPar, const real (*MyPos)[3], const real (*MyVel)[3],
                         real (*MyAcc)[3], real (*MyJerk)[3], const double Time )
{

// ===================================================================
// set the external gravity to represent the second star in a binary system
// ===================================================================
   const real Theta0 = M_PI;
   const real W      = sqrt( 0.25*NEWTON_G*M_ext/R_ext )/R_ext;

   real dv[3], dr[3], r, GM_r3, Temp;
   real r_ext[3], v_ext[3], Theta;

   Theta    = W*Time + Theta0;
   r_ext[0] = R_ext*cos( Theta );
   r_ext[1] = R_ext*sin( Theta );
   r_ext[2] = 0.0;
   v_ext[0] = -R_ext*W*sin( Theta );
   v_ext[1] = +R_ext*W*cos( Theta );
   v_ext[2] = 0.0;

#  pragma omp parallel for private( dv, dr, r, GM_r3, Temp ) schedule( runtime )
   for (int p=0; p<NPar; p++)
   {
      for (int d=0; d<3; d++)
      {
         dr[d] = r_ext[d] - MyPos[p][d];
         dv[d] = v_ext[d] - MyVel[p][d];
      }

      r     = SQRT( dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2] );
      GM_r3 = NEWTON_G*M_ext/(r*r*r);

      for (int d=0; d<3; d++)    MyAcc[p][d] += GM_r3*dr[d];

      if ( EXT_JERK )
      {
         Temp = -(real)3.0*( dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2] )/(r*r);

         for (int d=0; d<3; d++)    MyJerk[p][d] += GM_r3*( dv[d] + Temp*dr[d] );
      }
   } // for (int p=0; p<NPar; p++)
// ===================================================================

} // FUNCTION : Ext_AddAccFromFunc



//----------------------------------------------------------------------
// Function    :  Ext_AddAccFromFile
// Description :  Add external acceleration from the loaded file
//
// Note        :  1. External acceleration (and jerk) will be ADDED to the input array
//                2. External jerk is calculated only when EXT_JERK == true
//
// Parameter   :  NPar   : Number of particles to be calculated
//                MyPos  : Position array
//                MyVel  : Velocity array
//                MyAcc  : Acceleration array
//                MyJerk : Jerk array
//                Time   : Target physical time
//----------------------------------------------------------------------
void Ext_AddAccFromFile( const int NPar, const real (*MyPos)[3], const real (*MyVel)[3],
                         real (*MyAcc)[3], real (*MyJerk)[3], const double Time )
{

   const int  NSkip = (EXT_JERK) ? ( (EXT_ACC_DER==EXT_ACC_DER_QUAR) ? 2 : 1 ) : 0;    // we don't have dAcc at the boundary cells
   const real _dh   = (real)1.0 / EXT_DH;
   const real ONE   = (real)1.0;


// 1. acceleration : CIC scheme
   if ( EXT_PAR_INT == EXT_PAR_INT_CIC )
   {
#     pragma omp parallel
      {
//       1D array -> 3D array (somehow they must be put inside the OpenMP parallel section ...)
         real (* AccX)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccX;
         real (* AccY)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccY;
         real (* AccZ)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccZ;

         real (*dAccX)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccX;
         real (*dAccY)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccY;
         real (*dAccZ)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccZ;

         real (* Acc[3])[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = {  AccX,  AccY,  AccZ };
         real (*dAcc[3])[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = { dAccX, dAccY, dAccZ };


         real MyPosInBox[3];              // particle position relative to the array corner (in the unit of cell)
         real dr[3];                      // distance to the left cell (in the unit of cell)
         long idxL[3], idxR[3];           // array index of the left/right cell


#        pragma omp for schedule( runtime )
         for (long p=0; p<NPar; p++)
         {
//          get the nearest cell indices
            for (int d=0; d<3; d++)
            {
               MyPosInBox[d]  = ( MyPos[p][d] + EXT_CEN[d] )*_dh;
               dr        [d]  = MyPosInBox[d] - (real)0.5;
               idxL      [d]  = (int)dr[d];
               idxR      [d]  = idxL[d] + 1;
               dr        [d] -= (real)idxL[d];
            } // for (int d=0; d<3; d++)


//          skip particles lying outside the input array
            if (  idxL[0] < NSkip  ||  idxL[1] < NSkip  ||  idxL[2] < NSkip  ||
                  idxR[0] >= EXT_SIZE[0]-NSkip  ||  idxR[1] >= EXT_SIZE[1]-NSkip  ||  idxR[2] >= EXT_SIZE[2]-NSkip  )
            {
               Aux_Message( stderr, "particle position (%14.7e,%14.7e,%14.7e) lies outside the external acceleration array !!\n",
                            MyPos[p][0], MyPos[p][1], MyPos[p][2] );
               continue;
            }


//          add the external acceleration
            for (int d=0; d<3; d++)
            {
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxL[1] ][ idxL[0] ] * ( ONE-dr[2] ) * ( ONE-dr[1] ) * ( ONE-dr[0] );
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxL[1] ][ idxR[0] ] * ( ONE-dr[2] ) * ( ONE-dr[1] ) * (     dr[0] );
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxR[1] ][ idxL[0] ] * ( ONE-dr[2] ) * (     dr[1] ) * ( ONE-dr[0] );
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxL[1] ][ idxL[0] ] * (     dr[2] ) * ( ONE-dr[1] ) * ( ONE-dr[0] );
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxR[1] ][ idxR[0] ] * ( ONE-dr[2] ) * (     dr[1] ) * (     dr[0] );
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxR[1] ][ idxL[0] ] * (     dr[2] ) * (     dr[1] ) * ( ONE-dr[0] );
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxL[1] ][ idxR[0] ] * (     dr[2] ) * ( ONE-dr[1] ) * (     dr[0] );
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxR[1] ][ idxR[0] ] * (     dr[2] ) * (     dr[1] ) * (     dr[0] );
            }


//          add the external jerk
            if ( EXT_JERK )
            {
               for (int j=0; j<3; j++)                                  // [0/1/2] = [jerk_x/y/z]
               {
                  real der[3] = { (real)0.0, (real)0.0, (real)0.0 };    // [0/1/2] = [dx/dy/dz]

                  for (int d=0; d<3; d++)
                  {
                     der[d] += dAcc[j][ idxL[2] ][ idxL[1] ][ idxL[0] ][d] * ( ONE-dr[2] ) * ( ONE-dr[1] ) * ( ONE-dr[0] );
                     der[d] += dAcc[j][ idxL[2] ][ idxL[1] ][ idxR[0] ][d] * ( ONE-dr[2] ) * ( ONE-dr[1] ) * (     dr[0] );
                     der[d] += dAcc[j][ idxL[2] ][ idxR[1] ][ idxL[0] ][d] * ( ONE-dr[2] ) * (     dr[1] ) * ( ONE-dr[0] );
                     der[d] += dAcc[j][ idxR[2] ][ idxL[1] ][ idxL[0] ][d] * (     dr[2] ) * ( ONE-dr[1] ) * ( ONE-dr[0] );
                     der[d] += dAcc[j][ idxL[2] ][ idxR[1] ][ idxR[0] ][d] * ( ONE-dr[2] ) * (     dr[1] ) * (     dr[0] );
                     der[d] += dAcc[j][ idxR[2] ][ idxR[1] ][ idxL[0] ][d] * (     dr[2] ) * (     dr[1] ) * ( ONE-dr[0] );
                     der[d] += dAcc[j][ idxR[2] ][ idxL[1] ][ idxR[0] ][d] * (     dr[2] ) * ( ONE-dr[1] ) * (     dr[0] );
                     der[d] += dAcc[j][ idxR[2] ][ idxR[1] ][ idxR[0] ][d] * (     dr[2] ) * (     dr[1] ) * (     dr[0] );
                  }

                  MyJerk[p][j] += der[0]*MyVel[p][0] + der[1]*MyVel[p][1] + der[2]*MyVel[p][2];

               } // for (int j=0; j<3; j++)
            } // if ( EXT_JERK )
         } // for (long p=0; p<NPar; p++)
      } // pragma omp parallel
   } // if ( EXT_PAR_INT == EXT_PAR_INT_CIC )


// 2. acceleration : TSC scheme
   else if ( EXT_PAR_INT == EXT_PAR_INT_TSC )
   {
#     pragma omp parallel
      {
//       1D array -> 3D array (somehow they must be put inside the OpenMP parallel section ...)
         real (* AccX)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccX;
         real (* AccY)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccY;
         real (* AccZ)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccZ;

         real (*dAccX)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccX;
         real (*dAccY)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccY;
         real (*dAccZ)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccZ;

         real (* Acc[3])[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = {  AccX,  AccY,  AccZ };
         real (*dAcc[3])[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = { dAccX, dAccY, dAccZ };


         real MyPosInBox[3];              // particle position relative to the array corner (in the unit of cell)
         real dr[3];                      // distance to the center of the central cell (can be negative)
         long idxL[3], idxC[3], idxR[3];  // array index of the left/central/right cell


#        pragma omp for schedule( runtime )
         for (long p=0; p<NPar; p++)
         {
//          get the nearby cell indices
            for (int d=0; d<3; d++)
            {
               MyPosInBox[d] = ( MyPos[p][d] + EXT_CEN[d] )*_dh;
               idxC      [d] = (int)MyPosInBox[d];
               idxL      [d] = idxC[d] - 1;
               idxR      [d] = idxC[d] + 1;
               dr        [d] = MyPosInBox[d] - (real)idxC[d] - (real)0.5;
            }


//          skip particles lying outside the input array
            if (  idxL[0] < NSkip  ||  idxL[1] < NSkip  ||  idxL[2] < NSkip  ||
                  idxR[0] >= EXT_SIZE[0]-NSkip  ||  idxR[1] >= EXT_SIZE[1]-NSkip  ||  idxR[2] >= EXT_SIZE[2]-NSkip  )
            {
               Aux_Message( stderr, "particle position (%14.7e,%14.7e,%14.7e) lies outside the external acceleration array !!\n",
                            MyPos[p][0], MyPos[p][1], MyPos[p][2] );
               continue;
            }


//          define helpful macros
#           define WL(x) ( (real)0.5*((real)0.5-(x))*((real)0.5-(x)) )
#           define WR(x) ( (real)0.5*((real)0.5+(x))*((real)0.5+(x)) )
#           define WC(x) ( (real)0.75-(x)*(x) )


//          add the external acceleration
            for (int d=0; d<3; d++)
            {
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxL[1] ][ idxL[0] ] * WL(dr[2]) * WL(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxL[1] ][ idxC[0] ] * WL(dr[2]) * WL(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxL[1] ][ idxR[0] ] * WL(dr[2]) * WL(dr[1]) * WR(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxC[1] ][ idxL[0] ] * WL(dr[2]) * WC(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxC[1] ][ idxC[0] ] * WL(dr[2]) * WC(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxC[1] ][ idxR[0] ] * WL(dr[2]) * WC(dr[1]) * WR(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxR[1] ][ idxL[0] ] * WL(dr[2]) * WR(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxR[1] ][ idxC[0] ] * WL(dr[2]) * WR(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxL[2] ][ idxR[1] ][ idxR[0] ] * WL(dr[2]) * WR(dr[1]) * WR(dr[0]);

               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxL[1] ][ idxL[0] ] * WC(dr[2]) * WL(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxL[1] ][ idxC[0] ] * WC(dr[2]) * WL(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxL[1] ][ idxR[0] ] * WC(dr[2]) * WL(dr[1]) * WR(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxC[1] ][ idxL[0] ] * WC(dr[2]) * WC(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxC[1] ][ idxC[0] ] * WC(dr[2]) * WC(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxC[1] ][ idxR[0] ] * WC(dr[2]) * WC(dr[1]) * WR(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxR[1] ][ idxL[0] ] * WC(dr[2]) * WR(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxR[1] ][ idxC[0] ] * WC(dr[2]) * WR(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxC[2] ][ idxR[1] ][ idxR[0] ] * WC(dr[2]) * WR(dr[1]) * WR(dr[0]);

               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxL[1] ][ idxL[0] ] * WR(dr[2]) * WL(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxL[1] ][ idxC[0] ] * WR(dr[2]) * WL(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxL[1] ][ idxR[0] ] * WR(dr[2]) * WL(dr[1]) * WR(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxC[1] ][ idxL[0] ] * WR(dr[2]) * WC(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxC[1] ][ idxC[0] ] * WR(dr[2]) * WC(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxC[1] ][ idxR[0] ] * WR(dr[2]) * WC(dr[1]) * WR(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxR[1] ][ idxL[0] ] * WR(dr[2]) * WR(dr[1]) * WL(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxR[1] ][ idxC[0] ] * WR(dr[2]) * WR(dr[1]) * WC(dr[0]);
               MyAcc[p][d] += Acc[d][ idxR[2] ][ idxR[1] ][ idxR[0] ] * WR(dr[2]) * WR(dr[1]) * WR(dr[0]);
            } // for (int d=0; d<3; d++)


//          add the external jerk
            if ( EXT_JERK )
            {
               for (int j=0; j<3; j++)                                  // [0/1/2] = [jerk_x/y/z]
               {
                  real der[3] = { (real)0.0, (real)0.0, (real)0.0 };    // [0/1/2] = [dx/dy/dz]

                  for (int d=0; d<3; d++)
                  {
                     der[d] += dAcc[j][ idxL[2] ][ idxL[1] ][ idxL[0] ][d] * WL(dr[2]) * WL(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxL[1] ][ idxC[0] ][d] * WL(dr[2]) * WL(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxL[1] ][ idxR[0] ][d] * WL(dr[2]) * WL(dr[1]) * WR(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxC[1] ][ idxL[0] ][d] * WL(dr[2]) * WC(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxC[1] ][ idxC[0] ][d] * WL(dr[2]) * WC(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxC[1] ][ idxR[0] ][d] * WL(dr[2]) * WC(dr[1]) * WR(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxR[1] ][ idxL[0] ][d] * WL(dr[2]) * WR(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxR[1] ][ idxC[0] ][d] * WL(dr[2]) * WR(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxL[2] ][ idxR[1] ][ idxR[0] ][d] * WL(dr[2]) * WR(dr[1]) * WR(dr[0]);

                     der[d] += dAcc[j][ idxC[2] ][ idxL[1] ][ idxL[0] ][d] * WC(dr[2]) * WL(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxL[1] ][ idxC[0] ][d] * WC(dr[2]) * WL(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxL[1] ][ idxR[0] ][d] * WC(dr[2]) * WL(dr[1]) * WR(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxC[1] ][ idxL[0] ][d] * WC(dr[2]) * WC(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxC[1] ][ idxC[0] ][d] * WC(dr[2]) * WC(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxC[1] ][ idxR[0] ][d] * WC(dr[2]) * WC(dr[1]) * WR(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxR[1] ][ idxL[0] ][d] * WC(dr[2]) * WR(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxR[1] ][ idxC[0] ][d] * WC(dr[2]) * WR(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxC[2] ][ idxR[1] ][ idxR[0] ][d] * WC(dr[2]) * WR(dr[1]) * WR(dr[0]);

                     der[d] += dAcc[j][ idxR[2] ][ idxL[1] ][ idxL[0] ][d] * WR(dr[2]) * WL(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxL[1] ][ idxC[0] ][d] * WR(dr[2]) * WL(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxL[1] ][ idxR[0] ][d] * WR(dr[2]) * WL(dr[1]) * WR(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxC[1] ][ idxL[0] ][d] * WR(dr[2]) * WC(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxC[1] ][ idxC[0] ][d] * WR(dr[2]) * WC(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxC[1] ][ idxR[0] ][d] * WR(dr[2]) * WC(dr[1]) * WR(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxR[1] ][ idxL[0] ][d] * WR(dr[2]) * WR(dr[1]) * WL(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxR[1] ][ idxC[0] ][d] * WR(dr[2]) * WR(dr[1]) * WC(dr[0]);
                     der[d] += dAcc[j][ idxR[2] ][ idxR[1] ][ idxR[0] ][d] * WR(dr[2]) * WR(dr[1]) * WR(dr[0]);
                  } // for (int d=0; d<3; d++)

                  MyJerk[p][j] += der[0]*MyVel[p][0] + der[1]*MyVel[p][1] + der[2]*MyVel[p][2];

               } // for (int j=0; j<3; j++)
            } // if ( EXT_JERK )
         } // for (long p=0; p<NPar; p++)
      } // pragma omp parallel
   } // else if ( EXT_PAR_INT == EXT_PAR_INT_TSC )


   else
      Aux_Error( ERROR_INFO, "unsupported EXT_PAR_INT (%d) !!\n", EXT_PAR_INT );

} // FUNCTION : Ext_AddAccFromFile



//----------------------------------------------------------------------
// Function    :  Ext_LoadExtAcc
// Description :  Load external acceleration from files
//
// Note        :  1. External acceleration files are always named "EXT_ACCX/Y/Z"
//                2. External jerk will also be calculated when EXT_JERK == true
//
// Parameter   :
//----------------------------------------------------------------------
void Ext_LoadExtAcc()
{

   if ( MyRank == 0 )   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. load the external acceleration
   if ( MyRank == 0 )   Aux_Message( stdout, "   Load external acceleration ..." );


   const char *FileName[3] = { "EXT_ACCX", "EXT_ACCY", "EXT_ACCZ" };
   const long  ExpectSize  = EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2]*sizeof(real);

   long InputSize;

   for (int TRank=0; TRank<NGPU; TRank++)
   {
      if ( MyRank == TRank )
      {
         for (int d=0; d<3; d++)
         {
//          check file existence
            if ( MyRank == 0 )
            if ( !Aux_CheckFileExist(FileName[d]) )   Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName[d] );

            FILE *File_Acc = fopen( FileName[d], "rb" );

//          check file size
            if ( MyRank == 0 )
            {
               fseek( File_Acc, 0, SEEK_END );
               InputSize = ftell( File_Acc );

               if ( InputSize != ExpectSize )
                  Aux_Error( ERROR_INFO, "the size of the file <%s> is incorrect --> input = %ld <-> expect = %ld !!\n",
                             FileName[d], InputSize, ExpectSize );

               rewind( File_Acc );
            }

//          load data
            if      ( d == 0 )   fread( Ext_AccX, sizeof(real), (long)EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2], File_Acc );
            else if ( d == 1 )   fread( Ext_AccY, sizeof(real), (long)EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2], File_Acc );
            else if ( d == 2 )   fread( Ext_AccZ, sizeof(real), (long)EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2], File_Acc );

            fclose( File_Acc );
         } // for (int d=0; d<3; d++)
      } // if ( MyRank == TRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TRank=0; TRank<NGPU; TRank++)


   if ( MyRank == 0 )   Aux_Message( stdout, "done\n" );


// 2. calculate the derivatives
   if ( EXT_JERK )
   {
      if ( MyRank == 0 )   Aux_Message( stdout, "   Calculate the spatial derivatives of the external acceleration ..." );

//    initialize as zero (meaningful only for the outermost cells)
      for (long t=0; t<(long)EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2]*3; t++)
      {
         Ext_dAccX[t] = (real)0.0;
         Ext_dAccY[t] = (real)0.0;
         Ext_dAccZ[t] = (real)0.0;
      }


//    1D array -> 3D array
      real (* AccX)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccX;
      real (* AccY)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccY;
      real (* AccZ)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    )Ext_AccZ;

      real (*dAccX)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccX;
      real (*dAccY)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccY;
      real (*dAccZ)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = ( real(*)[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] )Ext_dAccZ;

      real (* Acc[3])[ EXT_SIZE[1] ][ EXT_SIZE[0] ]    = {  AccX,  AccY,  AccZ };
      real (*dAcc[3])[ EXT_SIZE[1] ][ EXT_SIZE[0] ][3] = { dAccX, dAccY, dAccZ };


//    calculate the x/y/z derivatives
      if ( EXT_ACC_DER == EXT_ACC_DER_QUAD )    // quadratic interpolation
      {
         const real _2dh = 0.5 / EXT_DH;
         int im, ip, jm, jp, km, kp;

         for (int d=0; d<3; d++)                {
         for (int k=1; k<EXT_SIZE[2]-1; k++)    {  km = k - 1;    kp = k + 1;
         for (int j=1; j<EXT_SIZE[1]-1; j++)    {  jm = j - 1;    jp = j + 1;
         for (int i=1; i<EXT_SIZE[0]-1; i++)    {  im = i - 1;    ip = i + 1;

            dAcc[d][k][j][i][0] = _2dh * ( Acc[d][k ][j ][ip] - Acc[d][k ][j ][im] );  // x derivative
            dAcc[d][k][j][i][1] = _2dh * ( Acc[d][k ][jp][i ] - Acc[d][k ][jm][i ] );  // y derivative
            dAcc[d][k][j][i][2] = _2dh * ( Acc[d][kp][j ][i ] - Acc[d][km][j ][i ] );  // z derivative

         }}}}
      } // if ( EXT_ACC_DER == EXT_ACC_DER_QUAD )


      else if ( EXT_ACC_DER == EXT_ACC_DER_QUAR )     // quartic interpolation
      {
         const real _12dh = 1.0 / (12.0*EXT_DH);
         int im1, ip1, jm1, jp1, km1, kp1;
         int im2, ip2, jm2, jp2, km2, kp2;

         for (int d=0; d<3; d++)                {
         for (int k=2; k<EXT_SIZE[2]-2; k++)    {  km1 = k - 1;   kp1 = k + 1;   km2 = k - 2;   kp2 = k + 2;
         for (int j=2; j<EXT_SIZE[1]-2; j++)    {  jm1 = j - 1;   jp1 = j + 1;   jm2 = j - 2;   jp2 = j + 2;
         for (int i=2; i<EXT_SIZE[0]-2; i++)    {  im1 = i - 1;   ip1 = i + 1;   im2 = i - 2;   ip2 = i + 2;

//          x derivative
            dAcc[d][k][j][i][0] = _12dh * ( - 1.0*Acc[d][k  ][j  ][ip2] + 1.0*Acc[d][k  ][j  ][im2]
                                            + 8.0*Acc[d][k  ][j  ][ip1] - 8.0*Acc[d][k  ][j  ][im1] );

//          y derivative
            dAcc[d][k][j][i][1] = _12dh * ( - 1.0*Acc[d][k  ][jp2][i  ] + 1.0*Acc[d][k  ][jm2][i  ]
                                            + 8.0*Acc[d][k  ][jp1][i  ] - 8.0*Acc[d][k  ][jm1][i  ] );

//          z derivative
            dAcc[d][k][j][i][2] = _12dh * ( - 1.0*Acc[d][kp2][j  ][i  ] + 1.0*Acc[d][km2][j  ][i  ]
                                            + 8.0*Acc[d][kp1][j  ][i  ] - 8.0*Acc[d][km1][j  ][i  ] );

         }}}}
      } // else if ( EXT_ACC_DER == EXT_ACC_DER_QUAR )


      else
         Aux_Error( ERROR_INFO, "unsupported EXT_ACC_DER (%d) !!\n", EXT_ACC_DER );


      if ( MyRank == 0 )   Aux_Message( stdout, "done\n" );
   } // if ( EXT_JERK )

} // FUNCTION : Ext_LoadExtAcc

