#include "Dori.h"




//---------------------------------------------------------------------------------------------------
// Function    :  CPU_Acc_Jerk
// Description :  Use CPU to calculate the gravitational acceleration and jerk of "I" particles from
//                all "J" particles ( J --> I )
//
// Note        :  Prefixes "I/J" for variables of "I/J" particles
//---------------------------------------------------------------------------------------------------
void CPU_Acc_Jerk( const int Nj, const real J_Mass[], const real J_Pos[][3], const real J_Vel[][3],
                   const int Ni, const real I_Pos[][3], const real I_Vel[][3],
                   real I_Acc[][3], real I_Jerk[][3], const real Eps2, const real G )
{

   real_acc Acc[3], Jerk[3];  // only the accumulation arrays are declared as "double" for the option "FLOAT8_ACC"
   real dr[3], dv[3], r2, _r, _r2, _r3, M_r3, Temp;

   for (int i=0; i<Ni; i++)
   {

//    initialize the accumulation arrays as zero
      for (int d=0; d<3; d++)
      {
         Acc [d] = (real_acc)0.0;
         Jerk[d] = (real_acc)0.0;
      }

      for (int j=0; j<Nj; j++)
      {
//       evaluate the acceleration and jerk
         for (int d=0; d<3; d++)    dr[d] = J_Pos[j][d] - I_Pos[i][d];

         r2   = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
#        ifdef SOFTEN
         r2  += Eps2;
#        endif
         _r   = 1.0 / SQRT(r2);
         _r2  = _r*_r;
         _r3  = _r2*_r;
         M_r3 = J_Mass[j]*_r3;

         for (int d=0; d<3; d++)    dv[d] = J_Vel[j][d] - I_Vel[i][d];

         Temp = -(real)3.0*( dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2] )*_r2;

         for (int d=0; d<3; d++)
         {
            Acc [d] += M_r3*dr[d];
            Jerk[d] += M_r3*( dv[d] + Temp*dr[d] );
         }
      } // for (int j=0; j<Nj; j++)

//    store the calculation results back to the input arrays
      for (int d=0; d<3; d++)
      {
         I_Acc [i][d] = Acc [d];
         I_Jerk[i][d] = Jerk[d];
      }
   } // for (int i=0; i<Ni; i++)

// multiply results by the gravitational constant if it's not unity
   if ( G != (real)1.0 )
   {
      for (int i=0; i<Ni; i++)
      for (int d=0; d<3; d++)
      {
         I_Acc [i][d] *= G;
         I_Jerk[i][d] *= G;
      }
   }

}



