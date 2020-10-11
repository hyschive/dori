#include "Dori.h"




//---------------------------------------------------------------------------------------------------
// Function    :  CPU_Pot
// Description :  Use CPU to calculate the gravitational potential of "I" particles from all
//                "J" particles ( J --> I )
//
// Note        :  Prefixes "I/J" for variables of "I/J" particles
//---------------------------------------------------------------------------------------------------
void CPU_Pot( const int Nj, real J_Mass[], real J_Pos[][3], const int Ni, real I_Pos[][3], real I_Pot[],
              const real Eps2, const real G )
{

   real_acc Pot;     // only the accumulation arrays are declared as "double" for the option "FLOAT8_ACC"
   real dr[3], r2;

   for (int i=0; i<Ni; i++)
   {

//    initialize the accumulation array as zero
      Pot = (real_acc)0.0;

      for (int j=0; j<Nj; j++)
      {
//       evaluate the gravitational potential
         for (int d=0; d<3; d++)    dr[d] = J_Pos[j][d] - I_Pos[i][d];

         r2   = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
#        ifdef SOFTEN
         r2  += Eps2;
#        endif

//       exclude contribution from itself
#        ifdef SOFTEN
         if ( r2 != Eps2 )
#        else
         if ( r2 != (real)0.0 )
#        endif
         Pot += -J_Mass[j] / SQRT(r2);
      } // for (int j=0; j<Nj; j++)

//    store the calculation results back to the input array
      I_Pot[i] = Pot;

    } // for (int i=0; i<Ni; i++)

// multiply results by the gravitational constant if it's not unity
   if ( G != (real)1.0 )
   {
      for (int i=0; i<Ni; i++)   I_Pot[i] *= G;
   }

}



