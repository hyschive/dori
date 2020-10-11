//---------------------------------------------------------------------------------------------------
// prefix "g" for pointer pointing to "Global" memory space in device
// prefix "s" for pointer pointing to "Shared" memory space in device
//
// prefix "I/J" for varialbe of "I/J" particles
// each thread calculates the potential of "(I_Base+tx)th" "I" particle
// from all "J" particles
//
// "I" particles are packed into "N/GroupSize_I" "IGroup"s with "GroupSize_I" particles within each group
// "J" particles are packed into "N/GroupSize_J" "JGroup"s with "GroupSize_J" particles within each group
//---------------------------------------------------------------------------------------------------

#include "Dori.h"

#define I_Start      bx * BLOCK_SIZE
#define J_Start      0
#define GroupSize_I  GRID_SIZE * BLOCK_SIZE
#define GroupSize_J  BLOCK_SIZE




__global__ void CUCAL_Pot( const int Nj, real gJ_Mass[], real gJ_Pos[][3],
                           const int Ni, real gI_Pos[][3], real gI_Pot[], const real Eps2 )
{

   const unsigned int tx = threadIdx.x;
   const unsigned int bx = blockIdx.x;

   __shared__ real sJ_Mass [BLOCK_SIZE];

   __shared__ real sJ_Pos_x[BLOCK_SIZE];
   __shared__ real sJ_Pos_y[BLOCK_SIZE];
   __shared__ real sJ_Pos_z[BLOCK_SIZE];


// (I/J)_Base : Base Address for (I/J)Group
   for (int I_Base=I_Start; I_Base<Ni; I_Base+=GroupSize_I)
   {
      real_acc Pot = (real_acc)0.0;

      int ii = I_Base+tx;
      int i  = ii%Ni;

      real I_Pos_x = gI_Pos[i][0];
      real I_Pos_y = gI_Pos[i][1];
      real I_Pos_z = gI_Pos[i][2];

      for (int J_Base=J_Start; J_Base<Nj; J_Base+=GroupSize_J)
      {
          int jj = J_Base+tx;
          int j  = jj%Nj;

          sJ_Mass [tx] = gJ_Mass[j];

          sJ_Pos_x[tx] = gJ_Pos [j][0];
          sJ_Pos_y[tx] = gJ_Pos [j][1];
          sJ_Pos_z[tx] = gJ_Pos [j][2];

          __syncthreads();

//       k : kth particle in JGroup
         for (int k=0; k<GroupSize_J; k++)
         {
#           ifndef N_IS_MULTIPLE_OF_BS
            int kk = J_Base+k;
#           endif


//          evaluate the gravitaional potential
//---------------------------------------------------------------------
            real dx = sJ_Pos_x[k] - I_Pos_x;
            real dy = sJ_Pos_y[k] - I_Pos_y;
            real dz = sJ_Pos_z[k] - I_Pos_z;

#           ifdef SOFTEN
            real R2  = dx*dx + Eps2;
#           else
            real R2  = dx*dx;
#           endif
                 R2 += dy*dy;
                 R2 += dz*dz;

            real mRinv = -(real)1.0 / SQRT(R2);

#           ifndef N_IS_MULTIPLE_OF_BS
            if ( kk < Nj )
            {
#           endif

//             exclude contribution from itself
#              ifdef SOFTEN
               if ( R2 != Eps2 )
#              else
               if ( R2 != (real)0.0 )
               #endif
               Pot += sJ_Mass[k]*mRinv;

#           ifndef N_IS_MULTIPLE_OF_BS
            }
#           endif

         } // for (int k=0; k<GroupSize_J; k++)

         __syncthreads();

      } // for (int J_Base=J_Start; J_Base<Nj; J_Base+=GroupSize_J)


      if ( ii < Ni )
      {
         gI_Pot [i] = Pot;
      }

   } // for (int I_Base=I_Start; I_Base<Ni; I_Base+=GroupSize_I)

}

