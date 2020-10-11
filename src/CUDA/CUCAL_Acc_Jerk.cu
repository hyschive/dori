//---------------------------------------------------------------------------------------------------
// prefix "g" for pointer pointing to "Global" memory space in device
// prefix "s" for pointer pointing to "Shared" memory space in device
//
// prefix "I/J" for varialbe of "I/J" particles
// each thread calculates the acceleration, and jerk of "(I_Base+tx)th" "I" particle
// from all "J" particles
//
// "I" particles are packed into "N/GroupSize_I" "IGroup"s with "GroupSize_I" particles within each group
// "J" particles are packed into "N/GroupSize_J" "JGroup"s with "GroupSize_J" particles within each group
//---------------------------------------------------------------------------------------------------

#include <cutil.h>
#include "Dori.h"

#define I_Start         bx * BLOCK_SIZE
#define J_Start         0
#define GroupSize_I     GRID_SIZE * BLOCK_SIZE
#define GroupSize_J     BLOCK_SIZE

#define SDATA(array, index)  CUT_BANK_CHECKER(array, index)




__global__ void CUCAL_Acc_Jerk( const int Nj, real gJ_Mass[], real gJ_Pos[][3], real gJ_Vel[][3],
                                const int Ni, real gI_Pos[][3], real gI_Vel[][3],
                                real gI_Acc[][3], real gI_Jerk[][3], const real Eps2 )
{

   const unsigned int tx = threadIdx.x;
   const unsigned int bx = blockIdx.x;

   __shared__ real sJ_Mass [BLOCK_SIZE];

   __shared__ real sJ_Pos_x[BLOCK_SIZE];
   __shared__ real sJ_Pos_y[BLOCK_SIZE];
   __shared__ real sJ_Pos_z[BLOCK_SIZE];

   __shared__ real sJ_Vel_x[BLOCK_SIZE];
   __shared__ real sJ_Vel_y[BLOCK_SIZE];
   __shared__ real sJ_Vel_z[BLOCK_SIZE];



// (I/J)_Base : Base Address for (I/J)Group
   for (int I_Base=I_Start; I_Base<Ni; I_Base+=GroupSize_I)
   {
      real_acc Acc [3] = { (real_acc)0.0, (real_acc)0.0, (real_acc)0.0 };
      real_acc Jerk[3] = { (real_acc)0.0, (real_acc)0.0, (real_acc)0.0 };

      int ii = I_Base + tx;
      int i  = ii % Ni;

      real I_Pos_x = gI_Pos[i][0];
      real I_Pos_y = gI_Pos[i][1];
      real I_Pos_z = gI_Pos[i][2];

      real I_Vel_x = gI_Vel[i][0];
      real I_Vel_y = gI_Vel[i][1];
      real I_Vel_z = gI_Vel[i][2];

      for (int J_Base=J_Start; J_Base<Nj; J_Base+=GroupSize_J)
      {
         int jj = J_Base + tx;
         int j  = jj % Nj;

         SDATA(sJ_Mass,tx)  = gJ_Mass[j];

         SDATA(sJ_Pos_x,tx) = gJ_Pos [j][0];
         SDATA(sJ_Pos_y,tx) = gJ_Pos [j][1];
         SDATA(sJ_Pos_z,tx) = gJ_Pos [j][2];

         SDATA(sJ_Vel_x,tx) = gJ_Vel [j][0];
         SDATA(sJ_Vel_y,tx) = gJ_Vel [j][1];
         SDATA(sJ_Vel_z,tx) = gJ_Vel [j][2];

         __syncthreads();


//       k : kth particle in JGroup
         for (int k=0; k<GroupSize_J; k++)
         {
#           ifndef N_IS_MULTIPLE_OF_BS
            int kk = J_Base+k;
#           endif

//          evaluate the gravitational acceleration and jerk
//---------------------------------------------------------------------
            real dx = SDATA(sJ_Pos_x,k)-I_Pos_x;
            real dy = SDATA(sJ_Pos_y,k)-I_Pos_y;
            real dz = SDATA(sJ_Pos_z,k)-I_Pos_z;

#           ifdef SOFTEN
            real R2  = dx*dx + Eps2;
#           else
            real R2  = dx*dx;
#           endif
                 R2 += dy*dy;
                 R2 += dz*dz;

            real Rinv   = (real)1.0 / SQRT(R2);
            real R2inv  = Rinv*Rinv;
            real R3inv  = R2inv*Rinv;
            real MR3inv = SDATA(sJ_Mass,k)*R3inv;

            real dVx = SDATA(sJ_Vel_x,k)-I_Vel_x;
            real dVy = SDATA(sJ_Vel_y,k)-I_Vel_y;
            real dVz = SDATA(sJ_Vel_z,k)-I_Vel_z;

            real dR_dot_dV = dx*dVx + dy*dVy + dz*dVz;
            real Temp      = -(real)3.0*dR_dot_dV*R2inv;


#           ifndef N_IS_MULTIPLE_OF_BS
            if ( kk < Nj )
            {
#           endif

               Acc[0] += MR3inv*dx;
               Acc[1] += MR3inv*dy;
               Acc[2] += MR3inv*dz;

               Jerk[0] += MR3inv*( dVx + Temp*dx );
               Jerk[1] += MR3inv*( dVy + Temp*dy );
               Jerk[2] += MR3inv*( dVz + Temp*dz );

#           ifndef N_IS_MULTIPLE_OF_BS
            }
#           endif

         } // for (int k=0; k<GroupSize_J; k++)

         __syncthreads();
      } // for (int J_Base=J_Start; J_Base<Nj; J_Base+=GroupSize_J)

      if ( ii < Ni )
      {
         gI_Acc [i][0] = Acc [0];
         gI_Acc [i][1] = Acc [1];
         gI_Acc [i][2] = Acc [2];

         gI_Jerk[i][0] = Jerk[0];
         gI_Jerk[i][1] = Jerk[1];
         gI_Jerk[i][2] = Jerk[2];
      }

   } // for (int I_Base=I_Start; I_Base<Ni; I_Base+=GroupSize_I)

}

