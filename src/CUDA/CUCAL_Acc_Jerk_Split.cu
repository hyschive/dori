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


// Variables Definition:
// NSplit             : number of times to split
// Segment            : a segment of threads containing all I particles and ceil to a multiple of block size
// Pseudo-Particle    : the pseudo IDs of I-particles when ceiling the size of segment to a multiple of block size
//                    : eg. Ni_perSeg = 128, Ni = 126 --> Pseudo-Particle Ids = 126,127
// Ni_perSeg          : number of I particles per segment
//                      eg. BLOCK_SIZE=128, Ni=129~256 --> Ni_perSeg = 128*2 = 256
//                      eg. BLOCK_SIZE=64,  Ni=129~192 --> Ni_perSeg =  64*3 = 192
// Ni_allSeg          : total number of I particles in all segments (including the pseudo-particles)
// NBlock_perSeg      : number of blocks in each segment
// NthSeg             : the segment being calculated (0, 1, ... NSplit-1)
// Nj_afterSplit_List : number of J-particles calculated by each segment
//---------------------------------------------------------------------------------------------------


#include "Dori.h"

#define I_Start         __umul24(bx, BLOCK_SIZE)
#define J_Start         __umul24(NthSeg, Nj_afterSplit_List[0])
#define I_End           Ni_allSeg
#define J_End           ( J_Start + Nj_afterSplit_List[NthSeg] )
#define GroupSize_I     GRID_SIZE * BLOCK_SIZE
#define GroupSize_J     BLOCK_SIZE




__global__ void CUCAL_Acc_Jerk_Split( const int Nj, real gJ_Mass[], real gJ_Pos[][3],
                                      real gJ_Vel[][3], const int Ni, real gI_Pos[][3], real gI_Vel[][3],
                                      real gI_Acc[][3], real gI_Jerk[][3], const real Eps2,
                                      const unsigned int NSplit, const unsigned int NBlock_perSeg,
                                      const unsigned int Ni_perSeg, const unsigned int Ni_allSeg,
                                      const unsigned int Nj_afterSplit_List[] )
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


   int I_Iter = 0;

// (I/J)_Base : Base Address for (I/J)Group
   for ( int I_Base=I_Start; I_Base<I_End; I_Base+=GroupSize_I, I_Iter++ )
   {
      real_acc Acc [3] = { (real_acc)0.0, (real_acc)0.0, (real_acc)0.0 };
      real_acc Jerk[3] = { (real_acc)0.0, (real_acc)0.0, (real_acc)0.0 };

      int NthSeg = (bx+I_Iter*GRID_SIZE) / NBlock_perSeg;

      int I_Seq  = I_Base+tx;
      int ii     = I_Seq%Ni_perSeg;
      int i      = ii%Ni;  // prevent from loading the pseudo-I-particles

      real I_Pos_x = gI_Pos[i][0];
      real I_Pos_y = gI_Pos[i][1];
      real I_Pos_z = gI_Pos[i][2];

      real I_Vel_x = gI_Vel[i][0];
      real I_Vel_y = gI_Vel[i][1];
      real I_Vel_z = gI_Vel[i][2];

      for (int J_Base=J_Start; J_Base<J_End; J_Base+=GroupSize_J)
      {
         int jj = J_Base+tx;
         int j  = jj%Nj;   // deal with the case that Nj is not a multiple of block size

         sJ_Mass [tx] = gJ_Mass[j];

         sJ_Pos_x[tx] = gJ_Pos [j][0];
         sJ_Pos_y[tx] = gJ_Pos [j][1];
         sJ_Pos_z[tx] = gJ_Pos [j][2];

         sJ_Vel_x[tx] = gJ_Vel [j][0];
         sJ_Vel_y[tx] = gJ_Vel [j][1];
         sJ_Vel_z[tx] = gJ_Vel [j][2];

         __syncthreads();

//       k : kth particle in JGroup
         for (int k=0; k<GroupSize_J; k++)
         {
#           ifndef N_IS_MULTIPLE_OF_BS
            int kk = J_Base+k;
#           endif


//          evaluate the gravitational acceleration and jerk
//---------------------------------------------------------------------
            real dx  = sJ_Pos_x[k] - I_Pos_x;
            real dy  = sJ_Pos_y[k] - I_Pos_y;
            real dz  = sJ_Pos_z[k] - I_Pos_z;

#           ifdef SOFTEN
            real R2  = dx*dx + Eps2;
#           else
            real R2  = dx*dx;
#           endif
                 R2 += dy*dy;
                 R2 += dz*dz;

            real Rinv      = (real)1.0 / SQRT(R2);
            real R2inv     = Rinv*Rinv;
            real R3inv     = R2inv*Rinv;
            real MR3inv    = sJ_Mass[k]*R3inv;

            real dVx       = sJ_Vel_x[k] - I_Vel_x;
            real dVy       = sJ_Vel_y[k] - I_Vel_y;
            real dVz       = sJ_Vel_z[k] - I_Vel_z;

            real dR_dot_dV = dx*dVx + dy*dVy + dz*dVz;
            real Temp      = -(real)3.0*dR_dot_dV*R2inv;

#           ifndef N_IS_MULTIPLE_OF_BS
            if ( kk < J_End )
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

      } // for (int J_Base=J_Start; J_Base<J_End; J_Base+=GroupSize_J)

      if ( ii < Ni )
      {
         const unsigned int SaveIndex = ii + NthSeg*Ni;

         gI_Acc [SaveIndex][0] = Acc [0];
         gI_Acc [SaveIndex][1] = Acc [1];
         gI_Acc [SaveIndex][2] = Acc [2];

         gI_Jerk[SaveIndex][0] = Jerk[0];
         gI_Jerk[SaveIndex][1] = Jerk[1];
         gI_Jerk[SaveIndex][2] = Jerk[2];
      }

   } // for ( int I_Base=I_Start; I_Base<I_End; I_Base+=GroupSize_I, I_Iter++ )

}

