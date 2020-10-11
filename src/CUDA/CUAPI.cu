#include "Dori.h"

__global__ void CUCAL_Acc_Jerk( const int Nj, real gJ_Mass[], real gJ_Pos[][3], real gJ_Vel[][3],
                                const int Ni, real gI_Pos[][3], real gI_Vel[][3],
                                real gI_Acc[][3], real gI_Jerk[][3], const real Eps2 );

__global__ void CUCAL_Acc_Jerk_Split( const int Nj, real gJ_Mass[], real gJ_Pos[][3],
                                      real gJ_Vel[][3], const int Ni, real gI_Pos[][3], real gI_Vel[][3],
                                      real gI_Acc[][3], real gI_Jerk[][3], const real Eps2,
                                      const unsigned int NSplit, const unsigned int NBlock_perSeg,
                                      const unsigned int Ni_perSeg, const unsigned int Ni_allSeg,
                                      const unsigned int Nj_afterSplit_List[] );

__global__ void CUCAL_Pot( const int Nj, real gJ_Mass[], real gJ_Pos[][3],
                           const int Ni, real gI_Pos[][3], real gI_Pot[], const real Eps2 );


static const dim3 Block_Dim(BLOCK_SIZE, 1, 1);
static const dim3 Grid_Dim(GRID_SIZE, 1, 1);

// declare all device pointers
static real  *dJ_Mass                     = NULL;
static real (*dJ_Pos )[3]                 = NULL;
static real (*dJ_Vel )[3]                 = NULL;
static real (*dI_Pos )[3]                 = NULL;
static real (*dI_Vel )[3]                 = NULL;
static real (*dI_Acc )[3]                 = NULL;
static real (*dI_Jerk)[3]                 = NULL;
static real  *dI_Pot                      = NULL;
static unsigned int *d_Nj_afterSplit_List = NULL;

// variables storing the splitting results in CUAPI_Get_Acc_Jerk_Split function
static real (*hI_Acc_Split)  [3]          = NULL;
static real (*hI_Jerk_Split) [3]          = NULL;
static unsigned int *h_Nj_afterSplit_List = NULL;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetDevice
// Description :  Set the active device
//
// Parameter   :  Mode :    -2 --> set automatically by CUDA (must work with the "compute-exclusive mode")
//                          -1 --> set by MPI ranks : SetDeviceID = MyRank % DeviceCount
//                       >=  0 --> set to "Mode"
//-------------------------------------------------------------------------------------------------------
void CUAPI_SetDevice( const int Mode )
{

   if ( MyRank == 0 )   printf( "Set GPU ... \n" );


// check
   if ( Mode < -2 )
   {
      fprintf( stderr, "ERROR : \"incorrect parameter %s = %d\" !!\n", "Mode", Mode );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      exit(-1);
   }


// get the hostname of each MPI process
   char Host[1024];
   gethostname( Host, 1024 );


// verify that there are GPU supporing CUDA
   int DeviceCount;
   CUDA_SAFE_CALL(  cudaGetDeviceCount( &DeviceCount )  );

   if ( DeviceCount == 0 )
   {
      fprintf( stderr, "ERROR : no devices supporting CUDA at rank %2d (host = %8s) !!\n", MyRank, Host );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      exit(-1);
   }


// set the device ID
   void **d_TempPtr = NULL;
   int SetDeviceID, GetDeviceID = 999;
   cudaDeviceProp DeviceProp;

   switch ( Mode )
   {
      case -2:
         CUDA_SAFE_CALL(  cudaMalloc( (void**) &d_TempPtr, sizeof(int) )  );  // to set the GPU ID
         CUDA_SAFE_CALL(  cudaFree( d_TempPtr )  );

//       make sure that the "exclusive" compute mode is adopted
         CUDA_SAFE_CALL(  cudaGetDevice( &GetDeviceID )  );
         CUDA_SAFE_CALL(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );

         if ( DeviceProp.computeMode != cudaComputeModeExclusive )
         {
            fprintf( stderr, "WARNING : the \"exclusive\" compute mode is NOT enabled for \"%s\" ",
                     "GPUID_SELECT = -2" );
            fprintf( stderr,           "at rank %2d (host = %8s) !!\n", MyRank, Host );
         }
         break;

      case -1:
         SetDeviceID = MyRank % DeviceCount;
         CUDA_SAFE_CALL(  cudaSetDevice( SetDeviceID )  );

         if ( NGPU > 1  &&  MyRank == 0 )
         {
            fprintf( stderr, "WARNING : please make sure that different MPI ranks will use different GPUs for " );
            fprintf( stderr,           "\"%s\" !!\n", "GPUID_SELECT = -1" );
         }
         break;

      default:
         SetDeviceID = Mode;

         if ( SetDeviceID < DeviceCount )
         {
            CUDA_SAFE_CALL(  cudaSetDevice( SetDeviceID )  );
         }

         else
         {
            fprintf( stderr, "ERROR : SetDeviceID (%d) >= DeviceCount (%d) at rank %2d (host = %8s) !!\n",
                     SetDeviceID, DeviceCount, MyRank, Host );
            fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,
                     __FUNCTION__  );
            exit(-1);
         }

         if ( NGPU > 1  &&  MyRank == 0 )
         {
            fprintf( stderr, "WARNING : please make sure that different MPI ranks will use different GPUs for " );
            fprintf( stderr,           "\"%s\" !!\n", "GPUID_SELECT >= 0" );
         }
         break;
   }


// check
   GetDeviceID = -999;
   CUDA_SAFE_CALL(  cudaGetDevice( &GetDeviceID )  );
   CUDA_SAFE_CALL(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );

// check
// (1) verify the device version
   if ( DeviceProp.major < 1 )
   {
      fprintf( stderr, "\n" );
      fprintf( stderr, "ERROR : device major version < 1 at rank %2d (host = %8s) !!\n", MyRank, Host );
      exit(-1);
   }

// (2) verify that the device ID is properly set
   if ( Mode != -2  &&  GetDeviceID != SetDeviceID )
   {
      fprintf( stderr, "ERROR : GetDeviceID (%d) != SetDeviceID (%d) at rank %2d (host = %8s) !!\n",
               GetDeviceID, SetDeviceID, MyRank, Host );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      exit(-1);
   }

// (3) verify that the adopted ID is accessible
   CUDA_SAFE_CALL(  cudaMalloc( (void**) &d_TempPtr, sizeof(int) )  );
   CUDA_SAFE_CALL(  cudaFree( d_TempPtr )  );

// (4) verify the capability of double precision
#  ifdef FLOAT8_ACC
   if ( DeviceProp.major < 2  &&  DeviceProp.minor < 3 )
   {
      fprintf( stderr, "ERROR : the GPU \"%s\" at rank %2d (host = %8s) does not support double precision !!\n",
               DeviceProp.name, MyRank, Host );
      fprintf( stderr, "        Please turn off the options \"FLOAT8\" and \"FLOAT8_ACC\" in the Makefile.\n" );
      exit(-1);
   }
#  endif

// (5) verify the GPU architecture
#  if   ( GPU_ARCH == FERMI )
   if ( DeviceProp.major != 2 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Fermi architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == KEPLER )
   if ( DeviceProp.major != 3 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Kepler architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == MAXWELL )
   if ( DeviceProp.major != 5 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Maxwell architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == PASCAL )
   if ( DeviceProp.major != 6 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Pascal architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == VOLTA )
   if ( DeviceProp.major != 7  &&  DeviceProp.minor != 0 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Volta architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  elif ( GPU_ARCH == TURING )
   if ( DeviceProp.major != 7  &&  DeviceProp.minor != 5 )
      Aux_Error( ERROR_INFO, "GPU \"%s\" with the compute capability %d.%d is incompatible with the Turing architecture !!\n"
                             "        --> Please reset GPU_ARCH in the Makefile properly\n",
                 DeviceProp.name, DeviceProp.major, DeviceProp.minor );

#  else
#  error : UNKNOWN GPU_ARCH !!
#  endif // GPU_ARCH

// (6) check if GRID_SIZE is a multiple of the number of multiprocesors
   if ( GRID_SIZE % DeviceProp.multiProcessorCount != 0 )
   {
      fprintf( stderr, "WARNING : GRID_SIZE (%d) is not a multiple of the number of multiprocessors (%d) !!\n",
               GRID_SIZE, DeviceProp.multiProcessorCount );
      fprintf( stderr, "          Rank %2d, host = %8s\n", MyRank, Host );
      fprintf( stderr, "          --> The performance may be deteriorated ...\n" );
   }

// (7) check if BLOCK_SIZE is a multiple of the warp size
   if ( BLOCK_SIZE % DeviceProp.warpSize != 0 )
   {
      fprintf( stderr, "WARNING : BLOCK_SIZE (%d) is not a multiple of the warp size (%d) !!\n",
               BLOCK_SIZE, DeviceProp.warpSize );
      fprintf( stderr, "          Rank %2d, host = %8s\n", MyRank, Host );
      fprintf( stderr, "          --> The performance may be deteriorated ...\n" );
   }


   if ( MyRank == 0 )   printf( "Set GPU ... done\n" );

}



//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_DiagnoseDevice
// Description :  Take a diagnosis of each GPU
//-------------------------------------------------------------------------------------------------------
void CUAPI_DiagnoseDevice()
{

   if ( MyRank == 0 )    printf( "Diagnose GPU ...\n" );


// get the hostname and PID of each process
   const int PID = getpid();
   char Host[1024];
   gethostname( Host, 1024 );


// get the number of devices
   int DeviceCount;
   CUDA_SAFE_CALL(  cudaGetDeviceCount( &DeviceCount )  );

   if ( DeviceCount == 0 )
   {
      fprintf( stderr, "ERROR : no devices supporting CUDA at rank %2d (host = %8s) !!\n", MyRank, Host );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      exit(-1);
   }


// get the device ID
   int GetDeviceID = 999;
   CUDA_SAFE_CALL(  cudaGetDevice( &GetDeviceID )  );


// load the device properties
   cudaDeviceProp DeviceProp;
   CUDA_SAFE_CALL(  cudaGetDeviceProperties( &DeviceProp, GetDeviceID )  );


// record the device properties
   const char FileName[] = "Record__Note";

   if ( MyRank == 0 )
   {
       FILE *Note = fopen( FileName, "a" );
       fprintf( Note, "Device Diagnosis\n" );
       fprintf( Note, "***********************************************************************************\n" );
       fclose( Note );
   }

   for (int YourTurn=0; YourTurn<NGPU; YourTurn++)
   {
      if ( MyRank == YourTurn )
      {
         int DriverVersion = 0, RuntimeVersion = 0;
         CUDA_SAFE_CALL(  cudaDriverGetVersion( &DriverVersion )  );
         CUDA_SAFE_CALL(  cudaRuntimeGetVersion( &RuntimeVersion )  );

         FILE *Note = fopen( FileName, "a" );
         fprintf( Note, "MPI rank = %3d, hostname = %10s, PID = %5d\n\n", MyRank, Host, PID );
         fflush( Note );

         GetCPUInfo( FileName );

         fprintf( Note, "\n" );
         fprintf( Note, "Number of GPUs                    : %d\n"    , DeviceCount );
         fprintf( Note, "GPU ID                            : %d\n"    , GetDeviceID );
         fprintf( Note, "GPU Name                          : %s\n"    , DeviceProp.name );
         fprintf( Note, "CUDA Driver Version               : %d.%d\n" , DriverVersion/1000, DriverVersion%100 );
         fprintf( Note, "CUDA Runtime Version              : %d.%d\n" , RuntimeVersion/1000, RuntimeVersion%100 );
         fprintf( Note, "CUDA Major Revision Number        : %d\n"    , DeviceProp.major );
         fprintf( Note, "CUDA Minor Revision Number        : %d\n"    , DeviceProp.minor );
         fprintf( Note, "Clock Rate                        : %f GHz\n", DeviceProp.clockRate/1.0e6);
         fprintf( Note, "Global Memory Size                : %d MB\n" , DeviceProp.totalGlobalMem/1024/1024 );
         fprintf( Note, "Constant Memory Size              : %d KB\n" , DeviceProp.totalConstMem/1024 );
         fprintf( Note, "Shared Memory Size per Block      : %d KB\n" , DeviceProp.sharedMemPerBlock/1024 );
         fprintf( Note, "Number of Registers per Block     : %d\n"    , DeviceProp.regsPerBlock );
         fprintf( Note, "Warp Size                         : %d\n"    , DeviceProp.warpSize );
         fprintf( Note, "Number of Multiprocessors:        : %d\n"    , DeviceProp.multiProcessorCount );
         fprintf( Note, "Number of Cores:                  : %d\n"    , DeviceProp.multiProcessorCount *
                                                                        (DeviceProp.major == 2 ? 32 : 8) );
         fprintf( Note, "Max Number of Threads per Block   : %d\n"    , DeviceProp.maxThreadsPerBlock );
         fprintf( Note, "Max Size of the Block X-Dimension : %d\n"    , DeviceProp.maxThreadsDim[0] );
         fprintf( Note, "Max Size of the Grid X-Dimension  : %d\n"    , DeviceProp.maxGridSize[0] );
         fprintf( Note, "Concurrent Copy and Execution     : %s\n"    , DeviceProp.deviceOverlap ? "Yes" : "No" );
#        if ( CUDART_VERSION >= 3000 )
         fprintf( Note, "Concurrent Kernel Execution       : %s\n"    , DeviceProp.concurrentKernels ? "Yes" :
                                                                                                       "No" );
#        endif
#        if ( CUDART_VERSION >= 3010 )
         fprintf( Note, "GPU has ECC Support Enabled       : %s\n"    , DeviceProp.ECCEnabled ? "Yes" : "No" );
#        endif

         fprintf( Note, "\n\n" );

         fclose( Note );
       }

       MPI_Barrier( MPI_COMM_WORLD );

   } // for (int YourTurn=0; YourTurn<NGPU; YourTurn++)

   if ( MyRank == 0 )
   {
      FILE *Note = fopen( FileName, "a" );
      fprintf( Note, "***********************************************************************************\n" );
      fclose( Note );

      printf( "Diagnose GPU ... done\n" );
   }

}



//-----------------------------------------------------------------------------------------
// Function    :  CUAPI_Get_Pot
// Description :  Calculate the potential of "I" particles from all "J" particles
//
// Note        :  Prefix "d" for pointers pointing to the "Device" memory space
//                Prefix "h" for pointers pointing to the "Host" memory space
//                Prefix "I/J" for varialbes of "I/J" particles
//-----------------------------------------------------------------------------------------
void CUAPI_Pot( const int Nj, const real hJ_Mass[], const real hJ_Pos[][3],
                const int Ni, const real hI_Pos[][3], real hI_Pot[], const real Eps2, const real G )
{

   const unsigned int J_Mem_Size = Nj * sizeof(real);
   const unsigned int I_Mem_Size = Ni * sizeof(real);


// copy data from host to device
//=========================================================================================
   CUDA_SAFE_CALL( cudaMemcpy( dJ_Mass, hJ_Mass,  J_Mem_Size, cudaMemcpyHostToDevice ) );
   CUDA_SAFE_CALL( cudaMemcpy( dJ_Pos,  hJ_Pos, 3*J_Mem_Size, cudaMemcpyHostToDevice ) );

   CUDA_SAFE_CALL( cudaMemcpy( dI_Pos,  hI_Pos, 3*I_Mem_Size, cudaMemcpyHostToDevice ) );


// execute the kernel to get potential
//=========================================================================================
   CUCAL_Pot <<< Grid_Dim, Block_Dim >>> ( Nj, dJ_Mass, dJ_Pos, Ni, dI_Pos, dI_Pot, Eps2 );
   CUT_CHECK_ERROR( "Kernel execution failed" );


// copy data from device to host
//=========================================================================================
   CUDA_SAFE_CALL( cudaMemcpy( hI_Pot, dI_Pot, I_Mem_Size, cudaMemcpyDeviceToHost ) );


// multiply results by the gravitational constant if it's not unity
//=========================================================================================
   if ( G != (real)1.0 )
   {
      for (int i=0; i<Ni; i++)   hI_Pot[i] *= G;
   }

}



//-----------------------------------------------------------------------------------------
// Function    :  CUAPI_Acc_Jerk
// Description :  Use GPU to calculate the acceleration and jerk of "I" particles from all "J" particles
//
//                prefix "d" for pointer pointing to "Device" memory space
//                prefix "h" for pointer pointing to "Host" memory space
//
//                prefix "I/J" for variables of "I/J" particles
//-----------------------------------------------------------------------------------------
void CUAPI_Acc_Jerk( const int Nj, const real hJ_Mass[], const real hJ_Pos[][3], const real hJ_Vel[][3],
                     const int Ni, const real hI_Pos[][3], const real hI_Vel[][3],
                     real hI_Acc[][3], real hI_Jerk[][3], const real Eps2, const bool Copy_J_MassPosVel,
                     const real G )
{

// determine the number of times for splitting
   unsigned int NSplit;

   if ( Ni > SPLIT_NMAX )              NSplit = 1;
   else if ( Ni > 0.5   *SPLIT_NMAX )  NSplit = 2;
   else if ( Ni > 0.25  *SPLIT_NMAX )  NSplit = 4;
   else if ( Ni > 0.125 *SPLIT_NMAX )  NSplit = 8;
   else if ( Ni > 0.0625*SPLIT_NMAX )  NSplit = 16;
   else                                NSplit = 32;

#  ifdef N_IS_MULTIPLE_OF_BS
   const int Max_NSplit = Nj / BLOCK_SIZE;
   if (NSplit > Max_NSplit)   NSplit = Max_NSplit;
#  endif

// enforce NSplit < Nj
   while ( NSplit > Nj )   NSplit /= 2;


   const unsigned int NBlock_perSeg = (Ni-1)/BLOCK_SIZE + 1;
   const unsigned int Ni_perSeg     = NBlock_perSeg*BLOCK_SIZE;
   const unsigned int Ni_allSeg     = Ni_perSeg*NSplit;
   const unsigned int Ni_afterSplit = Ni * NSplit;
   const unsigned int Nj_afterSplit = Nj / NSplit;

   for (int i=0; i<NSplit-1; i++)  h_Nj_afterSplit_List[i] = Nj_afterSplit;
   h_Nj_afterSplit_List[NSplit-1] = Nj - Nj_afterSplit*(NSplit-1);

   const unsigned int Out_Mem_Size = 3 * Ni_afterSplit * sizeof(real);
   const unsigned int J_Mem_Size   = Nj * sizeof(real);
   const unsigned int I_Mem_Size   = Ni * sizeof(real);


//  copy data from host to device
//=========================================================================================
   if ( Copy_J_MassPosVel )
   {
      CUDA_SAFE_CALL( cudaMemcpy( dJ_Mass, hJ_Mass,  J_Mem_Size, cudaMemcpyHostToDevice ) );
      CUDA_SAFE_CALL( cudaMemcpy( dJ_Pos,  hJ_Pos, 3*J_Mem_Size, cudaMemcpyHostToDevice ) );
      CUDA_SAFE_CALL( cudaMemcpy( dJ_Vel,  hJ_Vel, 3*J_Mem_Size, cudaMemcpyHostToDevice ) );
   }

      CUDA_SAFE_CALL( cudaMemcpy( dI_Pos,  hI_Pos, 3*I_Mem_Size, cudaMemcpyHostToDevice ) );
      CUDA_SAFE_CALL( cudaMemcpy( dI_Vel,  hI_Vel, 3*I_Mem_Size, cudaMemcpyHostToDevice ) );

      CUDA_SAFE_CALL( cudaMemcpy( d_Nj_afterSplit_List, h_Nj_afterSplit_List, 32*sizeof(int),
                                  cudaMemcpyHostToDevice ) );


// execute the kernel to get acceleration, and jerk
//=========================================================================================
   if ( NSplit == 1 )
   {
      CUCAL_Acc_Jerk <<< Grid_Dim, Block_Dim >>> ( Nj, dJ_Mass, dJ_Pos, dJ_Vel,
      Ni, dI_Pos, dI_Vel, dI_Acc, dI_Jerk, Eps2 );
   }
   else
   {
      CUCAL_Acc_Jerk_Split <<< Grid_Dim, Block_Dim >>> ( Nj, dJ_Mass, dJ_Pos, dJ_Vel,
                                                         Ni, dI_Pos, dI_Vel, dI_Acc, dI_Jerk, Eps2,
                                                         NSplit, NBlock_perSeg, Ni_perSeg, Ni_allSeg,
                                                         d_Nj_afterSplit_List );
   }
   CUT_CHECK_ERROR( "Kernel execution failed" );


// copy data from device to host
//=========================================================================================
   if ( NSplit == 1 )
   {
      CUDA_SAFE_CALL( cudaMemcpy( hI_Acc,  dI_Acc,  Out_Mem_Size, cudaMemcpyDeviceToHost ) );
      CUDA_SAFE_CALL( cudaMemcpy( hI_Jerk, dI_Jerk, Out_Mem_Size, cudaMemcpyDeviceToHost ) );
   }
   else
   {
      CUDA_SAFE_CALL( cudaMemcpy( hI_Acc_Split,  dI_Acc,  Out_Mem_Size, cudaMemcpyDeviceToHost ) );
      CUDA_SAFE_CALL( cudaMemcpy( hI_Jerk_Split, dI_Jerk, Out_Mem_Size, cudaMemcpyDeviceToHost ) );

      for (int i=0; i<Ni; i++)
      {
         for (int dim=0; dim<3; dim++)
         {
            hI_Acc [i][dim] = 0.0;
            hI_Jerk[i][dim] = 0.0;
         }

         for (int Split=0; Split<NSplit; Split++)
         {
            int j = i+Split*Ni;

            for (int dim=0; dim<3; dim++)
            {
               hI_Acc [i][dim] += hI_Acc_Split [j][dim];
               hI_Jerk[i][dim] += hI_Jerk_Split[j][dim];
            }
         }
      } // for (int i=0; i<Ni; i++)
   } // if ( NSplit == 1 ) ... else ...


// multiply results by the gravitational constant if it's not unity
//=========================================================================================
   if ( G != (real)1.0 )
   {
      for (int i=0; i<Ni; i++)
      for (int d=0; d<3; d++)
      {
         hI_Acc [i][d] *= G;
         hI_Jerk[i][d] *= G;
      }
   }

}



//-----------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate
// Description :  Allocate the GPU and PC memory
//-----------------------------------------------------------------------------------------
void CUAPI_MemAllocate()
{

   if ( MyRank == 0 )    fprintf( stdout, "Allocate GPU memory ... " );


// assume that the maximum number of pseudo I-particle is <= 2*SPLIT_NMAX
   const unsigned int Max_Pseudo_N = SPLIT_NMAX*2;
#  ifdef HYBRID_SCHEME
   const unsigned int I_Mem_Size   = TOTAL_N * sizeof(real);
   const unsigned int J_Mem_Size   = N * sizeof(real);
   const unsigned int Out_Mem_Size = ( TOTAL_N >= Max_Pseudo_N ) ?
                                     3*TOTAL_N*sizeof(real) : 3*Max_Pseudo_N*sizeof(real);
#  else
   const unsigned int I_Mem_Size   = N * sizeof(real);
   const unsigned int J_Mem_Size   = N * sizeof(real);
   const unsigned int Out_Mem_Size = ( N >= Max_Pseudo_N ) ?
                                     3*N*sizeof(real) : 3*Max_Pseudo_N*sizeof(real);
#  endif


// allocate the device memory
   CUDA_SAFE_CALL( cudaMalloc( (void**) &dJ_Mass,   J_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMalloc( (void**) &dJ_Pos,  3*J_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMalloc( (void**) &dJ_Vel,  3*J_Mem_Size ) );

   CUDA_SAFE_CALL( cudaMalloc( (void**) &dI_Pos,  3*I_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMalloc( (void**) &dI_Vel,  3*I_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMalloc( (void**) &dI_Acc,  Out_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMalloc( (void**) &dI_Jerk, Out_Mem_Size ) );

   CUDA_SAFE_CALL( cudaMalloc( (void**) &dI_Pot,    J_Mem_Size ) );

   CUDA_SAFE_CALL( cudaMalloc( (void**) &d_Nj_afterSplit_List, 32*sizeof(int) ) );


// allocate the host memory by CUDA
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &Mass,            J_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &Pos_Pred,      3*J_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &Vel_Pred,      3*J_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &Pos_LineUp,    3*I_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &Vel_LineUp,    3*I_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &Acc_local,     3*I_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &Jerk_local,    3*I_Mem_Size ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &hI_Acc_Split,  3*Max_Pseudo_N*sizeof(real) ) );
   CUDA_SAFE_CALL( cudaMallocHost( (void**) &hI_Jerk_Split, 3*Max_Pseudo_N*sizeof(real) ) );

   CUDA_SAFE_CALL( cudaMallocHost( (void**) &h_Nj_afterSplit_List, 32*sizeof(int) ) );


   if ( MyRank == 0 )    fprintf( stdout, "done\n" );

}



//-----------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree
// Description :  Free the memory previous allocated by the function "CUAPI_MemAllocate"
//-----------------------------------------------------------------------------------------
void CUAPI_MemFree()
{

   if ( MyRank == 0 )    fprintf( stdout, "Free GPU memory ... " );


// free the device memory
   if ( dJ_Mass               != NULL )   CUDA_SAFE_CALL(  cudaFree( dJ_Mass              )  );
   if ( dJ_Pos                != NULL )   CUDA_SAFE_CALL(  cudaFree( dJ_Pos               )  );
   if ( dJ_Vel                != NULL )   CUDA_SAFE_CALL(  cudaFree( dJ_Vel               )  );
   if ( dI_Pos                != NULL )   CUDA_SAFE_CALL(  cudaFree( dI_Pos               )  );
   if ( dI_Vel                != NULL )   CUDA_SAFE_CALL(  cudaFree( dI_Vel               )  );
   if ( dI_Acc                != NULL )   CUDA_SAFE_CALL(  cudaFree( dI_Acc               )  );
   if ( dI_Jerk               != NULL )   CUDA_SAFE_CALL(  cudaFree( dI_Jerk              )  );
   if ( dI_Pot                != NULL )   CUDA_SAFE_CALL(  cudaFree( dI_Pot               )  );
   if ( d_Nj_afterSplit_List  != NULL )   CUDA_SAFE_CALL(  cudaFree( d_Nj_afterSplit_List )  );

   dJ_Mass              = NULL;
   dJ_Pos               = NULL;
   dJ_Vel               = NULL;
   dI_Pos               = NULL;
   dI_Vel               = NULL;
   dI_Acc               = NULL;
   dI_Jerk              = NULL;
   dI_Pot               = NULL;
   d_Nj_afterSplit_List = NULL;


// free the host memory allocated by CUDA
   if ( Mass                  != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( Mass                 )  );
   if ( Pos_Pred              != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( Pos_Pred             )  );
   if ( Vel_Pred              != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( Vel_Pred             )  );
   if ( Pos_LineUp            != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( Pos_LineUp           )  );
   if ( Vel_LineUp            != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( Vel_LineUp           )  );
   if ( Acc_local             != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( Acc_local            )  );
   if ( Jerk_local            != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( Jerk_local           )  );
   if ( hI_Acc_Split          != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( hI_Acc_Split         )  );
   if ( hI_Jerk_Split         != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( hI_Jerk_Split        )  );
   if ( h_Nj_afterSplit_List  != NULL )   CUDA_SAFE_CALL(  cudaFreeHost( h_Nj_afterSplit_List )  );

   Mass                 = NULL;
   Pos_Pred             = NULL;
   Vel_Pred             = NULL;
   Pos_LineUp           = NULL;
   Vel_LineUp           = NULL;
   Acc_local            = NULL;
   Jerk_local           = NULL;
   hI_Acc_Split         = NULL;
   hI_Jerk_Split        = NULL;
   h_Nj_afterSplit_List = NULL;

   if ( MyRank == 0 )    fprintf( stdout, "done\n" );

}

