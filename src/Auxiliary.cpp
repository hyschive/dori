#include "Dori.h"
#include <sys/stat.h>




//----------------------------------------------------------------------
// Function    :  TakeNote
// Description :  Record the simulation information
//----------------------------------------------------------------------
void TakeNote( const double INIT_T, const double END_T, const long int INIT_STEP, const long int END_STEP,
               const double ENERGY_DT, const double MOMENTUM_DT, const double OUTPUT_DT, const double DT_DIAGNOSIS_DT,
               const int RESTART, const real INIT_E, const real INIT_L[3], const int INIT_DUMP_ID, const bool BINARY_OUTPUT,
               const bool CONST_INIT_DT, const int GPUID_SELECT )
{

   if ( MyRank == 0 )    fprintf( stdout, "Take note ...\n" );


   const char FileName[] = "Record__Note";

   FILE *Note = fopen( FileName, "a" );
   fprintf( Note, "\nSimulation Note\n" );
   fprintf( Note, "***********************************************************************************\n" );
   fclose( Note );

   if ( system("cat ./Input__NoteScript >> Record__Note") == -1 )
   {
      fprintf( stderr, "ERROR : the system call failed !!\n" );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      exit(-1);
   }

   Note = fopen( FileName, "a" );
   fprintf( Note, "***********************************************************************************\n" );
   fprintf( Note, "\n\n");


   fprintf( Note, "Simulation Parameters\n" );
   fprintf( Note, "***********************************************************************************\n" );
   fprintf( Note, "TOTAL_N             =   %d\n",        TOTAL_N         );
   fprintf( Note, "N                   =   %d\n",        N               );
   fprintf( Note, "NGPU                =   %d\n",        NGPU            );
   fprintf( Note, "OMP_NTHREAD         =   %d\n",        OMP_NTHREAD     );
   fprintf( Note, "EPS                 =   %13.7e\n",    SQRT(EPS_SQR)   );
   fprintf( Note, "NEWTON_G            =   %13.7e\n",    NEWTON_G        );
   fprintf( Note, "ETA                 =   %13.7e\n",    ETA             );
   fprintf( Note, "INIT_ETA            =   %13.7e\n",    INIT_ETA        );
   fprintf( Note, "END_T               =   %13.7e\n",    END_T           );
   fprintf( Note, "END_STEP            =   %ld\n",       END_STEP        );
   fprintf( Note, "MAX_DT              =   %20.14e\n",   MAX_DT          );
   fprintf( Note, "MIN_DT              =   %20.14e\n",   MIN_DT          );
   fprintf( Note, "CONST_INIT_DT       =   %d\n",        CONST_INIT_DT   );
   fprintf( Note, "BINARY_OUTPUT       =   %d\n",        BINARY_OUTPUT   );
   fprintf( Note, "BINARY_ORDER        =   %d\n",        BINARY_ORDER    );
   fprintf( Note, "OUTPUT_DT           =   %13.7e\n",    OUTPUT_DT       );
   fprintf( Note, "ENERGY_DT           =   %13.7e\n",    ENERGY_DT       );
   fprintf( Note, "MOMENTUM_DT         =   %13.7e\n",    MOMENTUM_DT     );
   fprintf( Note, "DT_DIAGNOSIS_DT     =   %13.7e\n",    DT_DIAGNOSIS_DT );
   fprintf( Note, "LOG_STEP            =   %ld\n",       LOG_STEP        );
   fprintf( Note, "RESTART             =   %d\n",        RESTART         );
   fprintf( Note, "INIT_METHOD         =   %d\n",        INIT_METHOD     );
   fprintf( Note, "INIT_T              =   %20.14e\n",   INIT_T          );
   fprintf( Note, "INIT_STEP           =   %ld\n",       INIT_STEP       );
   fprintf( Note, "INIT_DUMP_ID        =   %d\n",        INIT_DUMP_ID    );
   fprintf( Note, "INIT_E              =   %14.7e\n",    INIT_E          );
   fprintf( Note, "INIT_L[0]           =   %14.7e\n",    INIT_L[0]       );
   fprintf( Note, "INIT_L[1]           =   %14.7e\n",    INIT_L[1]       );
   fprintf( Note, "INIT_L[2]           =   %14.7e\n",    INIT_L[2]       );
   fprintf( Note, "SPLIT_NMAX          =   %d\n",        SPLIT_NMAX      );
   fprintf( Note, "GPUID_SELECT        =   %d\n",        GPUID_SELECT    );
   fprintf( Note, "GRAVITY_TYPE        =   %d\n",        GRAVITY_TYPE    );
   if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH ) {
   fprintf( Note, "EXT_METHOD          =   %d\n",        EXT_METHOD      );
   fprintf( Note, "EXT_JERK            =   %d\n",        EXT_JERK        );
   if ( EXT_METHOD == EXT_FILE ) {
   fprintf( Note, "EXT_SIZE[0]         =   %d\n",        EXT_SIZE[0]     );
   fprintf( Note, "EXT_SIZE[1]         =   %d\n",        EXT_SIZE[1]     );
   fprintf( Note, "EXT_SIZE[2]         =   %d\n",        EXT_SIZE[2]     );
   fprintf( Note, "EXT_CEN[0]          =   %14.7e\n",    EXT_CEN[0]      );
   fprintf( Note, "EXT_CEN[1]          =   %14.7e\n",    EXT_CEN[1]      );
   fprintf( Note, "EXT_CEN[2]          =   %14.7e\n",    EXT_CEN[2]      );
   fprintf( Note, "EXT_DH              =   %13.7e\n",    EXT_DH          );
   fprintf( Note, "EXT_PAR_INT         =   %d\n",        EXT_PAR_INT     );
   fprintf( Note, "EXT_ACC_DER         =   %d\n",        EXT_ACC_DER     ); }}
   fprintf( Note, "\n");
   fprintf( Note, "BLOCK_SIZE          =   %d\n",        BLOCK_SIZE      );
   fprintf( Note, "GRID_SIZE           =   %d\n",        GRID_SIZE       );
   fprintf( Note, "***********************************************************************************\n" );
   fprintf( Note, "\n\n");


   fprintf( Note, "Simulation Options\n" );
   fprintf( Note, "***********************************************************************************\n" );

#  ifdef GPU
   fprintf( Note, "GPU                 : ON\n" );
#  else
   fprintf( Note, "GPU                 : OFF\n" );
#  endif

#  ifdef SOFTEN
   fprintf( Note, "SOFTEN              : ON\n" );
#  else
   fprintf( Note, "SOFTEN              : OFF\n" );
#  endif

#  ifdef IMPLICIT_HERMITE
   fprintf( Note, "IMPLICIT_HERMITE    : ON\n" );
#  else
   fprintf( Note, "IMPLICIT_HERMITE    : OFF\n" );
#  endif

#  ifdef N_IS_MULTIPLE_OF_BS
   fprintf( Note, "N_IS_MULTIPLE_OF_BS : ON\n" );
#  else
   fprintf( Note, "N_IS_MULTIPLE_OF_BS : OFF\n" );
#  endif

#  ifdef HYBRID_SCHEME
   fprintf( Note, "HYBRID_SCHEME       : ON\n" );
#  else
   fprintf( Note, "HYBRID_SCHEME       : OFF\n" );
#  endif

#  ifdef FLOAT8_ACC
   fprintf( Note, "FLOAT8_ACC          : ON\n" );
#  else
   fprintf( Note, "FLOAT8_ACC          : OFF\n" );
#  endif

#  ifdef FLOAT8
   fprintf( Note, "FLOAT8              : ON\n" );
#  else
   fprintf( Note, "FLOAT8              : OFF\n" );
#  endif

#  ifdef SHARED_TIMESTEP
   fprintf( Note, "SHARED_TIMESTEP     : ON\n" );
#  else
   fprintf( Note, "SHARED_TIMESTEP     : OFF\n" );
#  endif

#  ifdef GPU
#  if   ( GPU_ARCH == FERMI )
   fprintf( Note, "GPU_ARCH            : FERMI\n" );
#  elif ( GPU_ARCH == KEPLER )
   fprintf( Note, "GPU_ARCH            : KEPLER\n" );
#  elif ( GPU_ARCH == MAXWELL )
   fprintf( Note, "GPU_ARCH            : MAXWELL\n" );
#  elif ( GPU_ARCH == PASCAL )
   fprintf( Note, "GPU_ARCH            : PASCAL\n" );
#  elif ( GPU_ARCH == VOLTA )
   fprintf( Note, "GPU_ARCH            : VOLTA\n" );
#  elif ( GPU_ARCH == TURING )
   fprintf( Note, "GPU_ARCH            : TURING\n" );
#  else
   fprintf( Note, "GPU_ARCH            : UNKNOWN\n" );
#  endif
#  endif

#  ifdef OPENMP
   fprintf( Note, "OPENMP              : ON\n" );
#  else
   fprintf( Note, "OPENMP              : OFF\n" );
#  endif

   fprintf( Note, "***********************************************************************************\n" );
   fprintf( Note, "\n\n" );


// record the parameters of OpenMP
#  ifdef OPENMP
   int omp_nthread, omp_chunk_size, omp_nested;
   omp_sched_t omp_schedule;

   omp_nested = omp_get_nested();
   omp_get_schedule( &omp_schedule, &omp_chunk_size );
#  pragma omp parallel
#  pragma omp master
   { omp_nthread = omp_get_num_threads(); }

   fprintf( Note, "Parameters of OpenMP\n" );
   fprintf( Note, "***********************************************************************************\n" );
   fprintf( Note, "OMP__NUM_THREADS          %d\n",      omp_nthread             );
   fprintf( Note, "OMP__SCHEDULE             %s\n",      ( omp_schedule == omp_sched_static  ) ? "STATIC"  :
                                                         ( omp_schedule == omp_sched_dynamic ) ? "DYNAMIC" :
                                                         ( omp_schedule == omp_sched_guided  ) ? "GUIDED"  :
                                                         ( omp_schedule == omp_sched_auto    ) ? "AUTO"    : "UNKNOWN" );
   fprintf( Note, "OMP__SCHEDULE_CHUNK_SIZE  %d\n",      omp_chunk_size          );
   fprintf( Note, "OMP__NESTED               %s\n",      ( omp_nested ) ? "ON" : "OFF" );
   fprintf( Note, "***********************************************************************************\n" );
   fprintf( Note, "\n\n");
#  endif // #ifdef OPENMP


   fclose( Note );


   if ( MyRank == 0 )    fprintf( stdout, "Take note ... done\n" );

}



//----------------------------------------------------------------------
// Function    :  dt_diagnosis
// Description :  Diagnose the time-step distribution
//----------------------------------------------------------------------
void dt_diagnosis()
{

   const char FileName[]        = "Record__TimeStep";
   static long int PreviousStep = -1;

   if ( Step == PreviousStep )   return;

   if ( MyRank == 0 )   fprintf( stdout, "Diagnose time-step ... " );


#  ifdef SHARED_TIMESTEP
   if ( MyRank == 0 )
   {
      FILE *File = fopen( FileName, "a" );

      if ( PreviousStep == -1 )  // first time of the invocation
         fprintf( File, "%9s    %20s    %20s\n", "Step", "Global_Time", "Time-Step" );

      fprintf( File, "%9ld    %20.14e    %20.14e\n", Step, Global_Time, dt[0] );

      fclose( File );
   }
#  else
   const int MinPow = (int)log2( MIN_DT );
   const int MaxPow = (int)log2( MAX_DT );
   const int NPow   = MaxPow - MinPow + 1;

   int count_local[NPow], count[NPow], Bin;

   for (int i=0; i<NPow; i++)
   {
      count_local[i] = 0;
      count      [i] = 0;
   }

   for (int i=0; i<N; i++)
   {
      Bin = MaxPow - (int)floor( log2(dt[i]) );

      if ( Bin >= 0  &&  Bin < NPow )  count_local[Bin]++;
      else
      {
         fprintf( stderr, "ERROR : out-of-range time-step (dt = %20.14e, Bin = %d, t = %20.14e, Rank = %d\n",
                  dt[i], Bin, Global_Time, MyRank );
         fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
         exit( -1 );
      }
   }

   MPI_Reduce( count_local, count, NPow, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MyRank == 0 )
   {
      FILE *File = fopen( FileName, "a" );

      if ( PreviousStep == -1 )  // first time of the invocation
      {
         char PowStr[3], Sign;
         int Pow;

         fprintf( File, "#%12s ", "Global_Time" );
         for (int Bin=0; Bin<NPow; Bin++)
         {
            Pow  = MaxPow - Bin;
            Sign = ( Pow >= 0 ) ? '+' : '-';

            sprintf( PowStr, "%d%d", abs(Pow)/10, abs(Pow)%10 );
            fprintf( File, "2^%c%2s ", Sign, PowStr );
         }
         fprintf( File, "\n" );
      }

      fprintf( File, "%13.7e ", Global_Time );
      for (int Bin=0; Bin<NPow; Bin++)    fprintf( File, "%5d ", count[Bin] );
      fprintf( File, "\n" );

      fclose( File );
   }


// check if we count all particles
   int Total = 0;
   for (int i=0; i<NPow; i++)    Total += count_local[i];
   if ( Total != N )
   {
      fprintf( stderr, "ERROR : missing particles (Total = %d, N = %d), t = %14.7e, Rank = %d\n",
               Total, N, Global_Time, MyRank );
      fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
      exit( -1 );
   }
#  endif // #ifdef ... else ... SHARED_TIMESTEP


   PreviousStep = Step;

   if ( MyRank == 0 )    fprintf( stdout, "done\n" );

}



//----------------------------------------------------------------------
// Function    :  Get_TotalEnergy
// Description :  Calculate the sum of the kinematic and potential energy
//
// Note        :  Including external potential
//----------------------------------------------------------------------
void Get_TotalEnergy( bool UseInputEgy, real INIT_E )
{

   const char FileName[] = "Record__Energy";

   double PE       = 0.0;
   double KE       = 0.0;
   double TotalE   = 0.0;
   double PE_local = 0.0;
   double KE_local = 0.0;

   static double   TotalE_Init  = 0.0;
   static long int PreviousStep = -1;


   if ( Step != PreviousStep )
   {
      if ( MyRank == 0 )    fprintf( stdout, "Record total energy ... " );


      BeginRing_Pot();

      for (int i=0; i<N; i++)
      {
         PE_local += Mass[i]*Pot[i];
         KE_local += Mass[i]*( Vel[i][0]*Vel[i][0] + Vel[i][1]*Vel[i][1] + Vel[i][2]*Vel[i][2] );
      }

      MPI_Reduce( &PE_local, &PE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &KE_local, &KE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

      KE *= 0.5;
      PE *= 0.5;  // because the potential energy has been double counted

//    add external potential
      if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH )
      {
         double PE_Ext       = 0.0;
         double PE_Ext_local = 0.0;
         real dr[3], r;

         for (int i=0; i<N; i++)
         {
            for (int d=0; d<3; d++)    dr[d] = Pos[i][d] - SOL_CEN[d];

            r             = SQRT( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );
            PE_Ext_local += Mass[i]*Ext_TotalPot( r, Global_Time, SOL_EXT_SC );
         }

         MPI_Reduce( &PE_Ext_local, &PE_Ext, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

         PE += PE_Ext;

      } // if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH )

      TotalE = PE + KE;

      if ( PreviousStep == -1 )  // first time of invocation
      {
         if ( UseInputEgy )   TotalE_Init = INIT_E;
         else                 TotalE_Init = TotalE;
      }

      if ( MyRank == 0 )
      {
         FILE *File = fopen( FileName, "a" );

         if ( PreviousStep == -1 )  // first time of invocation
            fprintf( File, "#%12s  %10s  %14s  %14s  %14s  %14s  %14s\n",
                     "Global_Time", "Step", "K.E.", "P.E", "Total E", "Abs Err", "Rel Err" );

         fprintf( File, "%13.7e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
                  Global_Time, Step, KE, PE, TotalE, TotalE-TotalE_Init, (TotalE-TotalE_Init)/TotalE_Init );

         fclose( File );
      }

      PreviousStep = Step;


      if ( MyRank == 0 )    fprintf( stdout, "done\n" );

   } // if ( Step != PreviousStep )

}



//----------------------------------------------------------------------
// Function    :  Get_TotalMomentum
// Description :  Calculate the total linear and angular momenta
//
// Note        :  1. Invoked by main()
//                2. Only compute angular momentum for now
//
// Parameter   :  UseInputMom : Use the input INIT_L to estimate errors
//                INIT_L      : Use by UseInputMom
//----------------------------------------------------------------------
void Get_TotalMomentum( bool UseInputMom, real INIT_L[3] )
{

   const char FileName[] = "Record__Momentum";

   static double   MomL_Init[3] = { 0.0, 0.0, 0.0 };
   static long int PreviousStep = -1;


   if ( Step != PreviousStep )
   {
      if ( MyRank == 0 )    fprintf( stdout, "Record total momentum ... " );


      double MomL[3]       = { 0.0, 0.0, 0.0 };
      double MomL_local[3] = { 0.0, 0.0, 0.0 };

      for (int i=0; i<N; i++)
      {
         real dr[3];
         for (int d=0; d<3; d++)    dr[d] = Pos[i][d] - SOL_CEN[d];

         MomL_local[0] += dr[1]*Vel[i][2] - dr[2]*Vel[i][1];
         MomL_local[1] += dr[2]*Vel[i][0] - dr[0]*Vel[i][2];
         MomL_local[2] += dr[0]*Vel[i][1] - dr[1]*Vel[i][0];
      }

      MPI_Reduce( MomL_local, MomL, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );


      if ( PreviousStep == -1 )  // first time of invocation
      {
         for (int d=0; d<3; d++)
         {
            if ( UseInputMom )   MomL_Init[d] = INIT_L[d];
            else                 MomL_Init[d] = MomL  [d];
         }
      }

      if ( MyRank == 0 )
      {
         FILE *File = fopen( FileName, "a" );

         if ( PreviousStep == -1 )  // first time of invocation
            fprintf( File, "#%12s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n",
                     "Global_Time", "Step", "MomL_x", "Abs Err", "Rel Err",
                                            "MomL_y", "Abs Err", "Rel Err",
                                            "MomL_z", "Abs Err", "Rel Err" );

         fprintf( File, "%13.7e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
                  Global_Time, Step, MomL[0], MomL[0]-MomL_Init[0], (MomL[0]-MomL_Init[0])/MomL_Init[0],
                                     MomL[1], MomL[1]-MomL_Init[1], (MomL[1]-MomL_Init[1])/MomL_Init[1],
                                     MomL[2], MomL[2]-MomL_Init[2], (MomL[2]-MomL_Init[2])/MomL_Init[2] );

         fclose( File );
      }

      PreviousStep = Step;


      if ( MyRank == 0 )    fprintf( stdout, "done\n" );
   } // if ( Step != PreviousStep )

} // FUNCTION : Get_TotalMomentum



//----------------------------------------------------------------------
// Function    :  OutputData
// Description :  Dump the mass, position, and velocity to file
//----------------------------------------------------------------------
void OutputData( const int Init_DumpID, const bool Binary_Output )
{

   static long int PreviousStep = -1;
   static int DumpID = Init_DumpID;
   FILE *File;

   if ( Step != PreviousStep )
   {
      char FileName[2][10]={"Data_"};

      FileName[1][0]  = 48 + DumpID/100000;
      FileName[1][1]  = 48 + DumpID%100000/10000;
      FileName[1][2]  = 48 + DumpID%10000/1000;
      FileName[1][3]  = 48 + DumpID%1000/100;
      FileName[1][4]  = 48 + DumpID%100/10;
      FileName[1][5]  = 48 + DumpID%10;

      strcat( FileName[0], FileName[1] );

      if ( MyRank == 0 )
      {
         FILE *File_Record = fopen( "Record__Dump", "a" );

         if ( DumpID == Init_DumpID )
            fprintf( File_Record, "#%7s    %20s    %9s\n", "DumpID", "Global_Time", "Step" );

         fprintf( File_Record, "%8d    %20.14e    %9ld\n", DumpID, Global_Time, Step );

         fclose( File_Record );
      }


      if ( MyRank == 0 )    fprintf( stdout, "Dump data (DumpID = %d) ...\n", DumpID );


      for (int YourTurn=0; YourTurn<NGPU; YourTurn++)
      {
         if ( MyRank == YourTurn )
         {
            if ( Binary_Output )
            {
               if ( MyRank == 0 )
               {
                  FILE *File_Check = fopen( FileName[0], "r" );
                  if ( File_Check != NULL )
                  {
                     fprintf( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n",
                              FileName[0] );
                     fclose( File_Check );
                  }

                  File = fopen( FileName[0], "wb" );

//                set file to the target size for BORDER_VP
                  if ( BINARY_ORDER == BORDER_VP )    ftruncate( fileno(File), TOTAL_N*7*sizeof(real) );

                  fclose( File );
               }

               if      ( BINARY_ORDER == BORDER_PV )
               {
                  File = fopen( FileName[0], "ab" );

                  for (int i=0; i<N; i++)
                  {
                     fwrite( &Mass[i],   sizeof(real), 1, File );

                     fwrite( &Pos[i][0], sizeof(real), 1, File );
                     fwrite( &Pos[i][1], sizeof(real), 1, File );
                     fwrite( &Pos[i][2], sizeof(real), 1, File );

                     fwrite( &Vel[i][0], sizeof(real), 1, File );
                     fwrite( &Vel[i][1], sizeof(real), 1, File );
                     fwrite( &Vel[i][2], sizeof(real), 1, File );
                  }
               } // if ( BINARY_ORDER == BORDER_PV )

               else if ( BINARY_ORDER == BORDER_VP )
               {
//                use "rb+" to avoid deleting the content of the file
                  File = fopen( FileName[0], "rb+" );

                  const long offset_v =  TOTAL_N*sizeof(real);
                  const long offset_p = MyRank*N*sizeof(real);

                  real *OneAtt = new real [N];

//                mass
                  for (int i=0; i<N; i++)    OneAtt[i] = Mass[i];
                  fseek( File, 0*offset_v + offset_p, SEEK_SET );
                  fwrite( OneAtt, sizeof(real), N, File );

//                x
                  for (int i=0; i<N; i++)    OneAtt[i] = Pos [i][0];
                  fseek( File, 1*offset_v + offset_p, SEEK_SET );
                  fwrite( OneAtt, sizeof(real), N, File );

//                y
                  for (int i=0; i<N; i++)    OneAtt[i] = Pos [i][1];
                  fseek( File, 2*offset_v + offset_p, SEEK_SET );
                  fwrite( OneAtt, sizeof(real), N, File );

//                z
                  for (int i=0; i<N; i++)    OneAtt[i] = Pos [i][2];
                  fseek( File, 3*offset_v + offset_p, SEEK_SET );
                  fwrite( OneAtt, sizeof(real), N, File );

//                vx
                  for (int i=0; i<N; i++)    OneAtt[i] = Vel [i][0];
                  fseek( File, 4*offset_v + offset_p, SEEK_SET );
                  fwrite( OneAtt, sizeof(real), N, File );

//                vy
                  for (int i=0; i<N; i++)    OneAtt[i] = Vel [i][1];
                  fseek( File, 5*offset_v + offset_p, SEEK_SET );
                  fwrite( OneAtt, sizeof(real), N, File );

//                vz
                  for (int i=0; i<N; i++)    OneAtt[i] = Vel [i][2];
                  fseek( File, 6*offset_v + offset_p, SEEK_SET );
                  fwrite( OneAtt, sizeof(real), N, File );

                  delete [] OneAtt;
               } // else if ( BINARY_ORDER == BORDER_VP )

               else
                  Aux_Error( ERROR_INFO, "unsupported BINARY_ORDER (%d) !!\n", BINARY_ORDER );

               fclose( File );
            } // if ( Binary_Output )

            else
            {
               if ( MyRank == 0 )
               {
                  FILE *File_Check = fopen( FileName[0], "r" );
                  if ( File_Check != NULL )
                  {
                     fprintf( stderr, "WARNING : the file \"%s\" already exists and will be overwritten !!\n",
                              FileName[0] );
                     fclose( File_Check );
                  }

                  File = fopen( FileName[0],"w" );

                  fprintf( File, "#%13s  %14s  %14s  %14s  %14s  %14s  %14s\n",
                           "Mass", "x", "y", "z", "Vx", "Vy", "Vz" );
               }
               else
                  File = fopen( FileName[0],"a" );

               for (int i=0; i<N; i++)
               {
                  fprintf( File, "%14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
                           Mass[i], Pos[i][0], Pos[i][1], Pos[i][2], Vel[i][0], Vel[i][1], Vel[i][2] );
               }

               fclose(File);
            } // if ( Binary_Output ) ... else ...

         } // if ( MyRank == YourTurn )

         MPI_Barrier( MPI_COMM_WORLD );
      } // for (int YourTurn=0; YourTurn<NGPU; YourTurn++)


//    verify the file size
      if ( Binary_Output  &&  MyRank == 0 )
      {
         long int ExpectSize, InputSize;

         FILE *FileJustWrite = fopen( FileName[0], "rb" );

         if ( FileJustWrite == NULL )
         {
            fprintf( stderr, "ERROR : the file \"%s\" does not exist !!\n", FileName[0] );
            exit(-1);
         }

         fseek( FileJustWrite, 0, SEEK_END );

         InputSize  = ftell( FileJustWrite );
         ExpectSize = TOTAL_N*7*sizeof(real);

         fclose( FileJustWrite );

         if ( InputSize > ExpectSize )
         {
            fprintf( stderr, "WARNING : the size of the output file <%s> is too large !!\n", FileName[0] );
            fprintf( stderr, "          Output = %ld bytes <-> Expect = %ld bytes ...\n", InputSize, ExpectSize );
            fprintf( stderr, "          file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
            MPI_Exit();
         }

         if ( InputSize < ExpectSize )
         {
            fprintf( stderr, "ERROR : the size of the output file <%s> is too small !!\n", FileName[0] );
            fprintf( stderr, "        Output = %ld bytes <-> Expect = %ld bytes ...\n", InputSize, ExpectSize );
            fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
            MPI_Exit();
         }
      } // if ( MyRank == 0 )


      if ( SOL_REC_DIS )
      {
         if ( N != 1 )     Aux_Error( ERROR_INFO, "SOL_REC_DIS must work with N=1 !!\n" );
         if ( NGPU != 1 )  Aux_Error( ERROR_INFO, "SOL_REC_DIS must work with NGPU=1 !!\n" );

         const int  ParID = 0;
         real dr[3], r;

         for (int d=0; d<3; d++)    dr[d] = Pos[ParID][d] - SOL_CEN[d];
         r = SQRT( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );

         const char FileName_Dis[] = "Record__Distance";
         FILE *File_Dis = fopen( FileName_Dis, "a" );

         if ( PreviousStep == -1 )  // first time of invocation
            fprintf( File_Dis, "#%14s  %10s  %15s\n", "Time [Gyr]", "Step", "Distance [kpc]" );

         fprintf( File_Dis, "%15.7e  %10ld  %15.7e\n", Global_Time, Step, r );

         fclose( File_Dis );
      } // if ( SOL_REC_DIS )


      if ( SOL_REC_HLR )
      {
//       compute the center of mass
         const double TotalM = Ext_GetCM();


//       collect all particles to the root rank
         real (*MassAll)   = ( MyRank == 0 ) ? new real [TOTAL_N]    : NULL;
         real (*PosAll)[3] = ( MyRank == 0 ) ? new real [TOTAL_N][3] : NULL;

         MPI_Gather( Mass,    N,   GAMER_MPI_REAL,
                     MassAll, N,   GAMER_MPI_REAL, 0, MPI_COMM_WORLD );

         MPI_Gather( Pos,     N*3, GAMER_MPI_REAL,
                     PosAll,  N*3, GAMER_MPI_REAL, 0, MPI_COMM_WORLD );


//       compute the half-light radius
//       --> strictly speaking, we are computing the 2D half-mass radius (i.e., the radius where the 2D projected
//           mass profile is half of the total mass)
         if ( MyRank == 0 )
         {
            double (*R2D) = new double [TOTAL_N];
            int    (*Idx) = new int    [TOTAL_N];
            double dx, dy, HLR, EnclosedMass=0.0;

            for (int i=0; i<TOTAL_N; i++)
            {
//             project along the z direction
               dx     = PosAll[i][0] - SOL_SC_CM[0];
               dy     = PosAll[i][1] - SOL_SC_CM[1];
               R2D[i] = sqrt( SQR(dx) + SQR(dy) );
            }

            Mis_Heapsort( TOTAL_N, R2D, Idx );

            for (int i=0; i<TOTAL_N; i++)
            {
               EnclosedMass += MassAll[ Idx[i] ];

               if ( EnclosedMass > 0.5*TotalM )
               {
                  HLR = R2D[i];
                  break;
               }
            }

            delete [] R2D;
            delete [] Idx;


//          output
            const char FileName_HLR[] = "Record__HalfLightRadius";
            FILE *File_HLR = fopen( FileName_HLR, "a" );

            if ( PreviousStep == -1 )  // first time of invocation
            {
               fprintf( File_HLR, "# Total mass  : %13.7e Msun\n", TotalM );
               fprintf( File_HLR, "# Time   unit : Gyr\n" );
               fprintf( File_HLR, "# Length unit : kpc\n" );
               fprintf( File_HLR, "\n" );
               fprintf( File_HLR, "#%12s  %10s  %13s  %14s  %14s  %14s  %14s\n",
                        "Time", "Step", "R_HL", "CM[x]", "CM[y]", "CM[z]", "R_CM" );
            }

            const double R_SC_Soliton = sqrt(  SQR( SOL_SC_CM[0]-SOL_CEN[0] ) +
                                               SQR( SOL_SC_CM[1]-SOL_CEN[1] ) +
                                               SQR( SOL_SC_CM[2]-SOL_CEN[2] )  );
            fprintf( File_HLR, "%13.7e  %10ld  %13.7e  %14.7e  %14.7e  %14.7e  %14.7e\n",
                     Global_Time, Step, HLR, SOL_SC_CM[0], SOL_SC_CM[1], SOL_SC_CM[2], R_SC_Soliton );

            fclose( File_HLR );
         } // if ( MyRank == 0 )

         delete [] MassAll;
         delete [] PosAll;
      } // if ( SOL_REC_HLR )


      if ( MyRank == 0 )    fprintf( stdout, "Dump data (DumpID = %d) ... done\n", DumpID );

      DumpID++;
      PreviousStep = Step;
   } // if ( Step != PreviousStep )

}



//----------------------------------------------------------------------
// Function    :  MemoryAllocate
// Description :  Allocate memory
//----------------------------------------------------------------------
void MemoryAllocate()
{

   if ( MyRank == 0 )    fprintf( stdout, "Allocate CPU memory ... " );


   t    = new double[N];
   dt   = new double[N];

   Pos  = new real[N][3];
   Vel  = new real[N][3];
   Acc  = new real[N][3];
   Jerk = new real[N][3];
   Pot  = new real[N];

   PlayerList     = new int   [N];
   dt_BeforeSync  = new double[N];

#  ifndef GPU
   const unsigned int Nj = N;
#  ifdef HYBRID_SCHEME
   const unsigned int Ni = TOTAL_N;
#  else
   const unsigned int Ni = N;
#  endif

   Mass        = new real [Nj];
   Pos_Pred    = new real [Nj][3];
   Vel_Pred    = new real [Nj][3];
   Pos_LineUp  = new real [Ni][3];
   Vel_LineUp  = new real [Ni][3];
   Acc_local   = new real [Ni][3];
   Jerk_local  = new real [Ni][3];
#  endif // #ifndef GPU

   if (  ( GRAVITY_TYPE == GRAVITY_EXTERNAL || GRAVITY_TYPE == GRAVITY_BOTH )  &&  EXT_METHOD == EXT_FILE  )
   {
      Ext_AccX = new real [ EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2] ];
      Ext_AccY = new real [ EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2] ];
      Ext_AccZ = new real [ EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2] ];

      if ( EXT_JERK )
      {
         Ext_dAccX = new real [ EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2]*3 ];  // 3 -> x/y/z derivatives
         Ext_dAccY = new real [ EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2]*3 ];
         Ext_dAccZ = new real [ EXT_SIZE[0]*EXT_SIZE[1]*EXT_SIZE[2]*3 ];
      }
   }

   if ( MyRank == 0 )    fprintf( stdout, "done\n" );

} // FUNCTION : MemoryAllocate



//----------------------------------------------------------------------
// Function    :  MemoryFree
// Description :  Deallocate memory
//----------------------------------------------------------------------
void MemoryFree()
{

   if ( MyRank == 0 )    fprintf( stdout, "Free CPU memory ... " );


   if ( t             != NULL )  delete [] t;
   if ( dt            != NULL )  delete [] dt;

   if ( Pos           != NULL )  delete [] Pos;
   if ( Vel           != NULL )  delete [] Vel;
   if ( Acc           != NULL )  delete [] Acc;
   if ( Jerk          != NULL )  delete [] Jerk;
   if ( Pot           != NULL )  delete [] Pot;

   if ( PlayerList    != NULL )  delete [] PlayerList;
   if ( dt_BeforeSync != NULL )  delete [] dt_BeforeSync;

#  ifndef GPU
   if ( Mass          != NULL )  delete [] Mass;
   if ( Pos_Pred      != NULL )  delete [] Pos_Pred;
   if ( Vel_Pred      != NULL )  delete [] Vel_Pred;
   if ( Pos_LineUp    != NULL )  delete [] Pos_LineUp;
   if ( Vel_LineUp    != NULL )  delete [] Vel_LineUp;
   if ( Acc_local     != NULL )  delete [] Acc_local;
   if ( Jerk_local    != NULL )  delete [] Jerk_local;
#  endif

   if ( Ext_AccX      != NULL )  delete [] Ext_AccX;
   if ( Ext_AccY      != NULL )  delete [] Ext_AccY;
   if ( Ext_AccZ      != NULL )  delete [] Ext_AccZ;

   if ( Ext_dAccX     != NULL )  delete [] Ext_dAccX;
   if ( Ext_dAccY     != NULL )  delete [] Ext_dAccY;
   if ( Ext_dAccZ     != NULL )  delete [] Ext_dAccZ;



   t              = NULL;
   dt             = NULL;
   Pos            = NULL;
   Vel            = NULL;
   Acc            = NULL;
   Jerk           = NULL;
   Pot            = NULL;
   PlayerList     = NULL;
   dt_BeforeSync  = NULL;

#  ifndef GPU
   Mass           = NULL;
   Pos_Pred       = NULL;
   Vel_Pred       = NULL;
   Pos_LineUp     = NULL;
   Vel_LineUp     = NULL;
   Acc_local      = NULL;
   Jerk_local     = NULL;
#  endif

   Ext_AccX       = NULL;
   Ext_AccY       = NULL;
   Ext_AccZ       = NULL;

   Ext_dAccX      = NULL;
   Ext_dAccY      = NULL;
   Ext_dAccZ      = NULL;


   if ( MyRank == 0 )    fprintf( stdout, "done\n" );

} // FUNCTION : MemoryFree



//----------------------------------------------------------------------
// Function    :  CheckParameter
// Description :  Verify some parameters loaded from the file "Input__Parameter"
//----------------------------------------------------------------------
void CheckParameter( const double INIT_T, const double END_T, const long int INIT_STEP, const long int END_STEP,
                     const double OUTPUT_DT, const double ENERGY_DT, const double MOMENTUM_DT )
{

   if ( MyRank == 0 )   printf( "Check parameter ...\n" );


#  ifndef SOFTEN
#     error ERROR : currently the program must use soften !!
#  endif

#  if ( defined FLOAT8  &&  !defined FLOAT8_ACC )
#     error ERROR : FLOAT8 must work with FLOAT8_ACC !!
#  endif

#  if ( defined OPENMP  &&  !defined _OPENMP )
#     error : ERROR : something is wrong in OpenMP, the macro "_OPENMP" is NOT defined !!
#  endif

   int NRank;
   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
   if ( NRank != NGPU )
   {
      fprintf( stderr, "ERROR : number of MPI ranks (%d) != number of GPUs (%d) !!\n", NRank, NGPU );
      exit(-1);
   }

   if ( TOTAL_N % NGPU != 0 )
   {
      fprintf( stderr, "ERROR : TOTAL_N (%d) must be a multiple of NGPU (%d) !!\n", TOTAL_N, NGPU );
      exit(-1);
   }

   if ( TOTAL_N != N*NGPU )
   {
      fprintf( stderr, "ERROR : TOTAL_N (%d) != N*NGPU (%d) !!\n", TOTAL_N, N*NGPU );
      exit(-1);
   }

   if ( N <= 0 )
   {
      fprintf( stderr, "ERROR : N (%d) < 0 !!\n", N );
      exit(-1);
   }

   if ( N < 32  &&  N <= SPLIT_NMAX/16 )
   {
      fprintf( stderr, "ERROR : N (%d) < 32 !!\n", N );
      fprintf( stderr, "        --> You should set SPLIT_NMAX smaller\n" );
      exit(-1);
   }

#  ifdef N_IS_MULTIPLE_OF_BS
   if ( N % BLOCK_SIZE != 0 )
   {
      fprintf( stderr, "ERROR : N (%d) is not a multiple of BLOCK_SIZE (%d) for the option ", N, BLOCK_SIZE );
      fprintf( stderr, "\"N_IS_MULTIPLE_OF_BS\" !!\n" );
      exit(-1);
   }
#  else
   if ( N % BLOCK_SIZE == 0 )
   {
      fprintf( stderr, "WARNING : N (%d) is a multiple of BLOCK_SIZE (%d)\n", N, BLOCK_SIZE );
      fprintf( stderr, "          --> You can turn on the option \"N_IS_MULTIPLE_OF_BS\" in the Makefile " );
      fprintf( stderr,           "to improve the performance !!\n" );
   }
#  endif

#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {
      if ( OMP_NTHREAD != omp_get_num_threads() )
         Aux_Message( stderr, "WARNING : OMP_NTHREAD (%d) != omp_get_num_threads (%d) at MPI_Rank %d !!\n",
                      OMP_NTHREAD, omp_get_num_threads(), MyRank );
   }
#  endif

   if ( INIT_T > END_T )
   {
      fprintf( stderr, "ERROR : INIT_T (%14.7e) > END_T (%14.7e) !!\n", INIT_T, END_T );
      exit(-1);
   }

   if ( INIT_STEP > END_STEP )
   {
      fprintf( stderr, "ERROR : INIT_STEP (%ld) > END_STEP (%ld) !!\n", INIT_STEP, END_STEP );
      exit(-1);
   }

   if ( BINARY_ORDER != BORDER_PV  &&  BINARY_ORDER != BORDER_VP )
      Aux_Error( ERROR_INFO, "unsupported BINARY_ORDER (%d) !!\n", BINARY_ORDER );

   if ( GRAVITY_TYPE != GRAVITY_SELF  &&  GRAVITY_TYPE != GRAVITY_EXTERNAL  &&  GRAVITY_TYPE != GRAVITY_BOTH )
      Aux_Error( ERROR_INFO, "unsupported GRAVITY_TYPE (%d) !!\n", GRAVITY_TYPE );

   if ( INIT_METHOD != INIT_FILE  &&  INIT_METHOD != INIT_FUNC )
      Aux_Error( ERROR_INFO, "unsupported INIT_METHOD (%d) !!\n", INIT_METHOD );

   if (  GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH  ) {
   if ( EXT_METHOD != EXT_FILE  &&  EXT_METHOD != EXT_FUNC )
      Aux_Error( ERROR_INFO, "unsupported EXT_METHOD (%d) !!\n", EXT_METHOD );

   if ( EXT_METHOD == EXT_FILE ) {
   if ( EXT_SIZE[0] <= 0  ||  EXT_SIZE[1] <= 0  ||  EXT_SIZE[2] <= 0 )
      Aux_Error( ERROR_INFO, "incorrect EXT_SIZE (%d, %d, %d) !!\n", EXT_SIZE[0], EXT_SIZE[1], EXT_SIZE[2] );

   if ( EXT_DH <= 0.0 )
      Aux_Error( ERROR_INFO, "EXT_DH = %14.7e <= 0.0 !!\n", EXT_DH );

   if ( EXT_PAR_INT != EXT_PAR_INT_CIC  &&  EXT_PAR_INT != EXT_PAR_INT_TSC )
      Aux_Error( ERROR_INFO, "unsupported EXT_PAR_INT (%d) !!\n", EXT_PAR_INT );

   if ( EXT_ACC_DER != EXT_ACC_DER_QUAD  &&  EXT_ACC_DER != EXT_ACC_DER_QUAR )
      Aux_Error( ERROR_INFO, "unsupported EXT_ACC_DER (%d) !!\n", EXT_ACC_DER );
   } // if ( EXT_METHOD == EXT_FILE )
   } // if (  GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH  )

#  ifndef SHARED_TIMESTEP
   if ( fmod(log2(MAX_DT), 1.0) != 0.0 )
   {
      fprintf( stderr, "ERROR : MAX_DT (%14.7e) is not a power of 2 !!\n", MAX_DT );
      exit(-1);
   }

   if ( fmod(INIT_T, MIN_DT) != 0.0 )
   {
      fprintf( stderr, "ERROR : INIT_T (%14.7e) is not a multiple of MIN_DT (%14.7e) !!\n", INIT_T, MIN_DT );
      exit(-1);
   }

   if ( OUTPUT_DT >= 0.0  &&  OUTPUT_DT < MAX_DT )
   {
      fprintf( stderr, "WARNING : setting OUTPUT_DT (%20.14e) < MAX_DT (%20.14e) will introduce extra\n",
               OUTPUT_DT, MAX_DT );
      fprintf( stderr, "          synchronization and hence deteriorate the overall performance !!\n" );
   }

   if ( ENERGY_DT >= 0.0  &&  ENERGY_DT < MAX_DT )
   {
      fprintf( stderr, "WARNING : setting ENERGY_DT (%20.14e) < MAX_DT (%20.14e) will introduce extra\n",
               ENERGY_DT, MAX_DT );
      fprintf( stderr, "          synchronization and hence deteriorate the overall performance !!\n" );
   }

   if ( MOMENTUM_DT >= 0.0  &&  MOMENTUM_DT < MAX_DT )
   {
      fprintf( stderr, "WARNING : setting MOMENTUM_DT (%20.14e) < MAX_DT (%20.14e) will introduce extra\n",
               MOMENTUM_DT, MAX_DT );
      fprintf( stderr, "          synchronization and hence deteriorate the overall performance !!\n" );
   }
#  endif


   if ( MyRank == 0 )   printf( "Check parameter ... done\n" );

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCPUInfo
// Description :  Record the CPU information
//
// Parameter   :  FileName : The name of the output file
//-------------------------------------------------------------------------------------------------------
void GetCPUInfo( const char *FileName )
{

   FILE *Note = fopen( FileName, "a" );
   char *line = NULL;
   size_t len = 0;
   char String[2][100];


// 1. get the CPU info
   const char *CPUInfo_Path = "/proc/cpuinfo";
   FILE *CPUInfo = fopen( CPUInfo_Path, "r" );

   if ( CPUInfo == NULL )
   {
      fprintf( stderr, "WARNING : the CPU information file \"%s\" does not exist !!\n", CPUInfo_Path );
      return;
   }

   while ( getline(&line, &len, CPUInfo) != -1 )
   {
      sscanf( line, "%s%s", String[0], String[1] );

      if (  strcmp( String[0], "model" ) == 0  &&  strcmp( String[1], "name" ) == 0  )
      {
         strncpy( line, "CPU Type  ", 10 );
         fprintf( Note, "%s", line );
         break;
      }
   }

   if ( line != NULL )
   {
      free( line );
      line = NULL;
   }

   fclose( CPUInfo );


// 2. get the memory info
   const char *MemInfo_Path = "/proc/meminfo";
   FILE *MemInfo = fopen( MemInfo_Path, "r" );

   if ( MemInfo == NULL )
   {
      fprintf( stderr, "WARNING : the memory information file \"%s\" does not exist !!\n", MemInfo_Path );
      return;
   }

   while ( getline(&line, &len, MemInfo) != -1 )
   {
      sscanf( line, "%s%s", String[0], String[1] );

      if (  strncmp( String[0], "MemTotal", 8 ) == 0  )
      {
         fprintf( Note, "Total Memory    : %d MB\n", atoi( String[1] )/1024 );
         break;
      }
   }

   if ( line != NULL )
   {
      free( line );
      line = NULL;
   }

   fclose( MemInfo );
   fclose( Note );

}



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Message
// Description :  Output the message and flush it
//
// Note        :  Use the variable argument lists provided in "cstdarg"
//
// Parameter   :  Type     : stdout/stderr/file stream
//                Format   : Output format
//                ...      : Arguments in vfprintf
//                           --> It is equivalent to call "fprintf( Type, Format, ... );   fflush( Type );"
//-------------------------------------------------------------------------------------------------------
void Aux_Message( FILE *Type, const char *Format, ... )
{

// flush all previous messages
   fflush( stdout ); fflush( stdout ); fflush( stdout );
   fflush( stderr ); fflush( stderr ); fflush( stderr );

   va_list Arg;
   va_start( Arg, Format );

   vfprintf( Type, Format, Arg );
   fflush( Type ); fflush( Type ); fflush( Type );

   va_end( Arg );

} // FUNCTION : Aux_Message



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Error
// Description :  Output the error messages and force the program to be terminated
//
// Note        :  Use the variable argument lists provided in "cstdarg"
//
// Parameter   :  File     : Name of the file where error occurs
//                Line     : Line number where error occurs
//                Func     : Name of the function where error occurs
//                Format   : Output format
//                ...      : Arguments in vfprintf
//-------------------------------------------------------------------------------------------------------
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... )
{

// flush all previous messages
   fflush( stdout ); fflush( stdout ); fflush( stdout );
   fflush( stderr ); fflush( stderr ); fflush( stderr );


// output error messages
   va_list Arg;
   va_start( Arg, Format );

   Aux_Message( stderr, "\n" );
   Aux_Message ( stderr, "********************************************************************************\n" );
   Aux_Message ( stderr, "ERROR : " );
   vfprintf    ( stderr, Format, Arg );
   Aux_Message ( stderr, "        Rank <%d>, file <%s>, line <%d>, function <%s>\n",
                 MyRank, File, Line, Func );
   Aux_Message ( stderr, "********************************************************************************\n" );
   Aux_Message( stderr, "\n" );

   va_end( Arg );


// terminate the program
   MPI_Exit();

} // FUNCTION : Aux_Error



//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_Exit
// Description :  Force the program to be terminated when any error occurs
//-------------------------------------------------------------------------------------------------------
void MPI_Exit()
{

// flush all previous messages
   fflush( stdout ); fflush( stdout ); fflush( stdout );
   fflush( stderr ); fflush( stderr ); fflush( stderr );

   Aux_Message( stderr, "\nProgram terminated with error ...... rank %d\n\n", MyRank );

   MPI_Abort( MPI_COMM_WORLD, 0 );

   exit(1);

} // FUNCTION : MPI_Exit



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_CheckFileExist
// Description :  Check whether or not the target file exists
//
// Note        :  Use the "stat" function to query the existence of the target file
//
// Parameter   :  FileName : Name of the target file
//
// Return      :  true/false <-> file exists/not exists
//-------------------------------------------------------------------------------------------------------
bool Aux_CheckFileExist( const char *FileName )
{

   struct stat Buf;
   return ( stat(FileName,&Buf) == 0 );

} // FUNCTION : Aux_CheckFileExist



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_IsFinite
// Description :  Check whether the input floating-point value is finite
//
// Note        :  1. Definition of "finite" --> not NaN, Inf, -Inf
//                2. Alternative to the built-in function isfinite() (or std::isfinite()) which is less portable
//
// Parameter   :  x : Floating-point value to be checked
//
// Return      :  1 : finite
//                0 : not finite
//-------------------------------------------------------------------------------------------------------
int Aux_IsFinite( const float x )
{

   if ( x != x  ||  x < -__FLT_MAX__  ||  x > __FLT_MAX__ )    return 0;
   else                                                        return 1;

} // FUNCTION : Aux_IsFinite



//-------------------------------------------------------------------------------------------------------
// Function overloading: double precision
//-------------------------------------------------------------------------------------------------------
int Aux_IsFinite( const double x )
{

   if ( x != x  ||  x < -__DBL_MAX__  ||  x > __DBL_MAX__ )    return 0;
   else                                                        return 1;

} // FUNCTION : Aux_IsFinite



template <typename T>
static void Heapsort_SiftDown( const int L, const int R, T Array[], int IdxTable[] );

//OPTIMIZATION : (1) quick sort  (2) try the "qsort" library
//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_Heapsort
// Description :  Use the Heapsort algorithm to sort the input array into ascending numerical order
//                --> An index table will also be constructed if "IdxTable != NULL"
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//                2. Overloaded with different types
//                   --> Explicit template instantiation is put in the end of this file
//
// Parameter   :  N        :  Size of Array
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Mis_Heapsort( const int N, T Array[], int IdxTable[] )
{

// initialize the IdxTable
   if ( IdxTable != NULL )
      for (int t=0; t<N; t++)    IdxTable[t] = t;

// heap creation
   for (int L=N/2-1; L>=0; L--)  Heapsort_SiftDown( L, N-1, Array, IdxTable );

// retirement-and-promotion
   T Buf;
   for (int R=N-1; R>0; R--)
   {
      Buf      = Array[R];
      Array[R] = Array[0];
      Array[0] = Buf;

      if ( IdxTable != NULL )
      {
         Buf         = IdxTable[R];
         IdxTable[R] = IdxTable[0];
         IdxTable[0] = Buf;
      }

      Heapsort_SiftDown( 0, R-1, Array, IdxTable );
   }

} // FUNCTION : Mis_Heapsort



//-------------------------------------------------------------------------------------------------------
// Function    :  Heapsort_SiftDown
// Description :  Sift-down process for the Heapsort algorithm
//
// Note        :  1. Ref : Numerical Recipes Chapter 8.3 - 8.4
//                2. Overloaded with different types
//
// Parameter   :  L        :  Left  range of the sift-down
//                R        :  Right range of the sift-down
//                Array    :  Array to be sorted into ascending numerical order
//                IdxTable :  Index table
//-------------------------------------------------------------------------------------------------------
template <typename T>
void Heapsort_SiftDown( const int L, const int R, T Array[], int IdxTable[] )
{

   int  Idx_up    = L;
   int  Idx_down  = 2*Idx_up + 1;
   T    Target    = Array[Idx_up];
   int  TargetIdx = ( IdxTable == NULL ) ? -1 : IdxTable[Idx_up];

   while ( Idx_down <= R )
   {
//    find the better employee
      if ( Idx_down < R  &&  Array[Idx_down+1] > Array[Idx_down] )   Idx_down ++;

//    terminate the sift-down process if the target (supervisor) is better than both its employees
      if ( Target >= Array[Idx_down] )    break;

//    otherwise, promote the better employee
      Array[Idx_up] = Array[Idx_down];
      if ( IdxTable != NULL )    IdxTable[Idx_up] = IdxTable[Idx_down];

//    prepare the next sift-down operation
      Idx_up   = Idx_down;
      Idx_down = 2*Idx_up + 1;
   }

// put target at its best position
   Array[Idx_up] = Target;
   if ( IdxTable != NULL )    IdxTable[Idx_up] = TargetIdx;

} // FUNCTION : Heapsort_SiftDown



// explicit template instantiation
template void Mis_Heapsort <int>    ( const int N, int    Array[], int IdxTable[] );
template void Mis_Heapsort <long>   ( const int N, long   Array[], int IdxTable[] );
template void Mis_Heapsort <ulong>  ( const int N, ulong  Array[], int IdxTable[] );
template void Mis_Heapsort <float>  ( const int N, float  Array[], int IdxTable[] );
template void Mis_Heapsort <double> ( const int N, double Array[], int IdxTable[] );

