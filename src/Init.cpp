#include "Dori.h"




//---------------------------------------------------------------------------
// Function    :  Init_MPI
// Description :  Initialize the MPI environment
//---------------------------------------------------------------------------
void Init_MPI( int argc, char *argv[] )
{

   int NRank;

   MPI_Init( &argc, &argv );
   MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );
   MPI_Comm_size( MPI_COMM_WORLD, &NRank );

// the direction of sending data is counterclockwise
// For example: NGPU = 10, MyRank = (0, 1, ...,8, 9) <--> RecvRank = (1, 2, ...,9, 0)
//                                                   <--> SendRank = (9, 0, ...,7, 8)
   RecvRank = (MyRank+1)%NRank;
   SendRank = (MyRank+NRank-1)%NRank;


   if ( MyRank == 0 )    fprintf( stdout, "Initialize MPI ... done\n" );

}



//----------------------------------------------------------------------
// Function    :  Init_t_dt_step
// Description :  Initialize t, dt and step
//----------------------------------------------------------------------
void Init_t_dt_step( const double INIT_T, const long int INIT_STEP, double &Energy_t, double &Output_t,
                     double &dt_diagnosis_t, const double ENERGY_DT, const double OUTPUT_DT,
                     const double DT_DIAGNOSIS_DT, const bool CONST_INIT_DT )
{

   if ( MyRank == 0 )    fprintf( stdout, "Initialize time-step ...\n" );


   double Min_dt = 1.e10;

// set initial t and step
   for (int i=0; i<N; i++)    t[i] = INIT_T;
   Global_Time = INIT_T;
   Step = INIT_STEP;

// initialize Energy_t, Output_t, and dt_diagnosis_t
   Energy_t       = Global_Time + ENERGY_DT;
   Output_t       = Global_Time + OUTPUT_DT;
   dt_diagnosis_t = Global_Time + DT_DIAGNOSIS_DT;


// set the initial time-step = INIT_ETA * Acc/Jerk
   for (int i=0; i<N; i++)
   {
      double a     = 0.0;
      double a_dot = 0.0;

      for (int dim=0; dim<3; dim++)
      {
          a     += Acc [i][dim]*Acc [i][dim];
          a_dot += Jerk[i][dim]*Jerk[i][dim];
      }

//    check if a_dot == 0 ( which will happen in the case that all initial velocity == 0 )
      if ( a_dot == 0.0 )
      {
         Aux_Message( stderr, "WARNING : initial jerk == 0 for the particle %d at rank %d!!\n", i, MyRank );
         Aux_Message( stderr, "          --> set dt = MIN_DT (%13.7e)\n", MIN_DT );

         dt[i] = MIN_DT;
      }

      else
         dt[i] = INIT_ETA * sqrt(a/a_dot);


#     ifndef SHARED_TIMESTEP
      double Pow = floor( log2(dt[i]) );
      dt[i] = pow( 2.0, Pow );

//    INIT_T must be a multiple of dt
      while ( fmod(INIT_T, dt[i]) != 0.0 )   dt[i] /= 2.0;

      if ( dt[i] < MIN_DT )   dt[i] = MIN_DT;
      if ( dt[i] > MAX_DT )   dt[i] = MAX_DT;

      if ( CONST_INIT_DT && dt[i] < Min_dt )   Min_dt = dt[i];
#     else
      if ( dt[i] < Min_dt )   Min_dt = dt[i];
#     endif
   }


// get the minimum dt from all processes
#  ifndef SHARED_TIMESTEP
   if ( CONST_INIT_DT )
   {
#  endif

      Min_dt = Get_Minimum_FromAllRanks( Min_dt );
      for (int i=0; i<N; i++)    dt[i] = Min_dt;

#  ifndef SHARED_TIMESTEP
   }
#  endif


   if ( MyRank == 0 )    fprintf( stdout, "Initialize time-step ... done\n" );

}



//---------------------------------------------------------------------------------------------------
// Function    :  Init_Particles
// Description :  Initialize the particles's mass, position, velocity, acceleration, and jerk;
//                Set the initial PlayerList
//---------------------------------------------------------------------------------------------------
void Init_Particles( const double INIT_T )
{

   if ( MyRank == 0 )    fprintf( stdout, "Initialize particle data ...\n" );


   if ( INIT_METHOD == INIT_FILE )
   {
      const int offset      = MyRank*N*7*sizeof(real);
      const char FileName[] = "INIT_CONDITION";

//    verify the file size
      if ( MyRank == 0 )
      {
         long int ExpectSize, InputSize;

         FILE *FileCheck = fopen( FileName, "rb" );

         if ( FileCheck == NULL )
         {
            fprintf( stderr, "ERROR : the file \"%s\" does not exist !!\n", FileName );
            exit(-1);
         }

         fseek( FileCheck, 0, SEEK_END );

         InputSize  = ftell( FileCheck );
         ExpectSize = TOTAL_N*7*sizeof(real);

         fclose( FileCheck );

         if ( InputSize > ExpectSize )
         {
            fprintf( stderr, "WARNING : the size of the input file <%s> is too large !!\n", FileName );
            fprintf( stderr, "          Input = %ld bytes <-> Expect = %ld bytes ...\n", InputSize, ExpectSize );
            fprintf( stderr, "          --> The extra data in the file will NOT be loaded ...\n" );
            fprintf( stderr, "          file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
         }

         if ( InputSize < ExpectSize )
         {
            fprintf( stderr, "ERROR : the size of the input file <%s> is too small !!\n", FileName );
            fprintf( stderr, "        Input = %ld bytes <-> Expect = %ld bytes ...\n", InputSize, ExpectSize );
            fprintf( stderr, "        file <%s>, line <%d>, function <%s>\n", __FILE__, __LINE__,  __FUNCTION__  );
            MPI_Exit();
         }
      } // if ( MyRank == 0 )


//    load data from the initial condition file
      for (int YourTurn=0; YourTurn<NGPU; YourTurn++)
      {
         if ( MyRank == YourTurn )
         {
            FILE *File = fopen( FileName, "rb" );

//          properly set the file position indicator for each process
            fseek( File, offset, SEEK_SET );

            for (int i=0; i<N; i++)
            {
               fread( &Mass[i],   sizeof(real), 1, File );

               fread( &Pos[i][0], sizeof(real), 1, File );
               fread( &Pos[i][1], sizeof(real), 1, File );
               fread( &Pos[i][2], sizeof(real), 1, File );

               fread( &Vel[i][0], sizeof(real), 1, File );
               fread( &Vel[i][1], sizeof(real), 1, File );
               fread( &Vel[i][2], sizeof(real), 1, File );
            }

            fclose(File);
         } // if ( MyRank == YourTurn )

         MPI_Barrier( MPI_COMM_WORLD );
      } // for (int YourTurn=0; YourTurn<NGPU; YourTurn++)
   } // if ( INIT_METHOD == INIT_FILE )


// set particle initial condition manually
   else if ( INIT_METHOD == INIT_FUNC )
   {
      Mass[0]    = 0.0;    // massless (doesn't matter since self-gravity is disabled)
      Pos [0][0] = SOL_RSC;
      Pos [0][1] = 0.0;
      Pos [0][2] = 0.0;
      Vel [0][0] = 0.0;
      Vel [0][1] = SQRT(  NEWTON_G*Ext_TotalEnclosedMass( SOL_RSC, INIT_T )/SOL_RSC  );
      Vel [0][2] = 0.0;
   } // else if ( INIT_METHOD == INIT_FUNC )


   else
   {
      fprintf( stderr, "ERROR : unsupported INIT_METHOD (%d) !!\n", INIT_METHOD );
      exit( EXIT_FAILURE );
   }


// set the initial PlayerList = all stars
   L = N;
   for (int i=0; i<N; i++)    PlayerList[i] = i;


// get initial acceleration and jerk
   int Temp_List[NGPU];
   for (int i=0; i<NGPU; i++)    Temp_List[i] = N;

#  ifdef HYBRID_SCHEME
   BeginHybrid_Acc_Jerk( INIT_T, Pos, Vel, Acc, Jerk, Temp_List );
#  else
   BeginRing_Acc_Jerk  ( INIT_T, Pos, Vel, Acc, Jerk, Temp_List, N );
#  endif


   if ( MyRank == 0 )    fprintf( stdout, "Initialize particle data ... done\n" );

}



//----------------------------------------------------------------------
// Function    :  ReadParameter
// Description :  Read parameters from the file "Input__Parameter"
//----------------------------------------------------------------------
void ReadParameter( double &INIT_T, double &END_T, long int &INIT_STEP, long int &END_STEP,
                    double &OUTPUT_DT, double &ENERGY_DT, double &DT_DIAGNOSIS_DT,
                    int &RESTART, real &INIT_E, int &INIT_DUMP_ID, bool &BINARY_OUTPUT, bool &CONST_INIT_DT,
                    int &GPUID_SELECT )
{

   if ( MyRank == 0 )    printf( "Load parameter ...\n" );


   char *input_line;
   char string[100];
   size_t len = 0;
   int Max_dt_Pow, Min_dt_Pow, temp_int;
   double eps;

   FILE *File = fopen( "./Input__Parameter", "r" );

   if ( File == NULL )
   {
     fprintf( stderr, "The file \"Input__Parameter\" does not exit !!\n" );
     exit(-1);
   }


// load parameters
   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &TOTAL_N,         string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &NGPU,            string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &OMP_NTHREAD,     string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &eps,             string );

#  ifdef FLOAT8
   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &NEWTON_G,        string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &ETA,             string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &INIT_ETA,        string );
#  else
   getline(&input_line, &len, File);
   sscanf( input_line, "%f%s",   &NEWTON_G,        string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%f%s",   &ETA,             string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%f%s",   &INIT_ETA,        string );
#  endif

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &END_T,           string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%ld%s",  &END_STEP,        string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &Max_dt_Pow,      string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &Min_dt_Pow,      string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &temp_int,        string );
   CONST_INIT_DT = temp_int;

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &temp_int,        string );
   BINARY_OUTPUT = temp_int;

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &OUTPUT_DT,       string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &ENERGY_DT,       string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &DT_DIAGNOSIS_DT, string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &RESTART,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,        string );
   INIT_METHOD = (OptInitMethod_t)temp_int;

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &INIT_T,          string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%ld%s",  &INIT_STEP,       string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &INIT_DUMP_ID,    string );

#  ifdef FLOAT8
   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &INIT_E,          string );
#  else
   getline(&input_line, &len, File);
   sscanf( input_line, "%f%s",   &INIT_E,          string );
#  endif

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &SPLIT_NMAX,      string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &GPUID_SELECT,    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,        string );
   GRAVITY_TYPE = (OptGravityType_t)temp_int;

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &temp_int,        string );
   EXT_JERK = temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,        string );
   EXT_METHOD = (OptExtMethod_t)temp_int;

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &EXT_SIZE[0],     string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &EXT_SIZE[1],     string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%d%s",   &EXT_SIZE[2],     string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &EXT_CEN[0],      string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &EXT_CEN[1],      string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &EXT_CEN[2],      string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &EXT_DH,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,        string );
   EXT_PAR_INT = (OptExtParInt_t)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,        string );
   EXT_ACC_DER = (OptExtAccDer_t)temp_int;

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &SOL_M22,         string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &SOL_RCORE,       string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &SOL_RSC,         string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &SOL_MSC,         string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &SOL_OSC_AMP,     string );

   getline(&input_line, &len, File);
   sscanf( input_line, "%lf%s",  &SOL_OSC_T,       string );

   fclose(File);


// set default parameters
   if ( NEWTON_G <= 0.0 )
   {
      NEWTON_G = 1.0;

      if ( MyRank == 0 )
         Aux_Message( stdout, "   NEWTON_G is set to the default value = %13.7e\n", NEWTON_G );
   }

   if ( EXT_METHOD == EXT_FILE )
   for (int d=0; d<3; d++)
   {
      if ( EXT_CEN[d] == -9999.9 )
      {
         EXT_CEN[d] = 0.5*EXT_SIZE[d]*EXT_DH;

         if ( MyRank == 0 )
            Aux_Message( stdout, "   EXT_CEN[%d] is set to the default value = %14.7e\n", d, EXT_CEN[d] );
      }
   }

   if (  ( GRAVITY_TYPE == GRAVITY_EXTERNAL || GRAVITY_TYPE == GRAVITY_BOTH )  &&  EXT_METHOD == EXT_FILE  )
   {
      if ( EXT_PAR_INT < 0 )
      {
         EXT_PAR_INT = EXT_PAR_INT_TSC;

         if ( MyRank == 0 )
            Aux_Message( stdout, "   EXT_PAR_INT is set to the default value = %d\n", EXT_PAR_INT );
      }

      if ( EXT_ACC_DER < 0 )
      {
         EXT_ACC_DER = EXT_ACC_DER_QUAR;

         if ( MyRank == 0 )
            Aux_Message( stdout, "   EXT_ACC_DER is set to the default value = %d\n", EXT_ACC_DER );
      }
   }

#  ifdef OPENMP
   const int OMP_Max_NThread = omp_get_max_threads();

   if ( OMP_NTHREAD <= 0 )
   {
      OMP_NTHREAD = OMP_Max_NThread;

      if ( MyRank == 0 )  Aux_Message( stdout, "   OMP_NTHREAD is set to the default value = %d\n", OMP_NTHREAD );
   }

   else if ( OMP_NTHREAD > OMP_Max_NThread   &&  MyRank == 0 )
   {
      Aux_Message( stderr, "WARNING : OMP_NTHREAD (%d) > omp_get_max_threads (%d) !!\n",
                   OMP_NTHREAD, OMP_Max_NThread );
   }

#  else
   if ( OMP_NTHREAD != 1  &&  MyRank == 0 )
      Aux_Message( stderr, "WARNING : OMP_NTHREAD is reset to 1 since \"OPENMP\" is not turned on !!\n" );

   OMP_NTHREAD = 1;
#  endif


// set related parameters
   N       = TOTAL_N / NGPU;
   MAX_DT  = pow( 2.0, (double)Max_dt_Pow );
   MIN_DT  = pow( 2.0, (double)Min_dt_Pow );
   EPS_SQR = real( eps*eps );


   if ( MyRank == 0 )    printf( "Load parameter ... done\n" );


} // FUNCTION : ReadParameter




#ifdef OPENMP
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_OpenMP
// Description :  Initialize OpenMP
//
// Note        :  Please set OMP_NTHREAD in advance for determing the number of OpenMP threads
//-------------------------------------------------------------------------------------------------------
void Init_OpenMP()
{

   if ( MyRank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// numbef of OMP threads
   omp_set_num_threads( OMP_NTHREAD );

// enable/disable nested parallelization
   omp_set_nested( false );

// schedule
   const int chunk_size = 1;
// omp_set_schedule( omp_sched_static,  chunk_size );
// omp_set_schedule( omp_sched_dynamic, chunk_size );
   omp_set_schedule( omp_sched_guided,  chunk_size );
// omp_set_schedule( omp_sched_auto,    chunk_size );


   if ( MyRank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : Init_OpenMP
#endif

