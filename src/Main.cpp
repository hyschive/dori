#include "Dori.h"

// global variables (all declared as extern variables in the file "Global.h")
int              TOTAL_N, N, NGPU, OMP_NTHREAD;
real             NEWTON_G, EPS_SQR, ETA, INIT_ETA;
int              SPLIT_NMAX;
long             LOG_STEP;

OptGravityType_t GRAVITY_TYPE;
OptInitMethod_t  INIT_METHOD;
OptBinaryOrder_t BINARY_ORDER;
OptExtMethod_t   EXT_METHOD;
OptExtParInt_t   EXT_PAR_INT;
OptExtAccDer_t   EXT_ACC_DER;
int              EXT_SIZE[3];
double           EXT_CEN[3], EXT_DH;
bool             EXT_JERK;


double   *t = NULL, *dt = NULL, *dt_BeforeSync = NULL, Global_Time, Next_Global_Time, MAX_DT, MIN_DT;
long int  Step;
int       MyRank, SendRank, RecvRank;
real     *Mass = NULL, (*Pos)[3] = NULL, (*Vel)[3] = NULL, (*Acc)[3] = NULL, (*Jerk)[3] = NULL, *Pot = NULL;
int      *PlayerList = NULL, L;

real    (*Pos_Pred)[3] = NULL, (*Vel_Pred)[3] = NULL, (*Pos_LineUp)[3] = NULL, (*Vel_LineUp)[3] = NULL;
real    (*Acc_local)[3] = NULL, (*Jerk_local)[3] = NULL;

real     *Ext_AccX = NULL, *Ext_AccY = NULL, *Ext_AccZ = NULL;
real     *Ext_dAccX = NULL, *Ext_dAccY = NULL, *Ext_dAccZ = NULL;




//----------------------------------------------------------------------
// Function    :  main
// Description :
//----------------------------------------------------------------------
int main( int argc, char* argv[] )
{

// local variables set by the file "Input__Parameter"
   double   INIT_T, END_T;
   long int INIT_STEP, END_STEP;
   int      INIT_DUMP_ID, GPUID_SELECT;
   double   OUTPUT_DT, ENERGY_DT, MOMENTUM_DT, DT_DIAGNOSIS_DT;
   int      RESTART;
   real     INIT_E    = (real)1.e10;
   real     INIT_L[3] = { (real)1.e10, (real)1.e10, (real)1.e10 };
   bool     BINARY_OUTPUT, CONST_INIT_DT;

   double   Energy_t       = 0.0;
   double   Momentum_t     = 0.0;
   double   Output_t       = 0.0;
   double   dt_diagnosis_t = 0.0;



// initialization
// ======================================================================================================
   Init_MPI( argc, argv );

   ReadParameter( INIT_T, END_T, INIT_STEP, END_STEP, OUTPUT_DT, ENERGY_DT, MOMENTUM_DT, DT_DIAGNOSIS_DT, RESTART,
                  INIT_E, INIT_L, INIT_DUMP_ID, BINARY_OUTPUT, CONST_INIT_DT, GPUID_SELECT );

#  ifdef OPENMP
   Init_OpenMP();
#  endif

   if ( MyRank == 0 )   CheckParameter( INIT_T, END_T, INIT_STEP, END_STEP, OUTPUT_DT, ENERGY_DT, MOMENTUM_DT );


   MemoryAllocate();

   if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH )
   {
      Ext_Init();

      if (  EXT_METHOD == EXT_FILE  )  Ext_LoadExtAcc();
   }

   if ( MyRank == 0 )   TakeNote( INIT_T, END_T, INIT_STEP, END_STEP, ENERGY_DT, MOMENTUM_DT, OUTPUT_DT, DT_DIAGNOSIS_DT,
                                  RESTART, INIT_E, INIT_L, INIT_DUMP_ID, BINARY_OUTPUT, CONST_INIT_DT, GPUID_SELECT );

#  ifdef GPU
   CUAPI_SetDevice( GPUID_SELECT );
   CUAPI_DiagnoseDevice();
   CUAPI_MemAllocate();
#  endif

   Init_Particles( INIT_T );
   Init_t_dt_step( INIT_T, INIT_STEP, Energy_t, Momentum_t, Output_t, dt_diagnosis_t, ENERGY_DT, MOMENTUM_DT, OUTPUT_DT, DT_DIAGNOSIS_DT,
                   CONST_INIT_DT );

   Get_Next_Global_Time();


   if ( ENERGY_DT >= 0.0 )    Get_TotalEnergy( RESTART, INIT_E );

   if ( MOMENTUM_DT >= 0.0 )  Get_TotalMomentum( RESTART, INIT_L );

   // do not output the initial condition when doing Restart
   if ( OUTPUT_DT >= 0.0 && !RESTART )    OutputData( INIT_DUMP_ID, BINARY_OUTPUT );

   if ( DT_DIAGNOSIS_DT >= 0.0 )    dt_diagnosis();


   double Timer_MPI = MPI_Wtime();  // MPI timer

// main loop
// ======================================================================================================
   while ( (Global_Time-END_T < -1.e-10) && (Step < END_STEP) )
   {

      if ( (MyRank == 0) && (Step%LOG_STEP == 0) )
         printf( "Time : %20.14e -> %20.14e,   Step = %9ld -> %9ld,   dt = %20.14e\n",
                 Global_Time, Next_Global_Time, Step, Step+1, Next_Global_Time-Global_Time );


//    determine the group of particles to be advanced to the Next_Global_Time
//    (for the shared time-step, we only have to set up the line-up starts list once during the initialization)
#     ifndef SHARED_TIMESTEP
      Get_LineUpStars();
#     endif


//    evolve line-up particles, set the present global time, get the new time-step and Next_Global_Time
      KickStars( false );

      Step ++;
      Global_Time = Next_Global_Time;
      Get_Next_Global_Time();


      if ( ENERGY_DT >= 0.0  &&  Global_Time >= Energy_t )
      {
         Synchronize();
         Get_TotalEnergy( false, 0.0 );
         Energy_t += ENERGY_DT;
      }

      if ( MOMENTUM_DT >= 0.0  &&  Global_Time >= Momentum_t )
      {
         Synchronize();
         Get_TotalMomentum( false, NULL );
         Momentum_t += MOMENTUM_DT;
      }

      if ( OUTPUT_DT >= 0.0  &&  Global_Time >= Output_t )
      {
          Synchronize();
          OutputData( INIT_DUMP_ID, BINARY_OUTPUT );
          Output_t += OUTPUT_DT;
      }

      if ( DT_DIAGNOSIS_DT >= 0.0  &&  Global_Time >= dt_diagnosis_t )
      {
         dt_diagnosis();
         dt_diagnosis_t += DT_DIAGNOSIS_DT;
      }

   } // while ( (Global_Time-END_T < -1.e-10) && (Step < END_STEP) )

   Synchronize();


// record the total simulation elapsed time
   MPI_Barrier( MPI_COMM_WORLD );
   Timer_MPI = MPI_Wtime() - Timer_MPI;


   if ( ENERGY_DT >= 0.0 )    Get_TotalEnergy( false, 0.0 );

   if ( MOMENTUM_DT >= 0.0 )  Get_TotalMomentum( false, NULL );

   if ( OUTPUT_DT >= 0.0 )    OutputData( INIT_DUMP_ID, BINARY_OUTPUT );

   if ( DT_DIAGNOSIS_DT >= 0.0 )    dt_diagnosis();


   if ( MyRank == 0 )
   {
      FILE *Note = fopen( "Record__Note", "a" );

      fprintf( Note, "\n" );
      fprintf( Note, "Total Steps                      : %ld\n" , Step );
      fprintf( Note, "Total Processing Time (by MPI)   : %lf s\n", Timer_MPI );
      fprintf( Note, "Time per Step         (by MPI)   : %lf s\n", Timer_MPI/Step );
      fprintf( Note, "\n");

      fclose( Note );
   }


   MemoryFree();

#  ifdef GPU
   CUAPI_MemFree();
#  endif

   MPI_Finalize();
   if ( MyRank == 0 )   fprintf( stdout, "Program terminated successfully\n" );
   exit(0);

}



