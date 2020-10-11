#include "Dori.h"

// global variables for Evolve subroutine
static real (* Acc_Pred)[3];
static real (*Jerk_Pred)[3];

static double (*a2_dt2)[3];
static double (*a3_dt3)[3];

#ifdef IMPLICIT_HERMITE
static real (*Vel0)[3];
#endif


// global variables for Ring scheme
static real  *Pot_local;

static real (*Other_Mass);
static real (*Other_Pos)[3];

static real (*SendBuffer_PosVel) [6];
static real (*RecvBuffer_PosVel) [6];
static real (*SendBuffer_AccJerk)[6];
static real (*RecvBuffer_AccJerk)[6];
static real (*SendBuffer_MassPos)[4];
static real (*RecvBuffer_MassPos)[4];

// global variables for BeginHybrid_Acc_Jerk subroutine
#ifdef HYBRID_SCHEME
static real (*Pos_LineUp_local) [3];
static real (*Vel_LineUp_local) [3];
#endif



//----------------------------------------------------------------------
// Function    :  Synchronize
// Description :  Advance all particles to the present Global_Time
//----------------------------------------------------------------------
void Synchronize()
{

// nothing to do for the share time-step scheme
#  ifdef SHARED_TIMESTEP
   return;
#  endif


// determine the group of particles to be advanced and their corresponding time-steps
   L = 0;

   for (int i=0; i<N; i++)
   {
      if ( t[i] != Global_Time )
      {
         PlayerList[L] = i;
         L++;
      }
   }


// perform synchronization if there is any rank with line-up particles
   int L_Max;
   MPI_Allreduce( &L, &L_Max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );

   if ( L_Max > 0 )
   {
//    set "Next_Global_Time = Global_Time" since "KickStars" will advance particles to the "Next_Global_Time"
//    store the Next_Global_Time before doning synchronization
      const double Next_Global_Time_BeforeSync = Next_Global_Time;
      Next_Global_Time = Global_Time;

//    advance all line-up particles to the present Global_Time
      KickStars( true );

//    restore the Next_Global_Time
      Next_Global_Time = Next_Global_Time_BeforeSync;
    }

} // FUNCTION : Synchronize



//----------------------------------------------------------------------
// Function    :  Get_Next_Global_Time
// Description :  Calculate the Next_Global_Time from the current PlayerList
//----------------------------------------------------------------------
void Get_Next_Global_Time()
{

   Next_Global_Time = 1.0e10;

   for (int j=0; j<L; j++)
   {
      int i = PlayerList[j];
      Next_Global_Time = fmin( t[i]+dt[i], Next_Global_Time );
   }

   Next_Global_Time = Get_Minimum_FromAllRanks( Next_Global_Time );

} // FUNCTION : Get_Next_Global_Time



//----------------------------------------------------------------------
// Function    :  Get_LineUpStars
// Description :  Determine the group of particles to be advanced to the Next_Global_Time
//----------------------------------------------------------------------
void Get_LineUpStars()
{
   L = 0;

   for (int i=0; i<N; i++)
   {
      if (  fabs( (t[i]+dt[i]-Next_Global_Time)/Next_Global_Time ) <= 1.e-10  )
      {
         PlayerList[L] = i;
         L++;
      }
   }
} // FUNCTION : Get_LineUpStars


//----------------------------------------------------------------------
// Function    :  Get_NewTimeStep
// Description :  Calculate the new time-step by Aarseth's formula
//----------------------------------------------------------------------
void Get_NewTimeStep( const int i, const int j, const double ds )
{

   double a0_sqr, a1_sqr, a2_dt2_sqr, a3_dt3_sqr, Next_dt;

   for (int dim=0; dim<3; dim++)   a2_dt2[j][dim] += a3_dt3[j][dim];

   a0_sqr =  Acc[i][0]* Acc[i][0] +  Acc[i][1]* Acc[i][1] +  Acc[i][2]* Acc[i][2];
   a1_sqr = Jerk[i][0]*Jerk[i][0] + Jerk[i][1]*Jerk[i][1] + Jerk[i][2]*Jerk[i][2];
   a2_dt2_sqr = a2_dt2[j][0]*a2_dt2[j][0] + a2_dt2[j][1]*a2_dt2[j][1] + a2_dt2[j][2]*a2_dt2[j][2];
   a3_dt3_sqr = a3_dt3[j][0]*a3_dt3[j][0] + a3_dt3[j][1]*a3_dt3[j][1] + a3_dt3[j][2]*a3_dt3[j][2];

   Next_dt = ETA * ds * sqrt(   ( sqrt( a0_sqr*a2_dt2_sqr ) + a1_sqr*ds*ds )  /
                                ( ds*sqrt( a1_sqr*a3_dt3_sqr ) + a2_dt2_sqr )   );

// acceleration may be zero if only external acceleration is considered
   if ( ! Aux_IsFinite(Next_dt) )
   {
      Aux_Message( stderr, "WARNING : strange time-step is found (Next_dt = %13.7e)\n", Next_dt );
      Aux_Message( stderr, "          ds = %13.7e, a0_sqr = %13.7e, a1_sqr = %13.7e, a2_dt2_sqr = %13.7e, a3_dt3_sqr = %13.7e\n",
                   ds, a0_sqr, a1_sqr, a2_dt2_sqr, a3_dt3_sqr );
      Aux_Message( stderr, "          --> set dt = MAX_DT (%13.7e) !!\n", MAX_DT );

      Next_dt = MAX_DT;
   }

#  ifdef SHARED_TIMESTEP

   dt[i] = Next_dt;

#  else // INDIVIDUAL_TIMESTEP

   double ds_pow2, ds2_pow2, Pow;

   if ( Next_dt < ds )
   {
      Pow   = floor( log2(Next_dt) );
      dt[i] = pow( 2.0, Pow );
   }
   else
   {
      ds_pow2  = pow(  2.0, floor( log2(ds) )  );  // in case of synchronization
      ds2_pow2 = 2.0*ds_pow2;

      if (  ( Next_dt >= ds2_pow2 )  &&  ( fmod(Next_Global_Time, ds2_pow2) == 0.0 )  )
         dt[i] = ds2_pow2;
      else
         dt[i] = ds_pow2;
   }

   if ( dt[i] < MIN_DT )   dt[i] = MIN_DT;
   if ( dt[i] > MAX_DT )   dt[i] = MAX_DT;

#  endif

} // FUNCTION : Get_NewTimeStep



//----------------------------------------------------------------------
// Function    :  BeginRing_Acc_Jerk
// Description :  Use the Ring Scheme to accumulate the acceleration and jerk from all other particles
//----------------------------------------------------------------------
void BeginRing_Acc_Jerk( const real Pos_Pred[][3], const real Vel_Pred[][3], real Acc_Pred[][3],
                         real Jerk_Pred[][3], const int L_List[], const int L_Max )
{


// 1. initialize the XXX_LineUp arrays which are used to store the information of particles being advanced
//    initialize the (Acc,Jerk)_Pred arrays as zero
#  pragma omp parallel for schedule( runtime )
   for (int j=0; j<L; j++)
   {
      const int i = PlayerList[j];

      for (int dim=0; dim<3; dim++)
      {
          Pos_LineUp[j][dim] = Pos_Pred[i][dim];
          Vel_LineUp[j][dim] = Vel_Pred[i][dim];

           Acc_Pred [j][dim] = (real)0.0;
          Jerk_Pred [j][dim] = (real)0.0;
      }
   }


// 2. self-gravity
   if ( GRAVITY_TYPE == GRAVITY_SELF  ||  GRAVITY_TYPE == GRAVITY_BOTH )
   {
      SendBuffer_PosVel  = new real[L_Max][6];
      RecvBuffer_PosVel  = new real[L_Max][6];
      SendBuffer_AccJerk = new real[L_Max][6];
      RecvBuffer_AccJerk = new real[L_Max][6];

      MPI_Request Req_PosVel [2];
      MPI_Request Req_AccJerk[2];

      int ProcessID;          // the part of data being processed (0, 1, ..., NGPU-1)
      int L_Process     = 0;  // number of line-up particles being processed in the current stage
      int Nsend_PosVel;       // number of particles whose position and velocity are being sent
      int Nrecv_PosVel;       // number of particles whose position and velocity are being received
      int Nsend_AccJerk = 0;  // number of particles whose acceleration and jerk are being sent
      int Nrecv_AccJerk = 0;  // number of particles whose acceleration and jerk are being received

      bool Copy         = true;


//    begin RING scheme
      for (int RingLoop=0; RingLoop<NGPU; RingLoop++)
      {
         ProcessID     = (MyRank+RingLoop) % NGPU;
         L_Process     = L_List[ ProcessID ];
         Nsend_PosVel  = L_Process;
         Nrecv_PosVel  = L_List[ (ProcessID+1) % NGPU ];
         Nsend_AccJerk = L_List[ (ProcessID+NGPU-1) % NGPU ];
         Nrecv_AccJerk = L_Process;

//       send data to SendRank and receive data from RecvRank
         if ( NGPU != 1 )
         {
            if ( RingLoop != 0 )       DataTransfer_Force ( Acc_Pred, Jerk_Pred, SendBuffer_AccJerk,
                                                            RecvBuffer_AccJerk, Req_AccJerk, Nsend_AccJerk,
                                                            Nrecv_AccJerk );
            if ( RingLoop != NGPU-1 )  DataTransfer_Force ( Pos_LineUp, Vel_LineUp, SendBuffer_PosVel,
                                                            RecvBuffer_PosVel, Req_PosVel,
                                                            Nsend_PosVel, Nrecv_PosVel );
         }

//       calculate acceleration and jerk of line-up particles and save them to Acc_local and Jerk_local
//--------------------------------------------------------------------------------------------------------------
         if ( L_Process != 0 )
         {
#           ifdef GPU
            CUAPI_Acc_Jerk( N, Mass, Pos_Pred, Vel_Pred, L_Process, Pos_LineUp, Vel_LineUp,
                            Acc_local, Jerk_local, EPS_SQR, Copy, NEWTON_G );
#           else
            CPU_Acc_Jerk  ( N, Mass, Pos_Pred, Vel_Pred, L_Process, Pos_LineUp, Vel_LineUp,
                            Acc_local, Jerk_local, EPS_SQR, NEWTON_G );
#           endif
            Copy = false;
         }
//--------------------------------------------------------------------------------------------------------------


//       wait until the DataTransfer_Force function is complete and copy the received data properly
         if ( NGPU != 1 )
         {
            if ( RingLoop != 0 )          Block_and_Copy__Force ( Acc_Pred, Jerk_Pred, RecvBuffer_AccJerk,
                                                                  Req_AccJerk, Nrecv_AccJerk );
            if ( RingLoop != NGPU-1 )     Block_and_Copy__Force ( Pos_LineUp, Vel_LineUp, RecvBuffer_PosVel,
                                                                  Req_PosVel, Nrecv_PosVel );
         }


//       add the data of "Acc_local and Jerk_local" to "Acc_Pred and Jerk_Pred"
         for (int j=0; j<L_Process; j++)
         {
            for (int dim=0; dim<3; dim++)
            {
                Acc_Pred[j][dim] +=  Acc_local[j][dim];
               Jerk_Pred[j][dim] += Jerk_local[j][dim];
            }
         }

       } // for (int RingLoop=0; RingLoop<NGPU; RingLoop++)

//    reset the number of particles to be sent and received for the final data trasnfer
      if ( NGPU != 1 )
      {
         Nsend_AccJerk = L_Process;
         Nrecv_AccJerk = L;

         DataTransfer_Force   ( Acc_Pred, Jerk_Pred, SendBuffer_AccJerk, RecvBuffer_AccJerk, Req_AccJerk,
                                Nsend_AccJerk, Nrecv_AccJerk );
         Block_and_Copy__Force( Acc_Pred, Jerk_Pred, RecvBuffer_AccJerk, Req_AccJerk, Nrecv_AccJerk );
      }

      delete [] SendBuffer_PosVel;
      delete [] RecvBuffer_PosVel;
      delete [] SendBuffer_AccJerk;
      delete [] RecvBuffer_AccJerk;
   } // if ( GRAVITY_TYPE == GRAVITY_SELF  ||  GRAVITY_TYPE == GRAVITY_BOTH )


// 3. external gravity
   if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH )
   {
      switch ( EXT_METHOD )
      {
         case EXT_FUNC :   Ext_AddAccFromFunc( L, Pos_LineUp, Vel_LineUp, Acc_Pred, Jerk_Pred );  break;
         case EXT_FILE :   Ext_AddAccFromFile( L, Pos_LineUp, Vel_LineUp, Acc_Pred, Jerk_Pred );  break;
         default       :   Aux_Error( ERROR_INFO, "unsupported EXT_METHOD (%d) !!\n", EXT_METHOD );
      }
   } // if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH )

} // FUNCTION : BeginRing_Acc_Jerk



//----------------------------------------------------------------------
// Function    :  BeginRing_Pot
// Description :  Use the Ring Scheme to accumulate the potential from all other particles
//----------------------------------------------------------------------
void BeginRing_Pot()
{

   Pot_local   = new real[N];

   Other_Mass  = new real[N];
   Other_Pos   = new real[N][3];

   SendBuffer_MassPos = new real[N][4];
   RecvBuffer_MassPos = new real[N][4];

   const int MemSize = N*sizeof(real);
   MPI_Request Req[2];


// initialize Pot arrays as zero
   for (int i=0; i<N; i++)    Pot[i] = (real)0.0;


// copy the data of its own to "Other_XXX" arrays to calculate the potential from  group
   memcpy( Other_Mass, Mass,  MemSize);
   memcpy( Other_Pos,  Pos, 3*MemSize);


// begin RING scheme
   for (int RingLoop=0; RingLoop<NGPU; RingLoop++)
   {

//    transfer the "Other_XXX" arrays to SecvRank and receive data from RecvRank
      if ( (NGPU != 1) && (RingLoop != NGPU-1) )
         DataTransfer_Pot( Other_Mass, Other_Pos, SendBuffer_MassPos, RecvBuffer_MassPos, Req );

//    calculate potential exerted by "Other" particles and save them to Pot_local
#     ifdef GPU
      CUAPI_Pot( N, Other_Mass, Other_Pos, N, Pos, Pot_local, EPS_SQR, NEWTON_G );
#     else
      CPU_Pot  ( N, Other_Mass, Other_Pos, N, Pos, Pot_local, EPS_SQR, NEWTON_G );
#     endif

//    add the data of "Pot_local" to "Pot"
      for (int i=0; i<N; i++)    Pot[i] += Pot_local[i];

//    wait until the DataTransfer function is complete and copy the received data to "Other_XXX" arrays
      if ( (NGPU != 1) && (RingLoop != NGPU-1) )
         Block_and_Copy__Pot( Other_Mass, Other_Pos, RecvBuffer_MassPos, Req );
   }

   delete [] Pot_local;

   delete [] Other_Mass;
   delete [] Other_Pos;

   delete [] SendBuffer_MassPos;
   delete [] RecvBuffer_MassPos;

} // FUNCTION : BeginRing_Pot



//---------------------------------------------------------------------------------------------------
// Function    :  KickStars
// Description :  Evolve all line-up particles to Next_Global_Time and calculate the new time-step
//---------------------------------------------------------------------------------------------------
void KickStars( const bool Sync )
{

   const double temp1 = 1.0/6.0;
   int L_List[NGPU];
   int L_Max = 0;

// get the number of line-up particles of different ranks
   Get_L_FromAllRanks( L_List );

// get the maximum number among L_List
   for (int i=0; i<NGPU; i++)    if ( L_List[i] > L_Max )   L_Max = L_List[i];

     Acc_Pred = new real [L_Max][3];
    Jerk_Pred = new real [L_Max][3];


// Predictor
//=====================================================================
#  pragma omp parallel for schedule( runtime )
   for (int i=0; i<N; i++)
   {
      const double ds = Next_Global_Time - t[i];

      for (int dim=0; dim<3; dim++)
      {
         Pos_Pred[i][dim] = Pos[i][dim] + Vel[i][dim]*ds + (real)0.5*Acc [i][dim]*ds*ds + temp1*Jerk[i][dim]*ds*ds*ds;
         Vel_Pred[i][dim] = Vel[i][dim] + Acc[i][dim]*ds + (real)0.5*Jerk[i][dim]*ds*ds;
      }
   }


// get the acceleration and jerk from predicted postion and velocity
//=====================================================================
#  ifdef HYBRID_SCHEME
   BeginHybrid_Acc_Jerk ( Pos_Pred, Vel_Pred, Acc_Pred, Jerk_Pred, L_List );
#  else
   BeginRing_Acc_Jerk   ( Pos_Pred, Vel_Pred, Acc_Pred, Jerk_Pred, L_List, L_Max );
#  endif


// Corrector
//=====================================================================
// (second,third) derivative of acceleration multiplied by dt^(2,3)
   a2_dt2 = new double [L][3];
   a3_dt3 = new double [L][3];

#  ifdef IMPLICIT_HERMITE
   const int MemSize = 3*N*sizeof(real);
   Vel0 = new real [N][3];
   memcpy( Vel0, Vel, MemSize );
#  endif


#  pragma omp parallel for schedule( runtime )
   for (int j=0; j<L; j++)
   {
      const int i     = PlayerList[j];
      const double ds = Next_Global_Time - t[i];

      for (int dim=0; dim<3; dim++)
      {
         a2_dt2[j][dim] =      -6.0*(      Acc[i][dim] -      Acc_Pred[j][dim] )
                           -     ds*( 4.0*Jerk[i][dim] + 2.0*Jerk_Pred[j][dim] );
         a3_dt3[j][dim] =      12.0*(      Acc[i][dim] -      Acc_Pred[j][dim] )
                           + 6.0*ds*(     Jerk[i][dim] +     Jerk_Pred[j][dim] );

#        ifdef IMPLICIT_HERMITE
         const double temp2 = 1.0/12.0;

         Vel[i][dim] +=     0.5*(  Acc[i][dim] +  Acc_Pred[j][dim] )*ds
                        + temp2*( Jerk[i][dim] - Jerk_Pred[j][dim] )*ds*ds;

         Pos[i][dim] +=     0.5*( Vel0[i][dim] +       Vel[i][dim] )*ds
                        + temp2*(  Acc[i][dim] -  Acc_Pred[j][dim] )*ds*ds;
#        else

         Pos[i][dim] = Pos_Pred[i][dim] + temp1*( 0.25*a2_dt2[j][dim] + 0.05*a3_dt3[j][dim] )*ds*ds;
         Vel[i][dim] = Vel_Pred[i][dim] + temp1*(      a2_dt2[j][dim] + 0.25*a3_dt3[j][dim] )*ds;
#        endif

//       update both the Acc and Jerk to the present time
         Acc[i][dim] =  Acc_Pred[j][dim];
        Jerk[i][dim] = Jerk_Pred[j][dim];
      } // for (int dim=0; dim<3; dim++)

//    update the physical time for all line-up particles
      t[i] = Next_Global_Time;

//    get the new time-step
      if ( Sync )    dt[i] -= ds;
      else           Get_NewTimeStep( i, j, ds );
   } // for (int j=0; j<L; j++)


// get the new time-step for the shared time-step scheme
#  ifdef SHARED_TIMESTEP
   double Min_dt = 1.e10;
   for (int i=0; i<N; i++)    if ( dt[i] < Min_dt )   Min_dt = dt[i];
   Min_dt = Get_Minimum_FromAllRanks( Min_dt );
   for (int i=0; i<N; i++)    dt[i] = Min_dt;
#  endif


   delete []  Acc_Pred;
   delete [] Jerk_Pred;

   delete [] a2_dt2;
   delete [] a3_dt3;

#  ifdef IMPLICIT_HERMITE
   delete [] Vel0;
#  endif

} // FUNCTION : KickStars



#ifdef HYBRID_SCHEME
//----------------------------------------------------------------------
// Function    :  BeginHybrid_Acc_Jerk
// Description :  Calculate the acceleration and jerk by hybrid scheme (Harfst, 2007)
//                --> The alterantive subroutine to BeginRing_Acc_Jerk
//----------------------------------------------------------------------
void BeginHybrid_Acc_Jerk( const real Pos_Pred[][3], const real Vel_Pred[][3], real Acc_Pred[][3],
                           real Jerk_Pred[][3], const int L_List[] )
{
   int RecvCounts_List[NGPU], Disp_List[NGPU];
   int L_Sum = 0;

   Pos_LineUp_local = new real [L][3];
   Vel_LineUp_local = new real [L][3];


// 1. initializae the XX_LineUp_local arrays
#  pragma omp parallel for schedule( runtime )
   for (int j=0; j<L; j++)
   {
      const int i = PlayerList[j];

      for (int dim=0; dim<3; dim++)
      {
          Pos_LineUp_local[j][dim] = Pos_Pred[i][dim];
          Vel_LineUp_local[j][dim] = Vel_Pred[i][dim];

           Acc_Pred       [j][dim] = (real)0.0;    // necessary when only external acceleration is considered
          Jerk_Pred       [j][dim] = (real)0.0;    // ...
      }
   }


// 2. self-gravity
   if ( GRAVITY_TYPE == GRAVITY_SELF  ||  GRAVITY_TYPE == GRAVITY_BOTH )
   {
//    get the sum of L_List
      for (int i=0; i<NGPU; i++)    L_Sum += L_List[i];

//    get the RecvCounts list
      for (int i=0; i<NGPU; i++)    RecvCounts_List[i] = 3*L_List[i];

//    get the displacement list
      Disp_List[0] = 0;
      for (int i=1; i<NGPU; i++)    Disp_List[i] = Disp_List[i-1] + RecvCounts_List[i-1];


//    gather the data of line-up particles from all ranks
      MPI_Allgatherv( Pos_LineUp_local, 3*L, GAMER_MPI_REAL, Pos_LineUp, RecvCounts_List, Disp_List, GAMER_MPI_REAL,
                      MPI_COMM_WORLD );
      MPI_Allgatherv( Vel_LineUp_local, 3*L, GAMER_MPI_REAL, Vel_LineUp, RecvCounts_List, Disp_List, GAMER_MPI_REAL,
                      MPI_COMM_WORLD );

      /*
      AllGather( Pos_LineUp_local, Pos_LineUp, RecvCounts_List, Disp_List );
      AllGather( Vel_LineUp_local, Vel_LineUp, RecvCounts_List, Disp_List );
      */


//    calculate the local acceleration and jerk by GPU
//---------------------------------------------------------------------------------------------------------
#     ifdef GPU
      CUAPI_Acc_Jerk( N, Mass, Pos_Pred, Vel_Pred, L_Sum, Pos_LineUp, Vel_LineUp, Acc_local, Jerk_local,
                      EPS_SQR, true, NEWTON_G );
#     else
      CPU_Acc_Jerk  ( N, Mass, Pos_Pred, Vel_Pred, L_Sum, Pos_LineUp, Vel_LineUp, Acc_local, Jerk_local,
                      EPS_SQR, NEWTON_G );
#     endif
//---------------------------------------------------------------------------------------------------------


//    get the total accleration and jerk by summing up all local values
      AllReduce(  Acc_local,  Acc_Pred, RecvCounts_List, Disp_List );
      AllReduce( Jerk_local, Jerk_Pred, RecvCounts_List, Disp_List );

      /*
      for (int MyTurn=0; MyTurn<NGPU; MyTurn++)
      {
         MPI_Reduce(  Acc_local[0]+Disp_List[MyTurn],  Acc_Pred, RecvCounts_List[MyTurn], GAMER_MPI_REAL, MPI_SUM,
                     MyTurn, MPI_COMM_WORLD );
         MPI_Reduce( Jerk_local[0]+Disp_List[MyTurn], Jerk_Pred, RecvCounts_List[MyTurn], GAMER_MPI_REAL, MPI_SUM,
                     MyTurn, MPI_COMM_WORLD );
      }
      */
   } // if ( GRAVITY_TYPE == GRAVITY_SELF  ||  GRAVITY_TYPE == GRAVITY_BOTH )


// 3. external gravity
   if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH )
   {
      switch ( EXT_METHOD )
      {
         case EXT_FUNC :   Ext_AddAccFromFunc( L, Pos_LineUp_local, Vel_LineUp_local, Acc_Pred, Jerk_Pred );  break;
         case EXT_FILE :   Ext_AddAccFromFile( L, Pos_LineUp_local, Vel_LineUp_local, Acc_Pred, Jerk_Pred );  break;
         default       :   Aux_Error( ERROR_INFO, "unsupported EXT_METHOD (%d) !!\n", EXT_METHOD );
      }
   } // if ( GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  GRAVITY_TYPE == GRAVITY_BOTH )


   delete [] Pos_LineUp_local;
   delete [] Vel_LineUp_local;

} // FUNCTION : BeginHybrid_Acc_Jerk
#endif // #ifdef HYBRID_SCHEME

