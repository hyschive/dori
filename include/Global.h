#ifndef __GLOBAL_VARIABLES_H__
#define __GLOBAL_VARIABLES_H__



// global variables set by "Parameter" file
extern int              TOTAL_N;
extern int              N;
extern int              NGPU;
extern int              OMP_NTHREAD;
extern real             NEWTON_G;
extern real             EPS_SQR;
extern real             ETA;
extern real             INIT_ETA;
extern double           MAX_DT;
extern double           MIN_DT;

extern OptGravityType_t GRAVITY_TYPE;
extern OptInitMethod_t  INIT_METHOD;
extern OptExtMethod_t   EXT_METHOD;
extern OptExtParInt_t   EXT_PAR_INT;
extern OptExtAccDer_t   EXT_ACC_DER;
extern int              EXT_SIZE[3];
extern double           EXT_CEN[3], EXT_DH;
extern bool             EXT_JERK;


// global time and step
extern double     Global_Time;
extern double     Next_Global_Time;
extern long int   Step;


// rank information
extern int MyRank, SendRank, RecvRank;


// particle data
extern real  *Mass;
extern real (*Pos) [3];
extern real (*Vel) [3];
extern real (*Acc) [3];
extern real (*Jerk)[3];
extern real  *Pot;

extern double *t;
extern double *dt;
extern double *dt_BeforeSync;

extern real (*Pos_Pred  )[3];
extern real (*Vel_Pred  )[3];
extern real (*Pos_LineUp)[3];
extern real (*Vel_LineUp)[3];
extern real (* Acc_local)[3];
extern real (*Jerk_local)[3];

extern real  *Ext_AccX;
extern real  *Ext_AccY;
extern real  *Ext_AccZ;
extern real  *Ext_dAccX;
extern real  *Ext_dAccY;
extern real  *Ext_dAccZ;


// information of the particles to be advanced
extern int *PlayerList;
extern int  L;    // length of the PlayerList


// parameters of splitting
extern int SPLIT_NMAX;


// for star cluster heating
extern double SOL_M22;
extern double SOL_RCORE;
extern double SOL_RSC;
extern double SOL_MSC;
extern double SOL_OSC_AMP;
extern double SOL_OSC_T;

extern double SOL_DENS;
extern double SOL_TSC;



#endif // #ifndef __GLOBAL_VARIABLES_H__


