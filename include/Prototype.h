#ifndef __PROTOTYPE_H__
#define __PROTOTYPE_H__



// Auxiliary.h
void TakeNote( const double INIT_T, const double END_T, const long int INIT_STEP, const long int END_STEP,
               const double ENERGY_DT, const double OUTPUT_DT, const double DT_DIAGNOSIS_DT,
               const int RESTART, const real INIT_E, const int INIT_DUMP_ID, const bool BINARY_OUTPUT,
               const bool CONST_INIT_DT, const int GPUID_SELECT );
void Get_TotalEnergy( bool UseInputEgy, real INIT_E );
void OutputData( const int Init_DumpID, const bool Binary_Output );
void MemoryAllocate();
void MemoryFree();
void dt_diagnosis();
void CheckParameter( const double INIT_T, const double END_T, const long int INIT_STEP, const long int END_STEP,
                     const double OUTPUT_DT, const double ENERGY_DT );
void GetCPUInfo( const char *FileName );
void Aux_Message( FILE *Type, const char *Format, ... );
void Aux_Error( const char *File, const int Line, const char *Func, const char *Format, ... );
void MPI_Exit();
bool Aux_CheckFileExist( const char *FileName );
int Aux_IsFinite( const float x );
int Aux_IsFinite( const double x );


// CUDA_API.h
void CUAPI_Acc_Jerk( const int Nj, const real hJ_Mass[], const real hJ_Pos[][3], const real hJ_Vel[][3],
                     const int Ni, const real hI_Pos[][3], const real hI_Vel[][3],
                     real hI_Acc[][3], real hI_Jerk[][3], const real Eps2, const bool Copy_J_MassPosVel,
                     const real G );
void CUAPI_Pot( const int Nj, const real hJ_Mass[], const real hJ_Pos[][3],
                const int Ni, const real hI_Pos[][3], real hI_Pot[], const real Eps2, const real G );
void CUAPI_SetDevice( const int Mode );
void CUAPI_DiagnoseDevice();
void CUAPI_MemAllocate();
void CUAPI_MemFree();


// DataTransfer.h
void DataTransfer_Force( const real Data1[][3], const real Data2[][3],
                         real SendBuffer[][6], real RecvBuffer[][6],
                         MPI_Request Req[], const int Nsend, const int Nrecv );
void DataTransfer_Pot( const real Other_Mass[], const real Other_Pos[][3],
                       real SendBuffer[][4], real RecvBuffer[][4], MPI_Request Req[] );
void Block_and_Copy__Force( real Data1[][3], real Data2[][3], const real RecvBuffer[][6],
                            MPI_Request Req[], const int Nrecv );
void Block_and_Copy__Pot( real Other_Mass[], real Other_Pos[][3], const real RecvBuffer[][4],
                          MPI_Request Req[] );
void Get_L_FromAllRanks( int L_List[] );
double Get_Minimum_FromAllRanks( double input );
#ifdef HYBRID_SCHEME
void AllGather( real Array_local[][3], real Array_total[][3], const int ArraySize[], const int Disp[] );
void AllReduce( real Array_partial[][3], real Array_sum[][3], const int ArraySize[], const int Disp[] );
#endif


// Evolve.h
void Get_NewTimeStep( const int i, const int j, const double ds );
void BeginRing_Acc_Jerk( const real Pos_Pred[][3], const real Vel_Pred[][3], real Acc_Pred[][3],
                         real Jerk_Pred[][3], const int L_List[], const int L_Max );
void BeginRing_Pot();
void KickStars( const bool Sync );
void Get_LineUpStars();
void Get_Next_Global_Time();
void Get_LineUpStars();
void Synchronize();
#ifdef HYBRID_SCHEME
void BeginHybrid_Acc_Jerk( const real Pos_Pred[][3], const real Vel_Pred[][3], real Acc_Pred[][3],
                           real Jerk_Pred[][3], const int L_List[] );
#endif


// Init.h
void Init_MPI( int argc, char *argv[] );
void Init_t_dt_step( const double INIT_T, const long int INIT_STEP, double &Energy_t, double &Output_t,
                     double &dt_diagnosis_t, const double ENERGY_DT, const double OUTPUT_DT,
                     const double DT_DIAGNOSIS_DT, const bool CONST_INIT_DT );
void Init_Particles();
void ReadParameter( double &INIT_T, double &END_T, long int &INIT_STEP, long int &END_STEP,
                    double &OUTPUT_DT, double &ENERGY_DT, double &DT_DIAGNOSIS_DT,
                    int &RESTART, real &INIT_E, int &INIT_DUMP_ID, bool &BINARY_OUTPUT, bool &CONST_INIT_DT,
                    int &GPUID_SELECT );
#ifdef OPENMP
void Init_OpenMP();
#endif


// CPU solvers
void CPU_Acc_Jerk( const int Nj, const real J_Mass[], const real J_Pos[][3], const real J_Vel[][3],
                   const int Ni, const real I_Pos[][3], const real I_Vel[][3],
                   real I_Acc[][3], real I_Jerk[][3], const real Eps2, const real G );
void CPU_Pot( const int Nj, real J_Mass[], real J_Pos[][3], const int Ni, real I_Pos[][3], real I_Pot[],
              const real Eps2, const real G );


// ExternalForce.h
void Ext_AddAccFromFunc( const int NPar, const real (*MyPos)[3], const real (*MyVel)[3],
                         real (*MyAcc)[3], real (*MyJerk)[3] );
void Ext_AddAccFromFile( const int NPar, const real (*MyPos)[3], const real (*MyVel)[3],
                         real (*MyAcc)[3], real (*MyJerk)[3] );
void Ext_LoadExtAcc();



#endif // #ifndef __PROTOTYPE_H__


