#include "Dori.h"

#ifdef HYBRID_SCHEME
static real (*Buffer)[3];
#endif




//---------------------------------------------------------------------------
// Function    :  DataTransfer_Force
// Description :  Transfer data between different processes when calculating force
//---------------------------------------------------------------------------
void DataTransfer_Force( const real Data1[][3], const real Data2[][3],
                         real SendBuffer[][6], real RecvBuffer[][6],
                         MPI_Request Req[], const int Nsend, const int Nrecv )
{

// copy data from Data1 and Data2 to SendBuffer
   for (int i=0; i<Nsend; i++)
   {
      SendBuffer[i][0] = Data1[i][0];
      SendBuffer[i][1] = Data1[i][1];
      SendBuffer[i][2] = Data1[i][2];

      SendBuffer[i][3] = Data2[i][0];
      SendBuffer[i][4] = Data2[i][1];
      SendBuffer[i][5] = Data2[i][2];
   }


// send data to SendRank
   MPI_Isend( SendBuffer, 6*Nsend, GAMER_MPI_REAL, SendRank, 0, MPI_COMM_WORLD, &Req[0] );

// receive data from RecvRank
   MPI_Irecv( RecvBuffer, 6*Nrecv, GAMER_MPI_REAL, RecvRank, 0, MPI_COMM_WORLD, &Req[1] );

}



//---------------------------------------------------------------------------
// Function    :  Block_and_Copy__Force
// Description :  Return after all data have been properly received and copied
//                ~ for calculating acceleration and jerk
//---------------------------------------------------------------------------
void Block_and_Copy__Force( real Data1[][3], real Data2[][3], const real RecvBuffer[][6],
                            MPI_Request Req[], const int Nrecv )
{
   MPI_Waitall( 2, Req, MPI_STATUSES_IGNORE );

// copy data from RecvBuffer to Data1 and Data2 arrays
   for (int i=0; i<Nrecv; i++)
   {
      Data1[i][0] = RecvBuffer[i][0];
      Data1[i][1] = RecvBuffer[i][1];
      Data1[i][2] = RecvBuffer[i][2];

      Data2[i][0] = RecvBuffer[i][3];
      Data2[i][1] = RecvBuffer[i][4];
      Data2[i][2] = RecvBuffer[i][5];
   }
}



//---------------------------------------------------------------------------
// Function    :  DataTransfer_Pot
// Description :  Transfer data between different processes when calculating potential
//---------------------------------------------------------------------------
void DataTransfer_Pot( const real Other_Mass[], const real Other_Pos[][3],
                       real SendBuffer[][4], real RecvBuffer[][4], MPI_Request Req[] )
{

// copy data from "Other_XXX" arrays to SendBuffer
   for (int i=0; i<N; i++)
   {
      SendBuffer[i][0] = Other_Mass[i];

      SendBuffer[i][1] = Other_Pos[i][0];
      SendBuffer[i][2] = Other_Pos[i][1];
      SendBuffer[i][3] = Other_Pos[i][2];
   }


   const int BufferCount = 4*N;

// send data to SendRank rank
   MPI_Isend( SendBuffer, BufferCount, GAMER_MPI_REAL, SendRank, 0, MPI_COMM_WORLD, &Req[0] );

// receive data from RecvRank
   MPI_Irecv( RecvBuffer, BufferCount, GAMER_MPI_REAL, RecvRank, 0, MPI_COMM_WORLD, &Req[1] );

}



//---------------------------------------------------------------------------
// Function    :  Block_and_Copy__Pot
// Description :  Return after all data has been properly received and copied
//                ~ for calculating potential
//---------------------------------------------------------------------------
void Block_and_Copy__Pot( real Other_Mass[], real Other_Pos[][3], const real RecvBuffer[][4],
                          MPI_Request Req[] )
{
   MPI_Waitall( 2, Req, MPI_STATUSES_IGNORE );

// copy data from RecvBuffer to "Other_XXX" arrays
   for (int i=0; i<N; i++)
   {
      Other_Mass[i]   = RecvBuffer[i][0];

      Other_Pos[i][0] = RecvBuffer[i][1];
      Other_Pos[i][1] = RecvBuffer[i][2];
      Other_Pos[i][2] = RecvBuffer[i][3];
   }
}



//---------------------------------------------------------------------------
// Function    :  Get_Minimum_FromAllRanks
// Description :  Get the minimum value of the input variable from all ranks
//---------------------------------------------------------------------------
double Get_Minimum_FromAllRanks( double input )
{
   double min;

   MPI_Reduce( &input, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
   MPI_Bcast( &min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   return min;
}



//---------------------------------------------------------------------------
// Function    :  Get_L_FromAllRanks
// Description :  Get the number of line-up particles of different ranks
//---------------------------------------------------------------------------
void Get_L_FromAllRanks( int L_List[] )
{
   /*
   L_List[MyRank] = L;

   for (int Rank=0; Rank<NGPU; Rank++)    MPI_Bcast( &L_List[Rank], 1, MPI_INT, Rank, MPI_COMM_WORLD );
   */

   MPI_Allgather( &L, 1, MPI_INT, L_List, 1, MPI_INT, MPI_COMM_WORLD );
}



#ifdef HYBRID_SCHEME
//---------------------------------------------------------------------------
// Function    :  AllGather
// Description :  An alternative subroutine to MPI_Allgatherv in BeginHybrid_Acc_Jerk subroutine
//---------------------------------------------------------------------------
void AllGather( real Array_local[][3], real Array_total[][3], const int ArraySize[], const int Disp[] )
{

   int Rank_Send, Rank_Recv, NSend, NRecv;
   MPI_Request Req [2];

   NSend = ArraySize[MyRank];

   memcpy( Array_total[0]+Disp[MyRank], Array_local, 4*NSend );

   for (int i=1; i<NGPU; i++)
   {
      Rank_Send = (MyRank+i)%NGPU;
      Rank_Recv = (MyRank+NGPU-i)%NGPU;
      NRecv     = ArraySize[Rank_Recv];

      MPI_Isend( Array_local, NSend, GAMER_MPI_REAL, Rank_Send, 0, MPI_COMM_WORLD, &Req[0] );
      MPI_Irecv( Array_total[0]+Disp[Rank_Recv], NRecv, GAMER_MPI_REAL, Rank_Recv, 0, MPI_COMM_WORLD, &Req[1] );

      MPI_Waitall( 2, Req, MPI_STATUSES_IGNORE );
   }

}



//---------------------------------------------------------------------------
// Function    :  AllReduce
// Description :  An alternative subroutine to MPI_Reduce in BeginHybrid_Acc_Jerk subroutine
//---------------------------------------------------------------------------
void AllReduce( real Array_partial[][3], real Array_sum[][3], const int ArraySize[], const int Disp[] )
{

    int Rank_Send, Rank_Recv, NSend, NRecv;
    MPI_Request Req [2];
    Buffer = new real [L][3];

    NRecv  = 3*L;

    memcpy( Array_sum, Array_partial[0]+Disp[MyRank], 4*NRecv );

   for (int i=1; i<NGPU; i++)
   {
      Rank_Send = (MyRank+i)%NGPU;
      Rank_Recv = (MyRank+NGPU-i)%NGPU;
      NSend     = ArraySize[Rank_Send];

      MPI_Isend( Array_partial[0]+Disp[Rank_Send], NSend, GAMER_MPI_REAL, Rank_Send, 0, MPI_COMM_WORLD, &Req[0] );
      MPI_Irecv( Buffer, NRecv, GAMER_MPI_REAL, Rank_Recv, 0, MPI_COMM_WORLD, &Req[1] );

      MPI_Waitall( 2, Req, MPI_STATUSES_IGNORE );

      for (int j=0; j<L; j++)
      for (int dim=0; dim<3; dim++)    Array_sum[j][dim] += Buffer[j][dim];
   }

   delete [] Buffer;

}
#endif // HYBRID_SCHEME
