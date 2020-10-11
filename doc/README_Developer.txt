
Version 1.0
-----------------
1. only the serial version
2. only works with CUDA.1.0


Version 1.9
-----------------
1. parallel version
2. use the modified ring schene to try to improve the performance
   --> failed, the performance is lower than the original ring scheme


Version 2.0
-----------------
1. parallel version
2. two scheme for acceleration and jerk calculation
    a. RING scheme
    b. HYBRID scheme (Harfst, 2007)
    --> use OPTION9 in Makefile to choose the desired scheme


Version 2.1
-----------------
1. define the MIN_DT in Parameter file
2. modify the sequence of input parameters in "Parameter" file
3. use "DumpID" (0,1,2,...) rather than "Step" (0,153,9895,...) as the names of output files
4. declare all variables related to "Step" as "long int" instead of "int"


Version 2.1.1
-----------------
1. use double precision to calculate the initial time-step


Version 2.1.2
-----------------
1. slightly modify the procedure of time-step evaluation
   --> may find the next global time only among the particles within PlayerList
2. the "Parameter" and "Init_Condition" files should be placed in the same directory as executable file
3. terminate the program if the SOFTEN option is not turned on
4. input the order of the minimum and maximum time-step
5. input EPS rather than EPS_SQR


Version 2.2
-----------------
1. correct the evaluation of potentail energy
   --> subtract the error induced by self interaction
2. move all initial check to an isolate function
   --> CheckParameter
3. when restarting the program, do NOT output the initial condition
4. output error message if "Parameter" or "Init_Condition" files do not exit

*** BUG ***
==============================================
--> can not deal with Nj%NSplit != 0
--> it's a bug in CUCAL_Acc_Jerk_Split.1.4.cu
==============================================


Version 2.3
-----------------
1. update to CUCAL_Acc_Jerk_Split.2.0.cu
   --> fix the bug in version 1.4
2. mv the function "CheckParameter" to Auxiliary.cpp
3. enforce that NSplit must < Nj


Version 2.4
-----------------
1. modify the Makefile
   --> use "-Xopencc -O3" for the CU_CAL_SPLIT files when the N_IS_POWER_OF_TWO option is turned off
   --> much better performance


Version 2.5
-----------------
1. work with 64-bit and 32-bit machines
2. work with CUDA version 1.1 instead of version 0.8


Version 2.5.1
-----------------
1. work with CUDA SDK 2.1


Version 2.6    06/23/2010
-------------------------
1. Support the shared time-step scheme
2. Support the double precision for the accumulation arrays
3. All "TAB"s are removed
4. A brief user manual is put in the directory "Doc"


Version 2.6.1  09/20/2010
-------------------------
1. Support the Fermi GPUs.
2. Support the double precision in GPU.
3. Suppot different GPU selection modes.


Version 2.6.2  09/20/2010
-------------------------
1. Fix the bug in the function "Get_NewTimeStep", in which the time-step may
   be incorrect when the synchronization is performed and "missing particles"
   may happen.
2. Use the function "floor" for the time-step estimation.
3. The format of the file "Record__TimeStep" is slightly modified.
4. The format of the file "Record__Energy" is modified.
5. Replace all "MPI_INTEGER" by "MPI_INT".
6. Add headers in the file "Record__TimeStep" and "Record__Dump".
7. Allow the size of the input file to be larger than what is expected
   (TOTAL_N*7*4 bytes).


Version 2.7.0  09/11/2014
-------------------------
1. Support adding external acceleration from
   (1) analytical function
   (2) input file
   --> GRAVITY_TYPE, EXT_METHOD
2. For external acceleration, support calculating jerk as well
   --> EXT_JERK, EXT_ACC_DER
3. Support different particles interpolation schemes (CIC/TSC)
   for calculating the external acceleration at the particles'
   position
   --> EXT_PAR_INT
4. Support double precision for either only the accumulation
   arrays (FLOAT8_ACC) or all floating-point variables (FLOAT8)
5. Support OpenMP parallelization
   --> OPENMP, OMP_NTHREAD
6. Support initializing particles by an analytical function
   --> INIT_METHOD
7. Provide the GNUPLOT plotting tool
   --> tool/Plot__Template_Par.1.0.gpt
8. Set GRID_SIZE=56 in FERMI
9. Fix the bug in "Synchronize"
   --> Perform synchronization if **any** rank has line-up particles



BUG
-------------------------------------------------------------
1. The evaluation of potenial energy will give a wrong number if the
   soften length is extremly small
   (because the self potential will occupy almost all available digits)



Unfinished work
-------------------------------------------------------------
2. try to use "structure" to store the particle data --> speed up the data transfer between host and GPU  ?
3. cannot turn off the SOFTEN option
4. Add Descriptions for all functions
5. Set BLOCK_SIZE and GRID_SIZE in run time
   --> In this way, we can adjust them according to the adopted GPU.
