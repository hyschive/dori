CUCAL_Acc_Jerk program
-----------------------

Version 1.1
1. can NOT work with arbitrary number of particles

2. sum the acceleration and jerk block by block to reduce the accumulative error


Version 1.2
1. work with arbitrary number of particles

2. does NOT sum the acceleration and jerk block by block to reduce the accumulative error


Version 1.4
1. define the "N_IS_POWER_OF_TWO" in Makefile. Define it when N is power of 2 to remove the "if ( kk < Nj )"
   conditional determination and improve the performance


Version 1.4.1
1. real  --> real_acc
   float --> real




CUCAL_Pot program
-----------------------

Version 1.1
1. can NOT work with arbitrary number of particles

2. sum the potential block by block to reduce the accumulative error


Version 1.2
1. work with arbitrary number of particles

2. does NOT sum the potential block by block to reduce the accumulative error


Version 1.3
1. define the "N_IS_POWER_OF_TWO" in Makefile. Define it when N is power of 2 to remove the "if ( kk < Nj )"
   conditional determination and improve the performance


Version 1.3.1
1. real  --> real_acc
   float --> real



CUCAL_Acc_Jerk_Split program
-----------------------

Version 1.0
1. split the force calculation into several segments


Version 1.1
1. move some calculations to host computer


Version 1.2
1. replace some integral multiplication by "__umul24" function to improve the performance


Version 1.4
1. define the "N_IS_POWER_OF_TWO" in Makefile. Define it when N is power of 2 to remove the "if ( kk < Nj )"
   conditional determination and improve the performance


Version 2.0.1
1. real  --> real_acc
   float --> real

